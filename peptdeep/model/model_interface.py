import os
import numpy as np
import pandas as pd
import torch
import yaml
import inspect
from tqdm import tqdm

import torch.multiprocessing as mp
import functools

from types import ModuleType

# from torch.optim.lr_scheduler import LambdaLR
from transformers.optimization import get_cosine_schedule_with_warmup

from zipfile import ZipFile
from typing import IO, Tuple, List, Union
from alphabase.yaml_utils import save_yaml, load_yaml
from alphabase.peptide.precursor import is_precursor_sorted

from peptdeep.settings import model_const
from peptdeep.utils import ( 
    logging, process_bar, get_device,
    get_available_device
)
from peptdeep.settings import global_settings

from peptdeep.model.featurize import (
    get_ascii_indices, get_batch_aa_indices,
    get_batch_mod_feature
)

# def get_cosine_schedule_with_warmup(
#     optimizer, num_warmup_steps, 
#     num_training_steps, num_cycles=0.5, 
#     last_epoch=-1
# ):
#     """
#     Create a schedule with a learning rate that decreases following the
#     values of the cosine function between 0 and `pi * cycles` after a warmup
#     period during which it increases linearly between 0 and 1.
#     """
#     def lr_lambda(current_step):
#         if current_step < num_warmup_steps:
#             return float(current_step) / max(1, num_warmup_steps)
#         progress = float(
#             current_step - num_warmup_steps
#         ) / max(1, num_training_steps - num_warmup_steps)
#         return max(0.0, 0.5 * (
#             1.0 + np.cos(np.pi * num_cycles * 2.0 * progress)
#         ))

#     return LambdaLR(optimizer, lr_lambda, last_epoch)

def append_nAA_column_if_missing(precursor_df):
    """
    Append a column containing the number of Amino Acids
    """
    if 'nAA' not in precursor_df.columns:
        precursor_df['nAA'] = precursor_df.sequence.str.len()
        precursor_df.sort_values('nAA', inplace=True)
        precursor_df.reset_index(drop=True,inplace=True)
    return precursor_df

class ModelInterface(object):
    """
    Provides standardized methods to interact
    with ml models. Inherit into new class and override
    the abstract (i.e. not implemented) methods.
    """
    
    def __init__(self,
        device:str='gpu',
        fixed_sequence_len:int = 0,
        min_pred_value:float = 0.0,
        **kwargs
    ):
        """
        Parameters
        ----------
        device : str, optional
            device type in 'get_available', 'cpu', 'mps', 'gpu' (or 'cuda'),
            by default 'gpu'
        
        fixed_sequence_len : int, optional
            See :attr:`fixed_sequence_len`, defaults to 0.

        min_pred_value : float, optional
            See :attr:`min_pred_value`, defaults to 0.0.
        """
        self.model:torch.nn.Module = None
        self.optimizer = None
        self.model_params:dict = {}
        self.set_device(device)
        self.fixed_sequence_len = fixed_sequence_len
        self.min_pred_value = min_pred_value

    @property
    def fixed_sequence_len(self)->int:
        """
        This attribute controls how to train and infer for variable-length sequences:
        
        - if the value is 0, all sequence tensors will be grouped by nAA and train/infer on same nAA in batch.
        - if the value is > 0: all sequence tensors will be padded by zeros to the fixed length.
        - if the value is < 0: in each batch, padded by zeros to max length of the batch.
        """
        return self._fixed_sequence_len

    @fixed_sequence_len.setter
    def fixed_sequence_len(self, seq_len:int):
        self._fixed_sequence_len = seq_len
        self.model_params['fixed_sequence_len'] = seq_len

    @property
    def min_pred_value(self)->float:
        """
        The predicted values cannot be smaller than this value.
        """
        return self._min_pred_value

    @min_pred_value.setter
    def min_pred_value(self, val:float):
        self._min_pred_value = val
        self.model_params['min_pred_value'] = val

    @property
    def device_type(self)->str:
        """Read-only"""
        return self._device_type
    
    @property
    def device(self)->torch.device:
        """Read-only"""
        return self._device

    @property
    def device_ids(self)->list:
        """Read-only"""
        return self._device_ids

    @property
    def target_column_to_predict(self)->str:
        return self._target_column_to_predict

    @target_column_to_predict.setter
    def target_column_to_predict(self, column:str):
        self._target_column_to_predict = column

    @property
    def target_column_to_train(self)->str:
        return self._target_column_to_train

    @target_column_to_train.setter
    def target_column_to_train(self, column:str):
        self._target_column_to_train = column

    def set_device(self, 
        device_type:str = 'gpu', 
        device_ids:list = []
    ):
        """
        Set the device (e.g. gpu (cuda), mps, cpu, ...) to be used for the model.

        Parameters
        ----------
        device_type : str, optional
            Device type, see :data:`peptdeep.utils.torch_device_dict`.
            It will check available devices using 
            :meth:`peptdeep.utils.get_available_device()` 
            if device_type=='get_available'.
            By default 'gpu'

        device_ids : list, optional
            List of int. Device ids for cuda/gpu (e.g. [1,3] for cuda:1,3). 
            By default []
        """
        self._device_ids = device_ids

        if device_type == 'get_available':
            self._device, self._device_type = get_available_device()
        else:
            self._device, self._device_type = get_device(
                device_type, device_ids
            )

        self._model_to_device()

    def _model_to_device(self):
        """ 
        Enable multiple GPUs using torch.nn.DataParallel.

        TODO It is better to use torch.nn.parallel.DistributedDataParallel, 
        but this may need more setups for models and optimizers.
        """
        if self.model is None: return
        if self.device_type != 'cuda':
            self.model.to(self.device)
        else:
            if (
                self.device_ids and len(self.device_ids) > 1
            ):
                self.model = torch.nn.DataParallel(self.model, self.device_ids)
            elif (
                not self.device_ids and torch.cuda.device_count()>1
            ):
                self.model = torch.nn.DataParallel(self.model)
            self.model.to(self.device)

    def build(self,
        model_class: torch.nn.Module,
        **kwargs
    ):
        """
        Builds the model by specifying the PyTorch module, 
        the parameters, the device, the loss function ...
        """
        self.model = model_class(**kwargs)
        self.model_params.update(**kwargs)
        self._model_to_device()
        self._init_for_training()

    def train_with_warmup(self,
        precursor_df: pd.DataFrame,
        *,
        batch_size=1024, 
        epoch=10, 
        warmup_epoch=5,
        lr=1e-4,
        verbose=False,
        verbose_each_epoch=False,
        **kwargs
    ):
        """
        Train the model according to specifications. Includes a warumup 
        phase with linear increasing and cosine decreasing for lr scheduling).
        """
        self._prepare_training(precursor_df, lr, **kwargs)

        lr_scheduler = self._get_lr_schedule_with_warmup(
            warmup_epoch, epoch
        )

        for epoch in range(epoch):
            if self.fixed_sequence_len == 0:
                batch_cost = self._train_one_epoch(
                    precursor_df, epoch,
                    batch_size, verbose_each_epoch,
                    **kwargs
                )
            else:
                batch_cost = self._train_one_epoch_by_padding_zeros(
                    precursor_df, epoch,
                    batch_size, verbose_each_epoch,
                    **kwargs
                )

            lr_scheduler.step()
            if verbose: print(
                f'[Training] Epoch={epoch+1}, lr={lr_scheduler.get_last_lr()[0]}, loss={np.mean(batch_cost)}'
            )
        
        torch.cuda.empty_cache()

    def train(self,
        precursor_df: pd.DataFrame,
        *,
        batch_size=1024, 
        epoch=10, 
        warmup_epoch:int=0,
        lr=1e-4,
        verbose=False,
        verbose_each_epoch=False,
        **kwargs
    ):
        """
        Train the model according to specifications.
        """

        if verbose: logging.info(
            f"Training with fixed sequence length: {self.fixed_sequence_len}"
        )

        if warmup_epoch > 0:
            self.train_with_warmup(
                precursor_df,
                batch_size=batch_size,
                epoch=epoch,
                warmup_epoch=warmup_epoch,
                lr=lr,
                verbose=verbose,
                verbose_each_epoch=verbose_each_epoch,
                **kwargs
            )
        else:
            self._prepare_training(precursor_df, lr, **kwargs)

            for epoch in range(epoch):
                if self.fixed_sequence_len == 0:
                    batch_cost = self._train_one_epoch(
                        precursor_df, epoch,
                        batch_size, verbose_each_epoch,
                        **kwargs
                    )
                else:
                    batch_cost = self._train_one_epoch_by_padding_zeros(
                        precursor_df, epoch,
                        batch_size, verbose_each_epoch,
                        **kwargs
                    )
                if verbose: print(f'[Training] Epoch={epoch+1}, Mean Loss={np.mean(batch_cost)}')
            
            torch.cuda.empty_cache()

    def predict(self,
        precursor_df:pd.DataFrame,
        *,
        batch_size:int=1024,
        verbose:bool=False,
        **kwargs
    )->pd.DataFrame:
        """
        The model predicts the properties based on the inputs it has been trained for.
        Returns the ouput as a pandas dataframe.
        """
        precursor_df = append_nAA_column_if_missing(precursor_df)
        self._pad_zeros_if_fixed_len(precursor_df)
        self._check_predict_in_order(precursor_df)
        self._prepare_predict_data_df(precursor_df,**kwargs)
        self.model.eval()

        _grouped = precursor_df.groupby('nAA')
        if verbose:
            batch_tqdm = tqdm(_grouped)
        else:
            batch_tqdm = _grouped
        with torch.no_grad():
            for nAA, df_group in batch_tqdm:
                for i in range(0, len(df_group), batch_size):
                    batch_end = i+batch_size
                    
                    batch_df = df_group.iloc[i:batch_end,:]

                    features = self._get_features_from_batch_df(
                        batch_df, **kwargs
                    )

                    if isinstance(features, tuple):
                        predicts = self._predict_one_batch(*features)
                    else:
                        predicts = self._predict_one_batch(features)

                    self._set_batch_predict_data(
                        batch_df, predicts, 
                        **kwargs
                    )

        torch.cuda.empty_cache()
        return self.predict_df

    def predict_mp(self,
        precursor_df:pd.DataFrame,
        *,
        batch_size:int=1024,
        mp_batch_size:int=100000,
        process_num:int=global_settings['thread_num'],
        **kwargs
    )->pd.DataFrame:
        """
        Predicting with multiprocessing is no GPUs are availible.
        Note this multiprocessing method only works for models those predict
        values within (inplace of) the precursor_df.
        """
        precursor_df = append_nAA_column_if_missing(precursor_df)

        if self.device_type != 'cpu':
            return self.predict(
                precursor_df, 
                batch_size=batch_size,
                verbose=False,
                **kwargs
            )
            
        _predict_func = functools.partial(self.predict, 
            batch_size=batch_size, verbose=False, **kwargs
        )

        def batch_df_gen(precursor_df, mp_batch_size):
            for i in range(0, len(precursor_df), mp_batch_size):
                yield precursor_df.iloc[i:i+mp_batch_size]

        self._check_predict_in_order(precursor_df)
        self._prepare_predict_data_df(precursor_df,**kwargs)

        print("Predicting with multiprocessing ...")
        self.model.share_memory()
        df_list = []
        with mp.get_context('spawn').Pool(process_num) as p:
            for ret_df in process_bar(p.imap(
                    _predict_func,
                    batch_df_gen(precursor_df, mp_batch_size),
                ), len(precursor_df)//mp_batch_size+1
            ):
                df_list.append(ret_df)

        self.predict_df = pd.concat(df_list)
        self.predict_df.reset_index(drop=True, inplace=True)
        
        return self.predict_df

    def save(self, filename:str):
        """
        Save the model state, the constants used, the code defining the model and the model parameters.
        """
        # TODO save tf.keras.Model
        dir = os.path.dirname(filename)
        if not dir: dir = './'
        if not os.path.exists(dir): os.makedirs(dir)
        torch.save(self.model.state_dict(), filename)
        with open(filename+'.txt','w') as f: f.write(str(self.model))
        save_yaml(filename+'.model_const.yaml', model_const)
        self._save_codes(filename+'.model.py')
        save_yaml(filename+'.param.yaml', self.model_params)

    def load(
        self,
        model_file: Tuple[str, IO],
        model_path_in_zip: str = None,
        **kwargs
    ):
        """
        Load a model specified in a zip file, a text file or a file stream.
        """
        # TODO load tf.keras.Model
        if isinstance(model_file, str):
            # We may release all models (msms, rt, ccs, ...) in a single zip file
            if model_file.lower().endswith('.zip'):
                self._load_model_from_zipfile(model_file, model_path_in_zip)
            else:
                self._load_model_from_pytorchfile(model_file)
        else:
            self._load_model_from_stream(model_file)

    def get_parameter_num(self):
        """
        Get total number of parameters in model.
        """
        return np.sum([p.numel() for p in self.model.parameters()])

    def build_from_py_codes(self,
        model_code_file_or_zip:str,
        code_file_in_zip:str=None,
        include_model_params_yaml:bool=True,
        **kwargs
    ):
        """
        Build the model based on a python file. Must contain a PyTorch 
        model implemented as 'class Model(...'
        """
        if model_code_file_or_zip.lower().endswith('.zip'):
            with ZipFile(model_code_file_or_zip, 'r') as model_zip:
                with model_zip.open(code_file_in_zip,'r') as f:
                    codes = f.read()
                if include_model_params_yaml:
                    with model_zip.open(
                        code_file_in_zip[:-len('model.py')]+'param.yaml',
                        'r'
                    ) as f:
                        params = yaml.load(f, yaml.FullLoader)
        else:
            with open(model_code_file_or_zip, 'r') as f:
                codes = f.read()
            if include_model_params_yaml:
                params = load_yaml(
                    model_code_file_or_zip[:-len('model.py')]+'param.yaml'
                )

        compiled_codes = compile(
            codes, 
            filename='model_file_py',
            mode='exec'
        )
        _module = ModuleType('_apd_nn_codes')
        #codes must contains torch model codes 'class Model(...'
        exec(compiled_codes, _module.__dict__)

        if include_model_params_yaml:
            for key, val in params.items():
                if key not in kwargs:
                    kwargs[key] = val

        self.model = _module.Model(**kwargs)
        self.model_params = kwargs
        self.model.to(self.device)
        self._init_for_training()

    def _init_for_training(self):
        """
        Set the loss function, and more attributes for different tasks.
        The default loss function is nn.L1Loss.
        """
        self.loss_func = torch.nn.L1Loss()

    def _as_tensor(self, 
        data:np.ndarray, 
        dtype:torch.dtype=torch.float32
    )->torch.Tensor:
        """Convert numerical np.array to pytorch tensor.
        The tensor will be stored in self.device

        Parameters
        ----------
        data : np.ndarray
            Numerical np.ndarray to be converted as a tensor
            
        dtype : torch.dtype, optional
            The dtype of the indices used for embedding should be `torch.long`. 
            Defaults to `torch.float32`

        Returns
        -------
        torch.Tensor
            The tensor stored in self.device
        """
        return torch.tensor(data, dtype=dtype, device=self.device)

    def _load_model_from_zipfile(self, model_file, model_path_in_zip):
        with ZipFile(model_file) as model_zip:
            with model_zip.open(model_path_in_zip,'r') as pt_file:
                self._load_model_from_stream(pt_file)

    def _load_model_from_pytorchfile(self, model_file):
        with open(model_file,'rb') as pt_file:
            self._load_model_from_stream(pt_file)

    def _load_model_from_stream(self, stream):
        (
            missing_keys, unexpect_keys 
        ) = self.model.load_state_dict(torch.load(
            stream, map_location=self.device),
            strict=False
        )
        if len(missing_keys) > 0:
            logging.warn(f"nn parameters {missing_keys} are MISSING while loading models in {self.__class__}")
        if len(unexpect_keys) > 0:
            logging.warn(f"nn parameters {unexpect_keys} are UNEXPECTED while loading models in {self.__class__}")

    def _save_codes(self, save_as):
        try:
            code = '''import torch\n'''
            code += '''import peptdeep.model.building_block as building_block\n'''
            code += '''from peptdeep.model.model_shop import *\n'''
            class_code = inspect.getsource(self.model.__class__)
            code += 'class Model' + class_code[class_code.find('('):]
            with open(save_as, 'w') as f:
                f.write(code)
        except (TypeError, ValueError, KeyError) as e:
            logging.info(f'Cannot save model source codes: {str(e)}')

    def _train_one_epoch_by_padding_zeros(self, 
        precursor_df, epoch, batch_size, verbose_each_epoch, 
        **kwargs
    ):
        """Training for an epoch by padding zeros"""
        batch_cost = []
        rnd_df = precursor_df.sample(frac=1)
        if verbose_each_epoch:
            batch_tqdm = tqdm(range(0, len(rnd_df), batch_size))
        else:
            batch_tqdm = range(0, len(rnd_df), batch_size)
        for i in batch_tqdm:
            batch_end = i+batch_size

            batch_df = rnd_df.iloc[i:batch_end,:]
            targets = self._get_targets_from_batch_df(
                batch_df, **kwargs
            )
            features = self._get_features_from_batch_df(
                batch_df, **kwargs
            )
            if isinstance(features, tuple):
                batch_cost.append(
                    self._train_one_batch(targets, *features)
                )
            else:
                batch_cost.append(
                    self._train_one_batch(targets, features)
                )
                
        if verbose_each_epoch:
            batch_tqdm.set_description(
                f'Epoch={epoch+1}, batch={len(batch_cost)}, loss={batch_cost[-1]:.4f}'
            )
        return batch_cost

    def _train_one_epoch(self, 
        precursor_df, epoch, batch_size, verbose_each_epoch, 
        **kwargs
    ):
        """Training for an epoch"""
        batch_cost = []
        _grouped = list(precursor_df.sample(frac=1).groupby('nAA'))
        rnd_nAA = np.random.permutation(len(_grouped))
        if verbose_each_epoch:
            batch_tqdm = tqdm(rnd_nAA)
        else:
            batch_tqdm = rnd_nAA
        for i_group in batch_tqdm:
            nAA, df_group = _grouped[i_group]
            # df_group = df_group.reset_index(drop=True)
            for i in range(0, len(df_group), batch_size):
                batch_end = i+batch_size

                batch_df = df_group.iloc[i:batch_end,:]
                targets = self._get_targets_from_batch_df(
                    batch_df, **kwargs
                )
                features = self._get_features_from_batch_df(
                    batch_df, **kwargs
                )
                if isinstance(features, tuple):
                    batch_cost.append(
                        self._train_one_batch(targets, *features)
                    )
                else:
                    batch_cost.append(
                        self._train_one_batch(targets, features)
                    )
                
            if verbose_each_epoch:
                batch_tqdm.set_description(
                    f'Epoch={epoch+1}, nAA={nAA}, batch={len(batch_cost)}, loss={batch_cost[-1]:.4f}'
                )
        return batch_cost

    def _train_one_batch(
        self, 
        targets:torch.Tensor, 
        *features,
    ):
        """Training for a mini batch"""
        self.optimizer.zero_grad()
        predicts = self.model(*features)
        cost = self.loss_func(predicts, targets)
        cost.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), 1.0)
        self.optimizer.step()
        return cost.item()

    def _predict_one_batch(self,
        *features
    ):
        """Predicting for a mini batch"""
        return self.model(
            *features
        ).cpu().detach().numpy()

    def _get_targets_from_batch_df(self,
        batch_df:pd.DataFrame, **kwargs,
    )->torch.Tensor:
        """Tell the `train()` method how to get target values from the `batch_df`.
           All sub-classes must re-implement this method.
           Use torch.tensor(np.array, dtype=..., device=self.device) to convert tensor.

        Parameters
        ----------
        batch_df : pd.DataFrame
            Dataframe of each mini batch.

        Returns
        -------
        torch.Tensor
            Target value tensor
        """
        return self._as_tensor(
            batch_df[self.target_column_to_train].values, 
            dtype=torch.float32
        )

    def _get_aa_indice_features_padding_zeros(
        self, batch_df:pd.DataFrame
    )->torch.LongTensor:
        """
        Get indices values of variable length sequences 
        using 128 ascii codes
        """
        if self.fixed_sequence_len < 0:
            max_len = batch_df.nAA.max()
        else:
            max_len = self.fixed_sequence_len
        return self._as_tensor(
            get_ascii_indices(
                batch_df['sequence'].apply(
                    lambda seq: seq + chr(0)*(max_len-len(seq))
                ).values.astype('U')
            ), 
            dtype=torch.long
        )

    def _get_aa_indice_features(
        self, batch_df:pd.DataFrame
    )->torch.LongTensor:
        """
        Get indices values for fixed length sequences 
        with 128 ascii codes.
        """
        return self._as_tensor(
            get_ascii_indices(
                batch_df['sequence'].values.astype('U')
            ), 
            dtype=torch.long
        )

    def _get_26aa_indice_features(
        self, batch_df:pd.DataFrame
    )->torch.LongTensor:
        """
        Get indices values for 26 upper-case letters (amino acids), 
        from 1 to 26. 0 is used for padding.
        """
        return self._as_tensor(
            get_batch_aa_indices(
                batch_df['sequence'].values.astype('U')
            ), 
            dtype=torch.long
        )
        
    def _get_features_from_batch_df(self,
        batch_df:pd.DataFrame, **kwargs,
    )->Union[torch.LongTensor, Tuple[torch.Tensor]]:
        """
        Any sub-class must re-implement this method:
        
        - Return `self._get_aa_features()` for sequence-level prediciton 
        - Return `self._get_aa_mod_features()` for modified sequence-level

        Parameters
        ----------
        batch_df : pd.DataFrame
            Batch of precursor dataframe.

        Returns
        -------
        Union[torch.LongTensor, Tuple[torch.Tensor]]: 
            A LongTensor if the sub-class call `self._get_aa_features(batch_df)` (default).
            Or a tuple of tensors if call `self._get_aa_mod_features(batch_df)`.
        """
        return self._get_aa_features(batch_df)

    def _get_aa_mod_features(self,
        batch_df:pd.DataFrame, **kwargs,
    )->Tuple[torch.Tensor]:
        return (
            self._get_aa_features(batch_df),
            self._get_mod_features(batch_df)
        )

    def _get_mod_features(
        self, batch_df:pd.DataFrame
    )->torch.Tensor:
        """
        Get modification features.
        """
        if self.fixed_sequence_len < 0:
            batch_df = batch_df.copy()
            batch_df['nAA'] = batch_df.nAA.max()
        return self._as_tensor(
            get_batch_mod_feature(batch_df)
        )

    def _get_aa_features(self, 
        batch_df:pd.DataFrame
    )->torch.LongTensor:
        """
        Get AA indices
        """
        if self.fixed_sequence_len == 0:
            return self._get_aa_indice_features(batch_df)
        else:
            return self._get_aa_indice_features_padding_zeros(batch_df)

    def _prepare_predict_data_df(self,
        precursor_df:pd.DataFrame, 
        **kwargs
    ):
        """
        This methods fills 0s in the column of 
        `self.target_column_to_predict` in `precursor_df`,
        and then does `self.predict_df=precursor_df`.
        """
        precursor_df[self.target_column_to_predict] = 0.0
        self.predict_df = precursor_df

    def _prepare_train_data_df(self,
        precursor_df:pd.DataFrame, 
        **kwargs
    ):
        """Changes to the training dataframe can be implemented here.

        Parameters
        ----------
        precursor_df : pd.DataFrame
            Dataframe containing the training data.
        """
        pass

    def _set_batch_predict_data(self,
        batch_df:pd.DataFrame,
        predict_values:np.ndarray,
        **kwargs
    ):
        """Set predicted values into `self.predict_df`.

        Parameters
        ----------
        batch_df : pd.DataFrame
            Dataframe of mini batch when predicting

        predict_values : np.array
            Predicted values
        """
        predict_values[predict_values<self._min_pred_value] = self._min_pred_value
        if self._predict_in_order:
            self.predict_df.loc[:,self.target_column_to_predict].values[
                batch_df.index.values[0]:batch_df.index.values[-1]+1
            ] = predict_values
        else:
            self.predict_df.loc[
                batch_df.index,self.target_column_to_predict
            ] = predict_values

    def _set_optimizer(self, lr):
        """Set optimizer"""
        self.optimizer = torch.optim.Adam(
            self.model.parameters(), lr=lr
        )

    def set_lr(self, lr:float):
        """Set learning rate"""
        if self.optimizer is None:
            self._set_optimizer(lr)
        else:
            for g in self.optimizer.param_groups:
                g['lr'] = lr

    def _get_lr_schedule_with_warmup(self, warmup_epoch, epoch):
        if warmup_epoch > epoch:
            warmup_epoch = epoch//2
        return get_cosine_schedule_with_warmup(
            self.optimizer, warmup_epoch, epoch
        )

    def _pad_zeros_if_fixed_len(self, precursor_df:pd.DataFrame):
        if self.fixed_sequence_len > 0:
            precursor_df.drop(
                index=precursor_df[
                    precursor_df.nAA>self.fixed_sequence_len
                ].index, 
                inplace=True,
            )
            precursor_df.reset_index(drop=True, inplace=True)
            precursor_df['nAA'] = self.fixed_sequence_len

    def _prepare_training(self, 
        precursor_df:pd.DataFrame, 
        lr:float, 
        **kwargs
    ):
        if 'nAA' not in precursor_df.columns:
            precursor_df['nAA'] = precursor_df.sequence.str.len()
        self._pad_zeros_if_fixed_len(precursor_df)
        self._prepare_train_data_df(precursor_df, **kwargs)
        self.model.train()

        self.set_lr(lr)

    def _check_predict_in_order(self, precursor_df:pd.DataFrame):
        if is_precursor_sorted(precursor_df):
            self._predict_in_order = True
        else:
            self._predict_in_order = False
