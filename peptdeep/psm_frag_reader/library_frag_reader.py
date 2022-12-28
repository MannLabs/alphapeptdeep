import pandas as pd
import numpy as np

from alphabase.io.psm_reader.dia_search_reader import (
    SpectronautReader
)

import alphabase.peptide.mobility as mobility

from peptdeep.psm_frag_reader.psm_frag_reader import (
    PSMReader_w_FragBase,
    psm_w_frag_reader_provider
)

class SpectronautMSMSReader(SpectronautReader, PSMReader_w_FragBase):
    def __init__(self,
        frag_types=['b','y','b_modloss','y_modloss'], 
        max_frag_charge=2,
        rt_unit='irt',
        **kwargs
    ):
        PSMReader_w_FragBase.__init__(self,
            frag_types = frag_types,
            max_frag_charge = max_frag_charge,
            **kwargs
        )

        SpectronautReader.__init__(self, rt_unit=rt_unit)

        self.min_allow_frag_num = 6
        self.groupby_raw_name = False
    
    @property
    def fragment_intensity_df(self):
        return self._fragment_intensity_df

    def _find_mapped_columns(self, lib_df):
        self.seq_col = None
        for col in self.column_mapping['sequence']:
            if col in lib_df.columns:
                self.seq_col = col
                break
        self.rt_col = None
        for col in self.column_mapping['rt']:
            if col in lib_df.columns:
                self.rt_col = col
                break
        self.mob_col = None
        for col in self.column_mapping['mobility']:
            if col in lib_df.columns:
                self.mob_col = col
                break
        self.raw_col = None
        if self.groupby_raw_name:
            if isinstance(self.column_mapping['raw_name'],str):
                if self.column_mapping['raw_name'] in lib_df.columns:
                    self.raw_col = self.column_mapping['raw_name']
            else:
                for col in self.column_mapping['raw_name']:
                    if col in lib_df.columns:
                        self.raw_col = col
                        break

    def _get_fragment_intensity(self, lib_df):

        frag_col_dict = dict(zip(
            self.charged_frag_types, 
            range(len(self.charged_frag_types))
        ))

        self._find_mapped_columns(lib_df)

        mod_seq_list = []
        seq_list = []
        charge_list = []
        rt_list = []
        mob_list = []
        frag_intens_list = []
        nAA_list = []
        raw_list = []

        group_cols = [self.mod_seq_column, self.seq_col, 'PrecursorCharge']

        if self.raw_col is not None:
            group_cols.append(self.raw_col)
        
        for keys, df_group in lib_df.groupby(
            group_cols
        ):
            if len(df_group) < self.min_allow_frag_num: continue
            if self.raw_col is None:
                mod_seq, seq, charge = keys
            else:
                mod_seq, seq, charge, raw = keys
            nAA = len(seq)
            intens = np.zeros(
                (nAA-1, len(self.charged_frag_types)),dtype=np.float32
            )
            for frag_type, frag_num, loss_type, frag_charge, inten in df_group[
                [
                    'FragmentType','FragmentNumber','FragmentLossType',
                    'FragmentCharge','RelativeIntensity'
                ]
            ].values:
                if frag_type in 'abc':
                    frag_num -= 1
                elif frag_type in 'xyz':
                    frag_num = nAA-frag_num-1
                else:
                    continue
                
                if loss_type == 'noloss':
                    frag_type = f'{frag_type}_z{frag_charge}'
                elif loss_type == 'H3PO4':
                    frag_type = f'{frag_type}_modloss_z{frag_charge}'
                else:
                    continue
                
                if frag_type not in frag_col_dict:
                    continue
                frag_col_idx = frag_col_dict[frag_type]
                intens[frag_num, frag_col_idx] = inten
            max_inten = np.max(intens)
            if max_inten <= 0: continue
            intens /= max_inten

            mod_seq_list.append(mod_seq)
            seq_list.append(seq)
            charge_list.append(charge)
            rt_list.append(df_group[self.rt_col].values[0])
            mob_list.append(df_group[self.mob_col].values[0])
            frag_intens_list.append(intens)
            nAA_list.append(nAA)
            if self.raw_col is not None:
                raw_list.append(raw)
        
        df = pd.DataFrame({
            self.mod_seq_column: mod_seq_list,
            self.seq_col: seq_list,
            'PrecursorCharge': charge_list,
            self.rt_col: rt_list,
            self.mob_col: mob_list,
        })

        if self.raw_col is not None:
            df[self.raw_col] = raw_list

        self._fragment_intensity_df = pd.DataFrame(
            np.concatenate(frag_intens_list),
            columns = self.charged_frag_types
        )

        indices = np.zeros(len(nAA_list)+1, dtype=np.int64)
        indices[1:] = np.array(nAA_list)-1
        indices = np.cumsum(indices)

        df['frag_start_idx'] = indices[:-1]
        df['frag_stop_idx'] = indices[1:]

        return df

    def _load_file(self, filename):
        df = pd.read_csv(filename, sep=self.csv_sep)
        self._find_mod_seq_column(df)

        df = self._get_fragment_intensity(df)

        return df

    def _post_process(self, 
        lib_df
    ):  
        self._psm_df['nAA'] = self._psm_df.sequence.str.len()
        self._psm_df[
            ['frag_start_idx','frag_stop_idx']
        ] = lib_df[['frag_start_idx','frag_stop_idx']]

        self.normalize_rt_by_raw_name()

        if (
            'mobility' in self._psm_df.columns
        ):
            self._psm_df['ccs'] = (
                mobility.mobility_to_ccs_for_df(
                    self._psm_df,
                    'mobility'
                )
            )
        
        self._psm_df = self._psm_df[
            ~self._psm_df.mods.isna()
        ].reset_index(drop=True)



psm_w_frag_reader_provider.register_reader('spectronaut', SpectronautMSMSReader)
