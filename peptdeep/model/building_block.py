import torch
import numpy as np
#BERT from huggingface
from transformers.models.bert.modeling_bert import BertEncoder

from peptdeep.settings import model_const
from peptdeep.settings import global_settings as settings

torch.set_num_threads(2)

mod_feature_size = len(model_const['mod_elements'])
max_instrument_num = model_const['max_instrument_num']
frag_types = settings['model']['frag_types']
max_frag_charge = settings['model']['max_frag_charge']
num_ion_types = len(frag_types)*max_frag_charge
aa_embedding_size = model_const['aa_embedding_size']

def aa_embedding(hidden_size):
    return torch.nn.Embedding(aa_embedding_size, hidden_size, padding_idx=0)

def ascii_embedding(hidden_size):
    return torch.nn.Embedding(128, hidden_size, padding_idx=0)

def aa_one_hot(aa_indices, *cat_others):
    aa_x = torch.nn.functional.one_hot(
        aa_indices, aa_embedding_size
    )
    return torch.cat((aa_x, *cat_others), 2)

def instrument_embedding(hidden_size):
    return torch.nn.Embedding(max_instrument_num, hidden_size)

def zero_param(*shape):
    return torch.nn.Parameter(torch.zeros(shape), requires_grad=False)

def xavier_param(*shape):
    x = torch.nn.Parameter(torch.empty(shape), requires_grad=False)
    torch.nn.init.xavier_uniform_(x)
    return x

def invert_attention_mask(
    encoder_attention_mask:torch.Tensor,
    dtype=torch.float32
)->torch.FloatTensor:
    """
    See `invert_attention_mask()` in https://github.com/huggingface/transformers/blob/main/src/transformers/modeling_utils.py#L737.
    Invert an attention mask (e.g., switches 0. and 1.).
    Args:
        encoder_attention_mask (`torch.Tensor`): An attention mask.
    Returns:
        `torch.Tensor`: The inverted attention mask.
    """
    if encoder_attention_mask.dim() == 3:
        encoder_extended_attention_mask = encoder_attention_mask[:, None, :, :]
    if encoder_attention_mask.dim() == 2:
        encoder_extended_attention_mask = encoder_attention_mask[:, None, None, :]
    encoder_extended_attention_mask = encoder_extended_attention_mask.to(dtype=dtype)  # fp16 compatibility
    encoder_extended_attention_mask = (1.0 - encoder_extended_attention_mask) * torch.finfo(dtype).min

    return encoder_extended_attention_mask


init_state = xavier_param

class SeqCNN_MultiKernel(torch.nn.Module):
    """
    Extract sequence features using `torch.nn.Conv1D` with 
    different kernel sizes (1(residue connection),3,5,7), 
    and then concatenate the outputs of these Conv1Ds.
    """
    def __init__(self, out_features:int):
        """
        Parameters
        ----------
        out_features : int
            Must be divided by 4.

        Raises
        ------
        ValueError
            "out_features must be divided by 4"
            
        """
        super().__init__()

        hidden = out_features//4
        if hidden*4 != out_features:
            raise ValueError('out_features must be divided by 4')

        self.cnn_short = torch.nn.Conv1d(
            hidden, hidden,
            kernel_size=3, padding=1
        )
        self.cnn_medium = torch.nn.Conv1d(
            hidden, hidden,
            kernel_size=5, padding=2
        )
        self.cnn_long = torch.nn.Conv1d(
            hidden, hidden,
            kernel_size=7, padding=3
        )

    def forward(self, x):
        x = x.transpose(1, 2)
        x1 = self.cnn_short(x)
        x2 = self.cnn_medium(x)
        x3 = self.cnn_long(x)
        return torch.cat((x, x1, x2, x3), dim=1).transpose(1,2)

#legacy
class SeqCNN(torch.nn.Module):
    """
    Extract sequence features using `torch.nn.Conv1D` with 
    different kernel sizes (1(residue connection),3,5,7), and then concatenate 
    the outputs of these Conv1Ds. The Output dim is 4*embedding_hidden.
    """
    def __init__(self, embedding_hidden):
        super().__init__()

        self.cnn_short = torch.nn.Conv1d(
            embedding_hidden, embedding_hidden,
            kernel_size=3, padding=1
        )
        self.cnn_medium = torch.nn.Conv1d(
            embedding_hidden, embedding_hidden,
            kernel_size=5, padding=2
        )
        self.cnn_long = torch.nn.Conv1d(
            embedding_hidden, embedding_hidden,
            kernel_size=7, padding=3
        )

    def forward(self, x):
        x = x.transpose(1, 2)
        x1 = self.cnn_short(x)
        x2 = self.cnn_medium(x)
        x3 = self.cnn_long(x)
        return torch.cat((x, x1, x2, x3), dim=1).transpose(1,2)


class Seq_Transformer(torch.nn.Module):
    """
    Using PyTorch built-in Transformer layers
    """
    def __init__(self,
        in_features,
        hidden_features,
        nheads=8,
        nlayers=2,
        dropout=0.1
    ):
        super().__init__()
        encoder_layers = torch.nn.TransformerEncoderLayer(
            in_features, nheads, hidden_features, dropout
        )
        self.transformer_encoder = torch.nn.TransformerEncoder(
            encoder_layers, nlayers
        )
        
    def forward(self, x):
        return self.transformer_encoder(x.permute(1,0,2)).permute(1,0,2)


class Hidden_Transformer(torch.nn.Module):
    """
    Transformer NN based on pytorch's built-in TransformerLayer class
    """
    def __init__(self, 
        hidden, hidden_expand=4,
        nheads=8, nlayers=4, dropout=0.1
    ):
        super().__init__()
        self.transormer = Seq_Transformer(
            hidden, hidden*hidden_expand, nheads=nheads, 
            nlayers=nlayers, dropout=dropout
        )
    def forward(self, x):
        return self.transormer(x)

class _Pseudo_Bert_Config:
    def __init__(self, 
        hidden_dim=256, 
        intermediate_size=1024,
        num_attention_heads=8,
        num_bert_layers=4,
        dropout=0.1,
        output_attentions=False,
    ):
        self.add_cross_attention = False
        self.chunk_size_feed_forward = 0
        self.is_decoder = False
        self.seq_len_dim = 1
        self.training = False
        self.hidden_act = "gelu"
        self.hidden_dropout_prob = dropout
        self.attention_probs_dropout_prob = dropout
        self.hidden_size = hidden_dim
        self.initializer_range = 0.02
        self.intermediate_size = intermediate_size
        self.layer_norm_eps = 1e-8
        self.num_attention_heads = num_attention_heads
        self.num_hidden_layers = num_bert_layers
        self.output_attentions = output_attentions

class Hidden_HFace_Transformer(torch.nn.Module):
    """
    Transformer NN based on HuggingFace's BertEncoder class
    """
    def __init__(self, 
        hidden_dim, hidden_expand=4,
        nheads=8, nlayers=4, dropout=0.1,
        output_attentions=False,
    ):
        super().__init__()
        self.config = _Pseudo_Bert_Config(
            hidden_dim=hidden_dim,
            intermediate_size=hidden_dim*hidden_expand,
            num_attention_heads=nheads,
            num_bert_layers=nlayers,
            dropout=dropout,
            output_attentions=False
        )
        self.output_attentions = output_attentions
        self.bert = BertEncoder(self.config)
    def forward(self, x:torch.Tensor, 
        attention_mask:torch.Tensor=None,
    )->tuple:
        """
        Parameters
        ----------
        x : torch.Tensor
            shape = (batch, seq_len, dim)
        attention_mask : torch.Tensor
            shape = (batch, seq_len), [0,1] tensor , 1=enable

        Returns
        -------
        (Tensor, [Tensor])
            out[0] is the hidden layer output, 
            and out[1] is the output attention 
            if self.output_attentions==True
        """
        if attention_mask is not None:
            attention_mask = invert_attention_mask(
                attention_mask, dtype=x.dtype
            )
        return self.bert(
            x,
            attention_mask=attention_mask,
            output_attentions=self.output_attentions,
            return_dict=False
        )
#legacy
HiddenBert = Hidden_HFace_Transformer

class HFace_Transformer_with_PositionalEncoder(torch.nn.Module):
    """
    HuggingFace transformer with a positional encoder in front.

    Parameters
    ----------
    hidden_dim : int
        Input and output feature dimension.

    hidden_expand : int, optional
        FFN hidden size = hidden*hidden_expand. Defaults to 4.

    nhead : int, optional
        Multi-head attention number. Defaults to 8.

    nlayers : int, optional
        Number of transformer layers. Defaults to 4.

    dropout : float, optional
        Dropout rate. Defaults to 0.1.

    output_attentions : bool, optional
        If output attention values. Defaults to False.

    max_len : int, optional
        Max input sequence length. Defaults to 200.
    """
    def __init__(self,
        hidden_dim:int, hidden_expand=4,
        nheads=8, nlayers=4, dropout=0.1,
        output_attentions=False,
        max_len=200,
    ):
        super().__init__()
        self.pos_encoder = PositionalEncoding(hidden_dim, max_len=max_len)
        self.bert = Hidden_HFace_Transformer(
            hidden_dim=hidden_dim, hidden_expand=hidden_expand,
            nheads=nheads, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )
    def forward(self, 
        x:torch.Tensor,
        attention_mask:torch.Tensor=None,
    )->tuple:
        """
        Parameters
        ----------
        x : torch.Tensor
            Input tensor

        Returns
        -------
        tuple
            Tensor: Output tensor.
            [Tensor]: Attention tensor, returned only if output_attentions is True.
        """
        x = self.pos_encoder(x)
        return self.bert(x, attention_mask)

class SeqLSTM(torch.nn.Module):
    """
    returns LSTM applied on sequence input
    """
    def __init__(self, in_features, out_features, 
                 rnn_layer=2, bidirectional=True
        ):
        super().__init__()

        if bidirectional:
            if out_features%2 != 0:
                raise ValueError("'out_features' must be able to be divided by 2")
            hidden = out_features//2
        else:
            hidden = out_features

        self.rnn_h0 = init_state(
            rnn_layer+rnn_layer*bidirectional,
            1, hidden
        )
        self.rnn_c0 = init_state(
            rnn_layer+rnn_layer*bidirectional,
            1, hidden
        ) 
        self.rnn = torch.nn.LSTM(
            input_size = in_features,
            hidden_size = hidden,
            num_layers = rnn_layer,
            batch_first = True,
            bidirectional = bidirectional,
        )

    def forward(self, x:torch.Tensor):
        h0 = self.rnn_h0.repeat(1, x.size(0), 1)
        c0 = self.rnn_c0.repeat(1, x.size(0), 1)
        x, _ = self.rnn(x, (h0,c0))
        return x

class SeqGRU(torch.nn.Module):
    """
    returns GRU applied on sequence input
    """
    def __init__(self, in_features, out_features, 
                 rnn_layer=2, bidirectional=True
        ):
        super().__init__()

        if bidirectional:
            if out_features%2 != 0:
                raise ValueError("'out_features' must be able to be divided by 2")
            # to make sure that output dim is out_features
            # as `bidirectional` will cat forward and reverse RNNs
            hidden = out_features//2
        else:
            hidden = out_features
        
        self.rnn_h0 = init_state(
            rnn_layer+rnn_layer*bidirectional, 
            1, hidden
        )
        self.rnn = torch.nn.GRU(
            input_size = in_features,
            hidden_size = hidden,
            num_layers = rnn_layer,
            batch_first = True,
            bidirectional = bidirectional,
        )

    def forward(self, x:torch.Tensor):
        h0 = self.rnn_h0.repeat(1, x.size(0), 1)
        x, _ = self.rnn(x, h0)
        return x

class SeqAttentionSum(torch.nn.Module):
    """
    apply linear transformation and tensor rescaling with softmax
    """
    def __init__(self, in_features):
        super().__init__()
        self.attn = torch.nn.Sequential(
            torch.nn.Linear(in_features, 1, bias=False),
            torch.nn.Softmax(dim=1),
        )
    
    def forward(self, x):
        attn = self.attn(x)
        return torch.sum(torch.mul(x, attn), dim=1)

class PositionalEncoding(torch.nn.Module):
    """
    transform sequence input into a positional representation
    """
    def __init__(self, out_features=128, max_len = 200):
        super().__init__()

        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(
                0, out_features, 2
            ) * (-np.log(max_len) / out_features)
        )
        pe = torch.zeros(1, max_len, out_features)
        pe[0, :, 0::2] = torch.sin(position * div_term)
        pe[0, :, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x):
        return x + self.pe[:,:x.size(1),:]

class PositionalEmbedding(torch.nn.Module):
    """
    transform sequence with the standard embedding function
    """
    def __init__(self, out_features=128, max_len=200):
        super().__init__()

        self.pos_emb = torch.nn.Embedding(
            max_len, out_features
        )

    def forward(self, x:torch.Tensor):
        return x + self.pos_emb(torch.arange(
            x.size(1), dtype=torch.long, device=x.device
        ).unsqueeze(0))

class Meta_Embedding(torch.nn.Module):
    """Encodes Charge state, Normalized Collision Energy (NCE) and Instrument for a given spectrum 
    into a 'meta' single layer network
    """
    def __init__(self,
        out_features,
    ):
        super().__init__()
        self.nn = torch.nn.Linear(
            max_instrument_num+1, out_features-1
        )

    def forward(self,
        charges, NCEs, instrument_indices,
    ):
        inst_x = torch.nn.functional.one_hot(
            instrument_indices, max_instrument_num
        )
        meta_x = self.nn(torch.cat((inst_x, NCEs), 1))
        meta_x = torch.cat((meta_x, charges), 1)
        return meta_x
#legacy
InputMetaNet = Meta_Embedding

class Mod_Embedding_FixFirstK(torch.nn.Module):
    """
    Encodes the modification vector in a single layer feed forward network, but not transforming the first k features
    """
    def __init__(self,
        out_features,
    ):
        super().__init__()
        self.k = 6
        self.nn = torch.nn.Linear(
            mod_feature_size-self.k, out_features-self.k,
            bias=False
        )

    def forward(self,
        mod_x,
    ):
        return torch.cat((
            mod_x[:,:,:self.k], 
            self.nn(mod_x[:,:,self.k:])
        ), 2)
#legacy
InputModNetFixFirstK = Mod_Embedding_FixFirstK

class AA_Mod_Embedding(torch.nn.Module):
    """
    Concatenates the AA (128 ASCII codes) embedding with the modifcation vector
    """
    def __init__(self,
        out_features,
        mod_feature_size = 8,
    ):
        super().__init__()
        self.mod_embedding = Mod_Embedding_FixFirstK(
            mod_feature_size
        )
        self.aa_embedding = ascii_embedding(
            out_features-mod_feature_size
        )
    def forward(self, aa_indices, mod_x):
        mod_x = self.mod_embedding(mod_x)
        aa_x = self.aa_embedding(aa_indices)
        return torch.cat((aa_x, mod_x), 2)
#legacy
InputAAEmbedding = AA_Mod_Embedding

class Mod_Embedding(torch.nn.Module):
    """
    Encodes the modification vector in a single layer feed forward network
    """
    def __init__(self,
        out_features,
    ):
        super().__init__()
        self.nn = torch.nn.Linear(
            mod_feature_size, out_features,
            bias=False
        )

    def forward(self,
        mod_x,
    ):
        return self.nn(mod_x)
#legacy
InputModNet = Mod_Embedding

class Input_26AA_Mod_PositionalEncoding(torch.nn.Module):
    """
    Encodes AA (26 AA letters) and modification vector
    """
    def __init__(self, out_features, max_len=200):
        super().__init__()
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        self.aa_emb = aa_embedding(
            out_features-mod_hidden
        )
        self.pos_encoder = PositionalEncoding(
            out_features, max_len
        )
        
    def forward(self, 
        aa_indices, mod_x
    ):
        mod_x = self.mod_nn(mod_x)
        x = self.aa_emb(aa_indices)
        return self.pos_encoder(torch.cat((x, mod_x), 2))
#legacy
AATransformerEncoding = Input_26AA_Mod_PositionalEncoding

class Input_AA_Mod_PositionalEncoding(torch.nn.Module):
    """
    Encodes AA (ASCII codes) and modification vector
    """
    def __init__(self, out_features, max_len=200):
        super().__init__()
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        self.aa_emb = ascii_embedding(
            out_features-mod_hidden
        )
        self.pos_encoder = PositionalEncoding(
            out_features, max_len
        )
        
    def forward(self, 
        aa_indices, mod_x
    ):
        mod_x = self.mod_nn(mod_x)
        x = self.aa_emb(aa_indices)
        return self.pos_encoder(torch.cat((x, mod_x), 2))

class Input_AA_Mod_Charge_PositionalEncoding(torch.nn.Module):
    """
    Embed AA (128 ASCII codes), modification, and charge state
    """
    def __init__(self, out_features, max_len=200):
        super().__init__()
        mod_hidden = 8
        self.charge_dim = 2
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        self.aa_emb = ascii_embedding(
            out_features-mod_hidden-self.charge_dim
        )
        self.pos_encoder = PositionalEncoding(
            out_features, max_len
        )
        
    def forward(self, 
        aa_indices, mod_x, charges
    ):
        mod_x = self.mod_nn(mod_x)
        x = self.aa_emb(aa_indices)
        charge_x = charges.unsqueeze(1).repeat(
            1, mod_x.size(1), self.charge_dim
        )
        return self.pos_encoder(torch.cat((x, mod_x,charge_x), 2))

class Input_26AA_Mod_LSTM(torch.nn.Module):
    """
    Applies an LSTM network to a AA (26 AA letters) sequence & modifications
    """
    def __init__(self,
        out_features,
        n_lstm_layers=1,
    ):
        super().__init__()
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        self.lstm = SeqLSTM(
            aa_embedding_size+mod_hidden,
            out_features,
            n_lstm_layers=n_lstm_layers, 
            bidirectional=True
        )
    def forward(self, aa_indices, mod_x):
        mod_x = self.mod_nn(mod_x)
        x = aa_one_hot(aa_indices, mod_x)
        return self.lstm(x)
#legacy
InputAALSTM = Input_26AA_Mod_LSTM

        
class Input_26AA_Mod_Meta_LSTM(torch.nn.Module):
    """
    Applies a LSTM network to a AA (26 AA letters) sequence and modifications,
    and concatenates with 'meta' information (charge, nce, instrument_indices) 
    """
    def __init__(self,
        out_features,
    ):
        super().__init__()
        meta_dim = 4
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        self.meta_nn = Meta_Embedding(meta_dim)
        self.nn = SeqLSTM(
            aa_embedding_size+mod_hidden,
            out_features-meta_dim,
            rnn_layer=1, bidirectional=True
        )
        
    def forward(self, 
        aa_indices, mod_x, charges, NCEs, instrument_indices
    ):
        mod_x = self.mod_nn(mod_x)
        x = aa_one_hot(aa_indices, mod_x)
        x = self.nn(x)
        meta_x = self.meta_nn(
            charges, NCEs, instrument_indices
        ).unsqueeze(1).repeat(1, mod_x.size(1), 1)
        return torch.cat((x, meta_x), 2)
#legacy
InputAALSTM_cat_Meta = Input_26AA_Mod_Meta_LSTM


class Input_26AA_Mod_Charge_LSTM(torch.nn.Module):
    """
    Applies a LSTM network to a AA (26 AA letters) sequence and modifications, 
    and concatenates with charge state information
    """
    def __init__(self,
        out_features,
    ):
        super().__init__()
        self.charge_dim = 2
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        self.nn = SeqLSTM(
            aa_embedding_size+mod_hidden,
            out_features-self.charge_dim,
            rnn_layer=1, bidirectional=True
        )
        
    def forward(self, aa_indices, mod_x, charges):
        mod_x = self.mod_nn(mod_x)
        x = aa_one_hot(aa_indices, mod_x)
        x = self.nn(x)
        charge_x = charges.unsqueeze(1).repeat(
            1, mod_x.size(1), self.charge_dim
        )
        return torch.cat((x, charge_x), 2)
#legacy
InputAALSTM_cat_Charge = Input_26AA_Mod_Charge_LSTM


class Seq_Meta_LSTM(torch.nn.Module):
    """
    Takes a hidden layer which processes the hidden tensor 
    as well as the 'meta' information of NCE, Instrument, Charge
    """
    def __init__(self,
        in_features,
        out_features,
    ):
        super().__init__()
        meta_dim = 4
        self.meta_nn = Meta_Embedding(meta_dim)
        self.nn = SeqLSTM(
            in_features+meta_dim,
            out_features,
            rnn_layer=1, bidirectional=False
        )
        
    def forward(self, x, charges, NCEs, instrument_indices):
        meta_x = self.meta_nn(
            charges, NCEs, instrument_indices
        ).unsqueeze(1).repeat(1, x.size(1), 1)
        return self.nn(torch.cat((x, meta_x), 2))
#legacy
OutputLSTM_cat_Meta = Seq_Meta_LSTM
    
class Seq_Meta_Linear(torch.nn.Module):
    """
    takes a hidden linear which processed the 'meta' information of NCE, Instrument, Charge
    """
    def __init__(self,
        in_features,
        out_features,
    ):
        super().__init__()
        meta_dim = 4
        self.meta_nn = Meta_Embedding(meta_dim)
        self.nn = torch.nn.Linear(
            in_features+meta_dim,
            out_features,
            bias=False
        )
        
    def forward(self, x, charges, NCEs, instrument_indices):
        meta_x = self.meta_nn(
            charges, NCEs, instrument_indices
        ).unsqueeze(1).repeat(1, x.size(1), 1)
        return self.nn(torch.cat((x, meta_x), 2))
#legacy
OutputLinear_cat_Meta = Seq_Meta_Linear

class Encoder_26AA_Mod_LSTM(torch.nn.Module):
    """
    Two LSTM layers on AA (26 AA letters) and modifications.
    """
    def __init__(self, out_features, n_lstm_layers=1):
        super().__init__()
        
        self.input_nn = Input_26AA_Mod_LSTM(out_features)
        self.nn = SeqLSTM(
            out_features, out_features, 
            rnn_layer=n_lstm_layers
        )

    def forward(self, aa_indices, mod_x):
        x = self.input_nn(aa_indices, mod_x)
        x = self.nn(x)
        return x

#legacy
Input_AA_LSTM_Encoder = Encoder_26AA_Mod_LSTM


class Encoder_26AA_Mod_CNN_LSTM(torch.nn.Module):
    """
    Encode AAs (26 AA letters) and modifications by CNN and LSTM layers
    """
    def __init__(self, out_features, n_lstm_layers=1):
        super().__init__()

        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)
        input_dim = aa_embedding_size+mod_hidden
        self.input_cnn = SeqCNN(input_dim)
        self.hidden_nn = SeqLSTM(
            input_dim*4, out_features, 
            rnn_layer=n_lstm_layers
        ) #SeqCNN outputs 4*input_dim

    def forward(self, aa_indices, mod_x):
        mod_x = self.mod_nn(mod_x)
        x = aa_one_hot(aa_indices, mod_x)
        x = self.input_cnn(x)
        x = self.hidden_nn(x)
        return x

#legacy
Input_AA_CNN_Encoder = Encoder_26AA_Mod_CNN_LSTM

class Encoder_26AA_Mod_CNN_LSTM_AttnSum(torch.nn.Module):
    """
    Encode AAs (26 AA letters) and modifications by CNN and LSTM layers, 
    then by 'SeqAttentionSum'.
    """
    def __init__(self, out_features, n_lstm_layers=2):
        super().__init__()
        
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)

        input_dim = aa_embedding_size+mod_hidden
        self.input_cnn = SeqCNN(input_dim)
        self.hidden_nn = SeqLSTM(
            input_dim*4, out_features, 
            rnn_layer=n_lstm_layers
        ) #SeqCNN outputs 4*input_dim
        self.attn_sum = SeqAttentionSum(out_features)

    def forward(self, aa_indices, mod_x):
        mod_x = self.mod_nn(mod_x)
        x = aa_one_hot(aa_indices, mod_x)
        x = self.input_cnn(x)
        x = self.hidden_nn(x)
        x = self.attn_sum(x)
        return x
#legacy
Input_AA_CNN_LSTM_Encoder = Encoder_26AA_Mod_CNN_LSTM_AttnSum

class Encoder_AA_Mod_CNN_LSTM_AttnSum(torch.nn.Module):
    """
    Encode AAs (128 ASCII codes) and modifications by CNN and LSTM layers, 
    and then by 'SeqAttentionSum'.
    """
    def __init__(self, out_features, n_lstm_layers=2):
        super().__init__()
        
        mod_hidden = 8
        input_dim = out_features//4
        self.aa_mod_embedding = AA_Mod_Embedding(
            input_dim, mod_feature_size=mod_hidden
        )
        self.input_cnn = SeqCNN(input_dim)
        self.hidden_nn = SeqLSTM(
            input_dim*4, out_features, 
            rnn_layer=n_lstm_layers
        ) #SeqCNN outputs 4*input_dim
        self.attn_sum = SeqAttentionSum(out_features)

    def forward(self, aa_indices, mod_x):
        x = self.aa_mod_embedding(aa_indices, mod_x)
        x = self.input_cnn(x)
        x = self.hidden_nn(x)
        x = self.attn_sum(x)
        return x


class Encoder_AA_Mod_Transformer(torch.nn.Module):
    """
    AAs (128 ASCII codes) and modifications embedded by Bert, 
    then encoded by 'SeqAttentionSum'.
    """
    def __init__(self,out_features,
        dropout=0.1,
        nlayers=4,
        output_attentions=False
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)
        
        self.input_nn = Input_AA_Mod_PositionalEncoding(out_features)

        self.output_attentions = output_attentions
        self.encoder = Hidden_HFace_Transformer(
            out_features, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )
    def forward(self, aa_indices, mod_x,
        attention_mask=None
    ):
        x = self.input_nn(
            aa_indices, mod_x
        )
        x = self.dropout(x)

        x = self.encoder(x, attention_mask)
        if self.output_attentions:
            self.attentions = x[1]
        else:
            self.attentions = None
        return x[0]

class Encoder_AA_Mod_Transformer_AttnSum(torch.nn.Module):
    """
    Encode AAs (128 ASCII codes) and modifications by transformers.
    """
    def __init__(self,out_features,
        dropout=0.1,
        nlayers=4,
        output_attentions=False
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)
        
        self.encoder_nn = Encoder_AA_Mod_Transformer(
            out_features, dropout=dropout, nlayers=nlayers,
            output_attentions=output_attentions 
        )
        self.attn_sum = SeqAttentionSum(out_features)
        
    def forward(self, aa_indices, mod_x):
        x = self.encoder_nn(aa_indices, mod_x)
        return self.dropout(self.attn_sum(x))

class Encoder_AA_Mod_Charge_Transformer(torch.nn.Module):
    """
    Encode AAs (128 ASCII codes), modifications and charge by transformers.
    """
    def __init__(self,out_features,
        dropout=0.1,
        nlayers=4,
        output_attentions=False
    ):
        super().__init__()
        
        self.dropout = torch.nn.Dropout(dropout)
        self.input_nn = Input_AA_Mod_Charge_PositionalEncoding(out_features)

        self.output_attentions = output_attentions
        self.encoder = Hidden_HFace_Transformer(
            out_features, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )
    def forward(self, aa_indices, mod_x, charges,
        attention_mask=None,
    ):
        x = self.input_nn(
            aa_indices, mod_x, charges
        )
        x = self.dropout(x)

        x = self.encoder(x, attention_mask)
        if self.output_attentions:
            self.attentions = x[1]
        else:
            self.attentions = None
        return x[0]

class Encoder_AA_Mod_Charge_Transformer_AttnSum(torch.nn.Module):
    """
    Encode AAs (128 ASCII codes), modifications and charge by transformers, 
    and then by 'SeqAttentionSum'
    """
    def __init__(self,out_features,
        dropout=0.1,
        nlayers=4,
        output_attentions=False
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)
        
        self.encoder_nn = Encoder_AA_Mod_Charge_Transformer(
            out_features, dropout=dropout, nlayers=nlayers,
            output_attentions=output_attentions 
        )
        self.attn_sum = SeqAttentionSum(out_features)
    def forward(self, aa_indices, mod_x, charges):
        x = self.encoder_nn(aa_indices, mod_x, charges)
        return self.dropout(self.attn_sum(x))

class Encoder_26AA_Mod_Charge_CNN_LSTM_AttnSum(torch.nn.Module):
    """
    Encode AAs (26 AA letters), modifications and charge by transformers, 
    and then by 'SeqAttentionSum'
    """
    def __init__(self, out_features):
        super().__init__()
        
        mod_hidden = 8
        self.mod_nn = Mod_Embedding_FixFirstK(mod_hidden)

        input_dim = aa_embedding_size+mod_hidden+1
        self.input_cnn = SeqCNN(input_dim)
        self.hidden_nn = SeqLSTM(
            input_dim*4, out_features, rnn_layer=2
        ) #SeqCNN outputs 4*input_dim
        self.attn_sum = SeqAttentionSum(out_features)

    def forward(self, aa_indices, mod_x, charges):
        mod_x = self.mod_nn(mod_x)
        x = aa_one_hot(
            aa_indices, mod_x, 
            charges.unsqueeze(1).repeat(1,mod_x.size(1),1)
        )
        x = self.input_cnn(x)
        x = self.hidden_nn(x)
        x = self.attn_sum(x)
        return x

#legacy
Input_AA_CNN_LSTM_cat_Charge_Encoder = Encoder_26AA_Mod_Charge_CNN_LSTM_AttnSum


class Decoder_LSTM(torch.nn.Module):
    """
    Decode with LSTM
    """
    def __init__(self, in_features, out_features):
        super().__init__()

        hidden = 128
        self.rnn = SeqLSTM(
            in_features, out_features,
            rnn_layer=1, bidirectional=False,
        )

        self.output_nn = torch.nn.Linear(
            hidden, out_features, bias=False
        )
    
    def forward(self, x:torch.tensor, output_len):
        x = self.rnn(
            x.unsqueeze(1).repeat(1,output_len,1)
        )
        x = self.output_nn(x)
        return x
#legacy
SeqLSTMDecoder = Decoder_LSTM

class Decoder_GRU(torch.nn.Module):
    """
    Decode with GRU
    """
    def __init__(self, in_features, out_features):
        super().__init__()

        hidden = 128
        self.rnn = SeqGRU(
            in_features, out_features,
            rnn_layer=1, bidirectional=False,
        )

        self.output_nn = torch.nn.Linear(
            hidden, out_features, bias=False
        )
    
    def forward(self, x:torch.tensor, output_len):
        x = self.rnn(
            x.unsqueeze(1).repeat(1,output_len,1)
        )
        x = self.output_nn(x)
        return x
#legacy
SeqGRUDecoder = Decoder_GRU

class Decoder_Linear(torch.nn.Module):
    """
    Decode w linear NN
    """
    def __init__(self, in_features, out_features):
        super().__init__()

        self.nn = torch.nn.Sequential(
            torch.nn.Linear(in_features, 64),
            torch.nn.PReLU(),
            torch.nn.Linear(64, out_features),
        )

    def forward(self, x):
        return self.nn(x)
#legacy
LinearDecoder = Decoder_Linear
