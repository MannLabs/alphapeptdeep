import torch
from typing import Optional, Tuple

class CustomTransformerEncoderLayer(torch.nn.Module):
    def __init__(self, config):
        super().__init__()
        self.config = config
        self.self_attn = BertAttention(config)
        self.intermediate = BertIntermediate(config)
        self.output = BertOutput(config)
        self.attn_output_weights = None
    def forward(
        self,
        hidden_states: torch.Tensor,
        src_mask: Optional[torch.FloatTensor] = None,
        is_causal: Optional[bool] = None,
        src_key_padding_mask: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor]:

        self_attention_outputs = self.self_attn(
            hidden_states=hidden_states,
            attention_mask=src_mask,
            output_attentions=self.config.output_attentions
        )

        attention_output, outputs = self_attention_outputs[0], self_attention_outputs[1:]
        self.attn_output_weights = outputs if self.config.output_attentions else None

        layer_output = self.output(self.intermediate(attention_output), attention_output)

        return layer_output

class BertAttention(torch.nn.Module):
    def __init__(self, config):
        super().__init__()
        assert config.batch_first == True, "`batch_first` must be   `True`."
        self.batch_first = config.batch_first

        self.self = torch.nn.MultiheadAttention(
            embed_dim=config.hidden_size,
            num_heads=config.num_attention_heads,
            dropout=config.hidden_dropout_prob,
            batch_first=self.batch_first,
        )
        self.dense = torch.nn.Linear(config.hidden_size, config.hidden_size)
        self.dropout = torch.nn.Dropout(config.hidden_dropout_prob)
        self.LayerNorm = torch.nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)

    def forward(
        self,
        hidden_states: torch.Tensor,
        attention_mask: Optional[torch.FloatTensor] = None,
        output_attentions: bool = False,
    ) -> Tuple[torch.Tensor]:

        self_outputs = self.self(
            query=hidden_states,
            key=hidden_states,
            value=hidden_states,
            need_weights=output_attentions,
            attn_mask=attention_mask,
        )

        attention_output = self.LayerNorm(
            self.dropout(self.dense(self_outputs[0]))
            + hidden_states
        )
        return (attention_output,) + self_outputs[1:]

class BertIntermediate(torch.nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dense = torch.nn.Linear(config.hidden_size, config.intermediate_size)

    def forward(self, hidden_states: torch.Tensor) -> torch.Tensor:
        return torch.nn.functional.gelu(self.dense(hidden_states))

class BertOutput(torch.nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dense = torch.nn.Linear(config.intermediate_size, config.hidden_size)
        self.LayerNorm = torch.nn.LayerNorm(config.hidden_size, eps=config.layer_norm_eps)
        self.dropout = torch.nn.Dropout(config.hidden_dropout_prob)

    def forward(self, hidden_states: torch.Tensor, input_tensor: torch.Tensor) -> torch.Tensor:
        hidden_states = self.dense(hidden_states)
        hidden_states = self.dropout(hidden_states)
        hidden_states = self.LayerNorm(hidden_states + input_tensor)
        return hidden_states

class CustomTransformerEncoder(torch.nn.TransformerEncoder):
    def __init__(self, config):
        customtransformerencoderlayer = CustomTransformerEncoderLayer(config)

        super().__init__(
            encoder_layer=customtransformerencoderlayer,
            num_layers=config.num_hidden_layers,
            norm=None,
            enable_nested_tensor=False,
            mask_check=False,            
        )
    def forward(
        self,
        src: torch.Tensor,
        attention_mask: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        return super().forward(src, mask=attention_mask)
