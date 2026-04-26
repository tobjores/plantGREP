import torch
from torch import nn
from torch import Tensor
from torch.nn.common_types import _size_1_t
from torch.nn.modules.utils import _single
from typing import Type, Union, List

# Helper function to turn single elements to a list of the desired size
def _to_list(element, length, accept_tuple = True) -> List:
    if accept_tuple and isinstance(element, tuple):
        element = list(element)
        element_type = 'tuple'
    else:
        element_type = 'list'
    
    element = element if isinstance(element, list) else [element]
    if len(element) == 1:
        element = element * length
    assert len(element) == length, f"Received {element_type} of wrong length. Length should be 1 or {length}, but is {len(element)}"
    return element

# Aggregation functions for the BiConv1d layer
def cat_tensors(*tensors):
    return torch.cat((tensors), dim = 1)

aggregation_fns = {
    'max': torch.max,
    'sum': torch.add,
    'cat': cat_tensors
}

# Define a bidirectional convolutional layer
class BiConv1d(nn.Conv1d):
    """
    Bidirectional 1D convolution for one-hot encoded DNA sequences
    The sequences are scanned twice, once with unmodified kernels and once with reverse-complemented kernels

    Args:
        kernels (int, optional): Number of convolving kernels. Default: 128
        kernel_size (int or tuple, optional): Size of the convolving kernel. Default: 7
        aggregation (str, optional): Method to aggregate the output from the forward and reverse-complemented kernels. 'max', 'sum', or 'cat'. Default: 'max'
    """

    def __init__(
        self,
        kernels: int = 128,
        kernel_size: _size_1_t = 7,
        aggregation: str = 'max',
        **kwargs
    ) -> None:
        kernel_size_ = _single(kernel_size)
        super().__init__(
            in_channels = 4,
            out_channels = kernels,
            kernel_size = kernel_size_,
            **kwargs
        )
        self.kernels = kernels
        assert aggregation in aggregation_fns, f"Invalid value for `aggregation`: {aggregation!r}; should be one of {list(aggregation_fns)}"
        self.aggregation = aggregation
        self.aggr_fn = aggregation_fns[self.aggregation]
        self.kwargs = {**kwargs}
        self.out_channels = 2 * self.kernels if self.aggregation == 'cat' else self.kernels

    def extra_repr(self):
        s = super().extra_repr() + ', aggregation={aggregation}'
        return s.format(**self.__dict__)

    def forward(self, x: Tensor) -> Tensor:
        out_fwd = self._conv_forward(x, self.weight, self.bias)
        out_rev = self._conv_forward(x, torch.flip(self.weight, dims = [1, 2]), self.bias)
        out = self.aggr_fn(out_fwd, out_rev)
        return out

# Define a convolutional module with a combination of a regular and a bidirectional convolutional layer
class semiBiConv(nn.Module):
    """
    Applies a regular and a bidirectional convolutional layer to the input and concatenates the results

    Args:
        ss_kernels (int, optional): Number of kernels in the regulator convolutional layer. Default: 32
        ds_kernels (int, optional): Number of kernels in the bidirectional convolutional layer. Default: 128
        kernel_size (int or list, optional): Size of the convolving kernels in the regular and bidirectional convolutional layers. Default: 7
        act_fn (torch.nn.Module or list, optional): Activation function used after the regular and bidirectional convolutional layer. Default: torch.nn.ReLU
        dropout_rate (float, list, optional): Probability for the dropout layer after the regular and bidirectional convolutional layer. Default: 0.2
        aggregation (str, optional): Method to aggregate the output from the forward and reverse-complemented kernels in the bidirectional convolutional layer. 'max', 'sum', or 'cat'. Default: 'max'
        padding (int, tuple or str, optional): Padding added to both sides of the input. Default: 'same'
    """

    def __init__(
        self,
        ss_kernels: int = 32,
        ds_kernels: int = 128,
        kernel_size: Union[int, List[int]] = 7,
        act_fn: Union[Type[nn.Module], List[Type[nn.Module]]] = nn.ReLU,
        dropout_rate: Union[float, List[float]] = 0.2,
        aggregation: str = 'max',
        padding: Union[str, _size_1_t] = 'same',
        **kwargs
    ) -> None:
        super().__init__()
        self.ss_kernels = ss_kernels if ss_kernels is not None else 0
        self.ds_kernels = ds_kernels if ds_kernels is not None else 0
        self.kernel_size = _to_list(kernel_size, 2)
        self.act_fn = _to_list(act_fn, 2)
        self.dropout_rate = _to_list(dropout_rate, 2)
        self.aggregation = aggregation
        self.padding = padding
        self.kwargs = {**kwargs}
        if self.ss_kernels == 0 and self.ds_kernels == 0:
            raise ValueError("At least one of `ss_kernels` and `ds_kernels` must be larger than 0")
        if self.padding != 'same' and self.kernel_size[0] != self.kernel_size[1]:
            raise ValueError("`padding` must be 'same' when using different kernel sizes for the regular and bidirectional convolutional layers")
        self.in_channels = 4
        self.out_channels = ss_kernels + 2 * self.ds_kernels if self.aggregation == 'cat' else self.ss_kernels + self.ds_kernels
        if self.ss_kernels > 0:
            self.ss_layer = nn.Conv1d(in_channels = 4, out_channels = self.ss_kernels, kernel_size = self.kernel_size[0], padding = self.padding, **self.kwargs)
            ss_net = [
                self.ss_layer,
                self.act_fn[0]()
            ]
            if self.dropout_rate[0] != 0.0 and self.dropout_rate[0] is not None:
                ss_net.append(nn.Dropout(p = self.dropout_rate[0]))
            self.ss_net = nn.Sequential(*ss_net)
        if self.ds_kernels > 0:
            self.ds_layer = BiConv1d(kernels = self.ds_kernels, kernel_size = self.kernel_size[1], aggregation = self.aggregation, padding = self.padding, **self.kwargs)
            ds_net = [
                self.ds_layer,
                self.act_fn[1]()
            ]
            if self.dropout_rate[1] != 0.0 and self.dropout_rate[1] is not None:
                ds_net.append(nn.Dropout(p = self.dropout_rate[1]))
            self.ds_net = nn.Sequential(*ds_net)

    def forward(self, x: Tensor) -> Tensor:
        if self.ss_kernels > 0 and self.ds_kernels > 0:
            out_ss = self.ss_net(x)
            out_ds = self.ds_net(x)
            out = torch.cat((out_ss, out_ds), dim = 1)
        elif self.ss_kernels > 0:
            out = self.ss_net(x)
        elif self.ds_kernels > 0:
            out = self.ds_net(x)
        return out
    
    def get_out_length(self, length_in: int) -> int:
        first_layer = self.ss_layer if self.ss_kernels > 0 else self.ds_layer
        length_out = (length_in + sum(first_layer._reversed_padding_repeated_twice) - int(first_layer.dilation[0]) * (first_layer.kernel_size[0] - 1) - 1) // first_layer.stride[0] + 1
        return length_out


# Define a dense layer (applies two succesive convolutional layers to the input and concatenates the results with the input):
class DenseLayer(nn.Module):
    """
    Applies two succesive convolutional layers and one dropout layer to the input and concatenates the results with the input

    Args:
        in_channels (int): Number of input channels
        growth_rate (int, optional): Number of final ouput channels. Default: 12
        R (int, optional): Multiplier of `growth_rate` to determine the number of ouput channels in the first convolutional layer. Default: 2
        kernel_size (int or list, optional): Size of the convolving kernels. Default: [1, 3]
        act_fn (torch.nn.Module, optional): Activation function used after each convolutional layer. Default: torch.nn.ReLU
        dropout_rate (float, optional): Probability for the dropout layer. Default: 0.2
        batchnorm (bool, optional): Include a BatchNorm1d layer. Default: False
    """
    def __init__(
        self,
        in_channels: int,
        growth_rate: int = 12,
        R: int = 2,
        kernel_size: Union[int, List[int]] = [1, 3],
        act_fn: Type[nn.Module] = nn.ReLU,
        dropout_rate: float = 0.2,
        batchnorm: bool = False
    ) -> None:
        super().__init__()
        self.in_channels = in_channels
        self.growth_rate = growth_rate
        self.R = R
        self.kernel_size = _to_list(kernel_size, 2)
        self.act_fn = act_fn
        self.dropout_rate = dropout_rate
        self.batchnorm = batchnorm
        self.int_out = self.growth_rate * self.R
        self.out_channels = growth_rate
        net = [
            nn.Conv1d(in_channels = self.in_channels, out_channels = self.int_out, kernel_size = self.kernel_size[0], stride = 1, padding = 'same'),
            act_fn(),
            nn.Conv1d(in_channels = self.int_out, out_channels = self.out_channels, kernel_size = self.kernel_size[1], stride = 1, padding = 'same'),
            act_fn()
        ]
        if self.batchnorm:
            net = [nn.BatchNorm1d(self.in_channels)] + net
        if self.dropout_rate != 0.0 and self.dropout_rate is not None:
            net.append(nn.Dropout(p = self.dropout_rate))
        self.net = nn.Sequential(*net)

    def forward(self, x: Tensor) -> Tensor:
        out = self.net(x)
        out = torch.cat((x, out), dim = 1)
        return out

# Define a dense block (combines multiple dense layers):
class DenseBlock(nn.Module):
    """
    Combines multiple dense layers

    Args:
        in_channels (int): Number of input channels
        layers (int, optional): Number of dense layers. Default: 6
        growth_rate (int, optional): Number of final ouput channels. Default: 12
        R (int, optional): Multiplier of `growth_rate` to determine the number of ouput channels in the first convolutional layer. Default: 2
        kernel_size (int or list, optional): Size of the convolving kernels. Default: [1, 3]
        act_fn (torch.nn.Module, optional): Activation function used after each convolutional layer. Default: torch.nn.ReLU
        dropout_rate (float, optional): Probability for the dropout layer. Default: 0.2
        batchnorm (bool, optional): Include a BatchNorm1d layer. Default: False
    """
    def __init__(
        self,
        in_channels: int,
        layers: int = 6,
        growth_rate: int = 12,
        R: int = 2,
        kernel_size: Union[int, List[int]] = [1, 3],
        act_fn: Type[nn.Module] = nn.ReLU,
        dropout_rate: float = 0.2,
        batchnorm: bool = False
    ) -> None:
        super().__init__()
        self.in_channels = in_channels
        self.layers = layers
        self.growth_rate = growth_rate
        self.R = R
        self.kernel_size = kernel_size
        self.act_fn = act_fn
        self.dropout_rate = dropout_rate
        self.batchnorm = batchnorm
        self.out_channels = self.in_channels + self.layers * self.growth_rate
        layer_list = []
        for layer_idx in range(self.layers):
            layer_list.append(
                DenseLayer(
                    in_channels = self.in_channels + layer_idx * self.growth_rate,
                    growth_rate = self.growth_rate,
                    R = self.R,
                    kernel_size = self.kernel_size,
                    act_fn = self.act_fn,
                    dropout_rate = self.dropout_rate,
                    batchnorm = self.batchnorm
                )
            )
        self.net = nn.Sequential(*layer_list)

    def forward(self, x: Tensor) -> Tensor:
        out = self.net(x)
        return out

# Define a transition layer (compresses the output from the previous DenseBlock):
class TransitionLayer(nn.Module):
    """
    Compresses the output from the previous layer by passing it through a convolutional and an average pooling layer

    Args:
        in_channels (int): Number of input channels
        compression_factor (int, optional): Factor by which the output channels from the previous layer are compressed. Default: 2
        act_fn (torch.nn.Module, optional): Activation function used for the convolutional layer. Default: torch.nn.ReLU
        pool_size (int, optional): Kernel size for the average pooling layer. Default: 2
        batchnorm (bool, optional): Include a BatchNorm1d layer. Default: False
    """
    def __init__(
        self,
        in_channels: int,
        compression_factor: int = 2,
        act_fn: Type[nn.Module] = nn.ReLU,
        pool_size: int = 2,
        batchnorm: bool = False
    ) -> None:
        super().__init__()
        assert compression_factor <= in_channels, "`compression_factor` cannot be larger than `in_channels`"
        self.in_channels = in_channels
        self.compression_factor = compression_factor
        self.act_fn = act_fn
        self.pool_size = pool_size
        self.batchnorm = batchnorm
        self.out_channels = self.in_channels // self.compression_factor
        net = [
            nn.Conv1d(
                in_channels = self.in_channels,
                out_channels = self.out_channels,
                kernel_size = 1,
                stride = 1,
                padding = 'same'
            ),
            act_fn(),
            nn.AvgPool1d(kernel_size = self.pool_size)
        ]
        if self.batchnorm:
            net = [nn.BatchNorm1d(self.in_channels)] + net
        self.net = nn.Sequential(*net)

    def forward(self, x: Tensor) -> Tensor:
        out = self.net(x)
        return out

# Define a DenseNet module:
class DenseNet(nn.Module):
    """
    Generates a PyTorch module with a basic DenseNet architecture
    (i.e. groups of dense layers separated by transition layers)
    
    Args:
        in_channels (int): Number of input channels
        blocks (int, optional): Number of DenseBlocks. Default = 4
        layers (int or list, optional): Number of dense layers in each DenseBlock. Default = 6
        growth_rate (int or list, optional): Growth rate for each DenseBlock. Default = 12
        R (int or list, optional): Multiplier of `growth_rate` for each DenseBlock. Default = 2
        kernel_size (int or list, optional): Size of the convolving kernels for each DenseBlock. Default: [1, 3]
        act_fn (torch.nn.Module, optional): Activation function used after each convolutional layer. Default: torch.nn.ReLU
        dropout_rate (float or list, optional): Probability for each dropout layer. Default: 0.2
        compression_factor (int or list, optional): Compression factor for each transition layer. Default: 2
        pool_size (int or list, optional): Kernel size for each average pooling layer. Default: 2
        omit_final_transition (bool, optional): Omit the transition layer after the last DenseBlock. Default: True
        batchnorm (str, optional): Include a BatchNorm1d layer in the dense layers or transition layers. Default: None
    """
    def __init__(
        self,
        in_channels: int,
        blocks: int = 4,
        layers: Union[int, List[int]] = 6,
        growth_rate: Union[int, List[int]] = 12,
        R: Union[int, List[int]] = 2,
        kernel_size: Union[int, List[int]] = [[1, 3]],
        act_fn: Union[int, List[int]] = torch.nn.ReLU,
        dropout_rate: Union[float, List[float]] = 0.2,
        compression_factor: Union[int, List[int]] = 2,
        pool_size: Union[int, List[int]] = 2,
        omit_final_transition: bool = True,
        batchnorm: Union[str, None, List[Union[str, None]]] = None
    ) -> None:
        super().__init__()
        self.in_channels = in_channels
        self.blocks = blocks
        self.layers = _to_list(layers, self.blocks)
        self.growth_rate = _to_list(growth_rate, self.blocks)
        self.R = _to_list(R, self.blocks)
        self.kernel_size = _to_list(kernel_size, self.blocks)
        self.act_fn = _to_list(act_fn, self.blocks)
        self.dropout_rate = _to_list(dropout_rate, self.blocks)
        self.compression_factor = _to_list(compression_factor, self.blocks)
        self.pool_size = _to_list(pool_size, self.blocks)
        self.omit_final_transition = omit_final_transition
        self.batchnorm = _to_list(batchnorm, self.blocks)
        bn_opts = [None, 'none', 'dense', 'transition', 'both']
        for bn in self.batchnorm:
            assert bn in bn_opts, f"Invalid value for `batchnorm`: {bn!r}; should be one of {bn_opts}"
        net = []
        c_in = in_channels
        for block_idx in range(blocks):
            net.append(
                DenseBlock(
                    in_channels = c_in,
                    layers = self.layers[block_idx],
                    growth_rate = self.growth_rate[block_idx],
                    R = self.R[block_idx],
                    kernel_size = self.kernel_size[block_idx],
                    act_fn = self.act_fn[block_idx],
                    dropout_rate = self.dropout_rate[block_idx],
                    batchnorm = self.batchnorm[block_idx] in ['dense', 'both']
                )
            )
            c_in = c_in + self.layers[block_idx] * self.growth_rate[block_idx]
            if block_idx < (self.blocks - 1) or not omit_final_transition:
                net.append(
                    TransitionLayer(
                        in_channels = c_in,
                        compression_factor = self.compression_factor[block_idx],
                        pool_size = self.pool_size[block_idx],
                        act_fn = self.act_fn[block_idx],
                        batchnorm = self.batchnorm[block_idx] in ['transition', 'both']
                    )
                )
                c_in = c_in // self.compression_factor[block_idx]
        
        self.net = nn.Sequential(*net)
        self.out_channels = c_in

    def forward(self, x: Tensor) -> Tensor:
        out = self.net(x)
        return out
    
    def get_out_length(self, length_in: int) -> int:
        trans_layers = self.blocks - 1 if self.omit_final_transition else self.blocks
        length_out = length_in // 2**trans_layers
        return length_out

# Define a function to build the BiDen model:
class BiDen(nn.Module):
    """
    A DenseNet model with a bidirectional convolution as the first layer

    Args:
        
    """

    def __init__(
        self,
        input_len: int,
        output_dim: int,
        conv_kwargs: dict = {},
        dn_kwargs: dict = {},
        fc_kwargs: dict = {}
    ) -> None:
        super().__init__()
        self.input_len = input_len
        self.output_dim = output_dim
        self.conv_kwargs = conv_kwargs
        self.dn_kwargs = dn_kwargs
        self.fc_kwargs = fc_kwargs
        self.convnet  = semiBiConv(**self.conv_kwargs)
        output_len = self.convnet.get_out_length(input_len)
        self.densenet = DenseNet(in_channels = self.convnet.out_channels, **self.dn_kwargs)
        output_len = self.densenet.get_out_length(output_len)
        self.fcnet = nn.Linear(in_features = self.densenet.out_channels * output_len, out_features = self.output_dim, **self.fc_kwargs)

    def forward(self, x: Tensor) -> Tensor:
        out = self.convnet(x)
        out = self.densenet(out)
        out = torch.flatten(out, start_dim = 1)
        out = self.fcnet(out)
        return out