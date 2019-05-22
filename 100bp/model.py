"""
Adapted from DeepSEA (Zhou & Troyanskaya, 2015).
"""
import numpy as np
import torch
import torch.nn as nn


class DeeperDeepSEA(nn.Module):
    """
    A deeper DeepSEA model architecture.

    Parameters
    ----------
    sequence_length : int
        The length of the sequences on which the model trains and and makes
        predictions.
    n_targets : int
        The number of targets (classes) to predict.

    Attributes
    ----------
    conv_net : torch.nn.Sequential
        The convolutional neural network component of the model.
    classifier : torch.nn.Sequential
        The linear classifier and sigmoid transformation components of the
        model.

    """

    def __init__(self, sequence_length, n_targets):
        super(DeeperDeepSEA, self).__init__()
        conv_kernel_size = 8
        pool_kernel_size = 4

        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(32, 32, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.BatchNorm1d(32),

            nn.Conv1d(32, 48, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(48, 48, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.BatchNorm1d(48),
            nn.Dropout(p=0.2),

            nn.Conv1d(48, 96, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Conv1d(96, 96, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(96),
            nn.Dropout(p=0.2))

        reduce_by = 2 * (conv_kernel_size - 1)
        pool_kernel_size = float(pool_kernel_size)
        self._n_channels = int(
            np.floor(
                (np.floor(
                    (sequence_length - reduce_by) / pool_kernel_size)
                 - reduce_by) / pool_kernel_size)
            - reduce_by)
        self.classifier = nn.Sequential(
            nn.Linear(96 * self._n_channels, n_targets),
            nn.ReLU(inplace=True),
            nn.BatchNorm1d(n_targets),
            nn.Linear(n_targets, n_targets),
            nn.Sigmoid())

    def forward(self, x):
        """
        Forward propagation of a batch.
        """
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 96 * self._n_channels)
        predict = self.classifier(reshape_out)
        return predict

def criterion():
    """
    Specify the appropriate loss function (criterion) for this
    model.

    Returns
    -------
    torch.nn._Loss
    """
    return nn.BCELoss()

def get_optimizer(lr):
    """
    Specify an optimizer and its parameters.

    Returns
    -------
    tuple(torch.optim.Optimizer, dict)
        The optimizer class and the dictionary of kwargs that should
        be passed in to the optimizer constructor.

    """
    return (torch.optim.SGD,
            {"lr": lr, "weight_decay": 1e-6, "momentum": 0.9})
