from time import time

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from metric_learn import SDML
from torch.optim import lr_scheduler
from torch.utils.data import Dataset

from distance_learning.datasets import Similarities, SiameseSimilarities, TripletSimilarities, ParameterSimulationPairs
from distance_learning.losses import ContrastiveLoss, TripletLoss
from distance_learning.networks import EmbeddingNet, SiameseNet, TripletNet, EmbeddingNetSummaryStatistics
from distance_learning.trainer import fit


def sdml_fit(samples, similarity_set, prior='covariance', balance_param=0.15):
    """Prior can be 'covariance', 'identity' or 'random'. 
    balance_param was used 0.5 in the first version of the paper, but it does not work here with such a large value. """

    n_samples = len(similarity_set)

    sdml = SDML(prior=prior, preprocessor=samples, verbose=True, balance_param=balance_param)

    pairs, Y = [], []
    for ind1 in range(n_samples):
        for ind2 in range(n_samples):
            pairs.append([ind1, ind2])
            if similarity_set[ind1, ind2]:
                Y.append(1)
            else:
                Y.append(-1)

    start = time()
    sdml.fit(pairs, Y)
    print("Fitting took {:.2f} seconds.".format(time() - start))

    return sdml


def contrastive_training(samples, similarity_set, cuda, batch_size=16, n_epochs=200, embedding_net=None,
                         positive_weight=None, load_all_data_GPU=True, margin=None, lr=None, optimizer=None,
                         scheduler=None, optimizer_kwargs={}, scheduler_kwargs={}, loader_kwargs={}):
    """If the dataset is small enough, we can speed up training by loading all on the GPU at beginning, by usinng
    load_all_data_GPU=True. It may crash if the dataset is too large. Note that in some cases using only CPU may still
    be quicker."""

    # Do all the setups

    # need to use the Similarities and SiameseSimilarities datasets

    similarities_dataset = Similarities(samples, similarity_set, "cuda" if cuda and load_all_data_GPU else "cpu")
    pairs_dataset = SiameseSimilarities(similarities_dataset, positive_weight=positive_weight)

    if cuda:
        if load_all_data_GPU:
            loader_kwargs_2 = {'num_workers': 0, 'pin_memory': False}
        else:
            loader_kwargs_2 = {'num_workers': 1, 'pin_memory': True}
    else:
        loader_kwargs_2 = {}

    pairs_train_loader = torch.utils.data.DataLoader(pairs_dataset, batch_size=batch_size, shuffle=True,
                                                     **loader_kwargs, **loader_kwargs_2)

    if margin is None:
        margin = 1.

    if embedding_net is None:
        embedding_net = EmbeddingNet()
    model_contrastive = SiameseNet(embedding_net)

    if cuda:
        model_contrastive.cuda()
    loss_fn = ContrastiveLoss(margin)

    if lr is None:
        lr = 1e-3

    if optimizer is None:  # default value
        optimizer = optim.Adam(embedding_net.parameters(), lr=lr)
    else:
        optimizer = optimizer(embedding_net.parameters(), lr=lr, **optimizer_kwargs)

    if scheduler is None:  # default value, i.e. a dummy scheduler
        scheduler = lr_scheduler.StepLR(optimizer, 8, gamma=1, last_epoch=-1)
    else:
        scheduler = scheduler(optimizer, **scheduler_kwargs)

    # now train:
    fit(pairs_train_loader, model_contrastive, loss_fn, optimizer, scheduler, n_epochs, cuda)

    return embedding_net


def triplet_training(samples, similarity_set, cuda, batch_size=16, n_epochs=400, embedding_net=None,
                     load_all_data_GPU=True, margin=None, lr=None, optimizer=None, scheduler=None, optimizer_kwargs={},
                     scheduler_kwargs={}, loader_kwargs={}):
    """If the dataset is small enough, we can speed up training by loading all on the GPU at beginning, by usinng
    load_all_data_GPU=True. It may crash if the dataset is too large. Note that in some cases using only CPU may still
    be quicker."""

    # Do all the setups

    # need to use the Similarities and TripletSimilarities datasets

    similarities_dataset = Similarities(samples, similarity_set, "cuda" if cuda and load_all_data_GPU else "cpu")
    triplets_dataset = TripletSimilarities(similarities_dataset)

    if cuda:
        if load_all_data_GPU:
            loader_kwargs_2 = {'num_workers': 0, 'pin_memory': False}
        else:
            loader_kwargs_2 = {'num_workers': 1, 'pin_memory': True}
    else:
        loader_kwargs_2 = {}

    triplets_train_loader = torch.utils.data.DataLoader(triplets_dataset, batch_size=batch_size, shuffle=True,
                                                        **loader_kwargs, **loader_kwargs_2)

    if margin is None:
        margin = 1.

    if embedding_net is None:
        embedding_net = EmbeddingNet()
    model_triplet = TripletNet(embedding_net)

    if cuda:
        model_triplet.cuda()
    loss_fn = TripletLoss(margin)

    if lr is None:
        lr = 1e-3

    if optimizer is None:  # default value
        optimizer = optim.Adam(embedding_net.parameters(), lr=lr)
    else:
        optimizer = optimizer(embedding_net.parameters(), lr=lr, **optimizer_kwargs)

    if scheduler is None:  # default value, i.e. a dummy scheduler
        scheduler = lr_scheduler.StepLR(optimizer, 8, gamma=1, last_epoch=-1)
    else:
        scheduler = scheduler(optimizer, **scheduler_kwargs)

    # now train:
    fit(triplets_train_loader, model_triplet, loss_fn, optimizer, scheduler, n_epochs, cuda)

    return embedding_net


def FP_nn_training(samples, param1, param2, cuda, batch_size=1, n_epochs=50, embedding_net=None, load_all_data_GPU=True,
                   lr=None, optimizer=None, scheduler=None, optimizer_kwargs={}, scheduler_kwargs={}, loader_kwargs={}):
    """If the dataset is small enough, we can speed up training by loading all on the GPU at beginning, by usinng
    load_all_data_GPU=True. It may crash if the dataset is too large. Note that in some cases using only CPU may still
    be quicker."""

    # Do all the setups

    target = np.concatenate((param1.reshape(-1, 1), param2.reshape(-1, 1)), axis=1)
    dataset_FP_nn = ParameterSimulationPairs(samples, target, "cuda" if cuda and load_all_data_GPU else "cpu")

    if cuda:
        if load_all_data_GPU:
            loader_kwargs_2 = {'num_workers': 0, 'pin_memory': False}
        else:
            loader_kwargs_2 = {'num_workers': 1, 'pin_memory': True}
    else:
        loader_kwargs_2 = {}

    data_loader_FP_nn = torch.utils.data.DataLoader(dataset_FP_nn, batch_size=batch_size, shuffle=True, **loader_kwargs,
                                                    **loader_kwargs_2)

    if embedding_net is None:
        embedding_net = EmbeddingNetSummaryStatistics()

    if cuda:
        embedding_net.cuda()
    loss_fn = nn.MSELoss(reduction="mean")

    if lr is None:
        lr = 1e-3

    if optimizer is None:  # default value
        optimizer = optim.Adam(embedding_net.parameters(), lr=lr)
    else:
        optimizer = optimizer(embedding_net.parameters(), lr=lr, **optimizer_kwargs)

    if scheduler is None:  # default value, i.e. a dummy scheduler
        scheduler = lr_scheduler.StepLR(optimizer, 8, gamma=1, last_epoch=-1)
    else:
        scheduler = scheduler(optimizer, **scheduler_kwargs)

    # now train:
    fit(data_loader_FP_nn, embedding_net, loss_fn, optimizer, scheduler, n_epochs, cuda)

    return embedding_net
