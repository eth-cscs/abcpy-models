import torch
# add to the pythonpath the current folder (otherwise it does not find the script files)
import os, sys

sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))

from algorithms import contrastive_training, triplet_training, FP_nn_training
from utilities import read_simulations, compute_similarity_matrix

cuda = torch.cuda.is_available()
device = "cuda" if cuda else "cpu"

print(device)


def save_net(path, net):
    torch.save(net.state_dict(), path)


param1, param2, output = read_simulations('training_data.txt')

# we use the whole dataset now, we do not keep any test example
epsilon = 0.6
similarity_set = compute_similarity_matrix(param1, param2, epsilon)

print("Contrastive")
embedding_net_contrastive = contrastive_training(output, similarity_set, cuda, positive_weight=0.4, n_epochs=400)
save_net("../saved-networks/contrastive.pth", embedding_net_contrastive.cpu())

print("Triplet")
embedding_net_triplet = triplet_training(output, similarity_set, cuda, n_epochs=800)
save_net("../saved-networks/triplet.pth", embedding_net_triplet.cpu())

cuda = False
print("FP-nn")
embedding_net_FP_nn = FP_nn_training(output, param1, param2, cuda, n_epochs=400)
save_net("../saved-networks/FP_nn.pth", embedding_net_FP_nn.cpu())
