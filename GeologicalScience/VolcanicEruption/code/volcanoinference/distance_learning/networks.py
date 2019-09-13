import torch.nn as nn
import torch.nn.functional as F


class EmbeddingNet(nn.Module):
    """"""

    def __init__(self):
        super(EmbeddingNet, self).__init__()
        # put some fully connected layers:
        self.fc1 = nn.Linear(72, 100)
        self.fc2 = nn.Linear(100, 80)
        self.fc3 = nn.Linear(80, 40)
        self.fc4 = nn.Linear(40, 15)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = self.fc4(x)
        return x

    def get_embedding(self, x):
        return self.forward(x)


class SiameseNet(nn.Module):
    """From https://github.com/adambielski/siamese-triplet"""

    def __init__(self, embedding_net):
        super(SiameseNet, self).__init__()
        self.embedding_net = embedding_net

    def forward(self, x1, x2):
        output1 = self.embedding_net(x1)
        output2 = self.embedding_net(x2)
        return output1, output2

    def get_embedding(self, x):
        return self.embedding_net(x)


class TripletNet(nn.Module):
    """From https://github.com/adambielski/siamese-triplet"""

    def __init__(self, embedding_net):
        super(TripletNet, self).__init__()
        self.embedding_net = embedding_net

    def forward(self, x1, x2, x3):
        output1 = self.embedding_net(x1)
        output2 = self.embedding_net(x2)
        output3 = self.embedding_net(x3)
        return output1, output2, output3

    def get_embedding(self, x):
        return self.embedding_net(x)


class EmbeddingNetSummaryStatistics(nn.Module):
    """ """

    def __init__(self):
        super(EmbeddingNetSummaryStatistics, self).__init__()
        # put some fully connected layers:
        self.fc1 = nn.Linear(72, 80)
        self.fc2 = nn.Linear(80, 40)
        self.fc3 = nn.Linear(40, 15)
        self.fc4 = nn.Linear(15, 2)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = self.fc4(x)
        return x

    def get_embedding(self, x):
        return self.forward(x)
