from abcpy.statistics import Statistics
import numpy as np
import torch


def save_net(path, net):
    torch.save(net.state_dict(), path)
    
def load_net(path, network_class):
    net = network_class()
    net.load_state_dict(torch.load(path))
    return net.eval()  # call the network to eval model. Needed with batch normalization and dropout layers. 

class NeuralEmbeddingStatistics(Statistics):
    """
    """

    def __init__(self, net, degree=1, cross=False):  # are these default values OK?
        self.net = net
        self.degree = degree
        self.cross = cross
        
    def statistics(self, data):      
        if isinstance(data, list):
            if np.array(data).shape == (len(data),):
                if len(data) == 1:
                    data = np.array(data).reshape(1, 1)
                data = np.array(data).reshape(len(data), 1)
            else:
                data = np.concatenate(data).reshape(len(data), -1)
        else:
            raise TypeError('Input data should be of type list, but found type {}'.format(type(data)))
        
        # simply apply the network transformation. 
        result = self.net(torch.from_numpy(data.astype("float32"))).cpu().detach().numpy()
        
        # Expand the data with polynomial expansion            
        result = self._polynomial_expansion(result)
        return np.array(result)

class NeuralEmbeddingStatisticsFromFile(NeuralEmbeddingStatistics):
    """
    """

    def __init__(self, path_to_net_state_dict, network_class, degree=1, cross=False):
        self.network_class = network_class
        self.net = load_net(path_to_net_state_dict, self.network_class)
        self.degree = degree
        self.cross = cross