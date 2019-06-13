import math
import random

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as spstat
import tensorflow as tf
from abcpy.continuousmodels import Uniform
from abcpy.probabilisticmodels import (Continuous, InputConnector,
                                       ProbabilisticModel)
from numpy.random import RandomState

random_state = 25

class PollakDemand(ProbabilisticModel, Continuous): 
    """ Economic model that describes the relationship between demand and price
        :param list parameters: Contains the probabilistic models and hyperparameters from which the model derives.
                                Final value in list explains what the variable of interest is: return or markup
        :param string markup_or_return: Indicates whether we want to perform inference on Markup or Returns
    """
    eval_param = ""

    def __init__ (self, parameters, markup_or_return ,name="Pollak"):
        self.eval_param = markup_or_return
        if not isinstance(parameters, list):
            raise TypeError('Input of Pollak model is of type list')
        # We expect input of type parameters = [gamma, delta, sigma, ln_loc, ln_scale, ln_scale_underlying]
        if len(parameters) != 6:
            raise RuntimeError('Input list must be length 5, containing [gamma, delta, sigma, ln_loc, ln_scale, ln_scale_underlying]')
        
        input_connector = InputConnector.from_list(parameters)
        super().__init__(input_connector,name)

    def _check_input(self, input_values):
        sigma = input_values[2]
        ln_scale= input_values[4]
        ln_scale_underlying = input_values[5]
        eval_param = self.eval_param
        
        # Check whether input has correct type or format
        if len(input_values) != 6:
            raise ValueError('Number of parameters of Pollak Generative Model mut be 6')
        
        # Check whether input is from correct domain
        if  sigma<=0 or ln_scale<=0 or ln_scale_underlying<=0:
            return False
        
        if eval_param not in ['return', 'markup']:
            return False
        
        return True
    
    def _check_output(self, values):
        if not isinstance(values[0], np.ndarray):
            raise ValueError('Output of the normal distribution is always a numnber.')
        return True
    
    def get_output_dimension(self):
        return 1

    def forward_simulate(self, input_values, rng= np.random.RandomState()  ):
        
        #Extract input parameters
        gamma = input_values[0]
        delta = input_values[1]
        sigma = input_values[2]
        ln_loc = input_values[3]
        ln_scale = input_values[4]
        ln_scale_underlying = input_values[5]
        eval_param = self.eval_param
        
        #Do the actual forward simulation
        arr_estimated_output = self.generative_step(gamma, delta, sigma, ln_loc, ln_scale, ln_scale_underlying)
        
        if eval_param == 'return':
            #calculate a sample of revenue
            arr_eval_param = ( ( (arr_estimated_output - gamma)/delta )**sigma )* arr_estimated_output
        
        elif eval_param == 'markup':
            #calculate a sample of markups - use formula for marginal revenue
            arr_eval_param = ( ((arr_estimated_output-gamma)/delta)**(1/sigma +1) ) * (delta*sigma) * (arr_estimated_output*(sigma-1)-gamma*sigma)**-1
        
        # Format the ouput to obey API
        result = [np.array([x]) for x in arr_eval_param]
        return result

    def generative_step(self, gamma, delta, sigma, ln_loc, ln_scale, ln_scale_underlying):
        
        #Drawing 10000 samples form a log-normal distirbution for \phi with parameter \ln_loc, \ln_scale, ln_scale_underlying
        #This is the sample size for the generated output x during a forward simulation step of the ABC algorithm
        sample_size_phi= 10000
        arr_phi = spstat.lognorm.rvs(s=ln_scale_underlying, loc=ln_loc , scale=ln_scale, size=sample_size_phi, random_state=random_state)

        #Creating array to hold phi and associated params e.g. gamma, delta, sigma
        arr_params = np.array([gamma, delta, sigma]).repeat(sample_size_phi).reshape(3,-1).T 
        arr_data  = np.concatenate( [arr_phi.reshape(-1,1), arr_params] ,axis=1 ) 

        #Learns the inv function, and use it estimate values of x corresponding to our array of \phi, \gamma, \delta, \sigma
        #This is the sample size used by the neural network to estimate the inv function x(\phi)
        sample_size_output=2500
        arr_estimated_output = self.approximate_inv_fun(gamma, delta, sigma, sample_size_output, arr_data.astype(np.float32))
        
        plt.plot(arr_estimated_output.tolist()[:200], arr_phi.tolist()[:200])
        plt.savefig('estimated_function.pdf')

        return arr_estimated_output

    def approximate_inv_fun(self, gamma, delta, sigma, sample_size, data_to_predict):

        #draw samples from X as Uniform
        #TODO:The upper and lower bound need to be changed to reflect real data
        print('Generating a sample of x from uniform and \phi from x. \nSample to be used to calculate local inverse function x from \phi given values for gamma, delta and sigma')
        arr_x = np.random.uniform(low=6, high=8, size=sample_size)

        #calc a sample of phi from x
        arr_phi = ( ((arr_x-gamma)/delta)**((1+sigma)/sigma) ) * (delta*sigma)/(arr_x*(sigma-1)-(gamma*sigma))

        # plt.plot(arr_x, arr_phi)
        # plt.show()

        #creating NN model to estimate inv_func x(\phi)
        arr_vars = np.array([gamma, delta, sigma]).repeat(sample_size).reshape(3,-1).T 
        arr_data  = np.concatenate( [arr_phi.reshape(-1,1), arr_vars, arr_x.reshape(-1,1)], axis=1 ) #passing in x, gamma, delta, sigma as vars
        arr_data = np.split(arr_data, indices_or_sections= [int(.6*sample_size), int(.8*sample_size) ]  )
        
        train = arr_data[0]
        val = arr_data[1]
        test = arr_data[2]

        x_train = train[:, :-1]
        y_train = train[:, -1]
        x_val = val[:, :-1]
        y_val = val[:, -1] 
        x_test = test[:, :-1]
        y_test = test[:, -1]

        tf.set_random_seed(25)
        arr_estimated_output = self.train_neural_network(x_train, y_train, x_val, y_val, x_test, y_test, data_to_predict)

        return arr_estimated_output

    def train_neural_network(self, x_train, y_train, x_val, y_val, x_test, y_test, x_to_predict):
        """ Routine to train neural network"""

        # Define the learn rate and batch size and epochs
        learn_rate = 0.0005
        batch_size = 64
        hm_epochs = 30

        # Define input / output placeholders
        x = tf.placeholder('float', [None, 4],name='input')
        y = tf.placeholder('float')

        # Define other Tensors
        prediction = self.neural_network_model(x)
        cost = tf.reduce_mean(tf.square(prediction - y))
        optimizer = tf.train.AdamOptimizer(learn_rate).minimize(cost)

        #cycles feed forward + backprop
        with tf.Session() as sess:
            sess.run(tf.global_variables_initializer())
            print('Initialized, Proceeding to train Network and Estimate Inverse Function')

            # Train in each epoch with whole data
            for epoch in range(hm_epochs+1):
                for step in range(len(y_train)//batch_size): 
                    for inputX, inputY in self.get_batch(x_train, y_train, batch_size):
                        _, l = sess.run([optimizer, cost], feed_dict={x:inputX, y:inputY } )                        
                if epoch %5 ==0:
                    epoch_loss_train = sess.run(cost, feed_dict={x:x_train, y:y_train} )
                    epoch_loss_val  = sess.run(cost, feed_dict={ x:x_val, y:y_val } )
                    print("Epoch {0}\t Train loss: {1}\t Val loss: {2}".format(epoch, epoch_loss_train , epoch_loss_val ) )
            
            _ , epoch_loss_val  = sess.run([prediction, cost], feed_dict={ x:x_test, y:y_test } )
            print('Loss on Test Set:', epoch_loss_val )

            print('Using Learnt Inverse Function to estimate output x from phi')
            arr_estimated_output = sess.run([prediction], {x:x_to_predict} )
            arr_estimated_output = np.concatenate(arr_estimated_output).ravel()
            return arr_estimated_output
        
    def neural_network_model(self, data):
        """ Routine to compute the Neural Network """
        
        # Define the number of nodes in each layer
        n_nodes_hl1 = 16
        n_nodes_hl2 = 8
        n_nodes_hl3 = 8
        n_nodes_hl4 = 8

        n_input = 4
        n_classes = 1

        # Define the layers
        hidden_1_layer = {'weights': tf.Variable(name='w_h1', initial_value=tf.random_normal([n_input, n_nodes_hl1], stddev=0.05)),
                          'biases': tf.Variable(name='b_h1', initial_value=tf.random_normal([n_nodes_hl1], stddev=0.05) )}

        hidden_2_layer = {'weights': tf.Variable(name='w_h2', initial_value=tf.random_normal([n_nodes_hl1, n_nodes_hl2], stddev=0.05)),
                          'biases': tf.Variable(name='b_h2', initial_value=tf.random_normal([n_nodes_hl2], stddev=0.05))}

        hidden_3_layer = {'weights': tf.Variable(name='w_h3', initial_value=tf.random_normal([n_nodes_hl2,n_nodes_hl3], stddev=0.05)),
                          'biases': tf.Variable(name='b_h3', initial_value=tf.random_normal([n_nodes_hl3], stddev=0.05))}

        hidden_4_layer = {'weights': tf.Variable(name='w_h4', initial_value=tf.random_normal([n_nodes_hl3,n_nodes_hl4], stddev=0.05)),
                          'biases': tf.Variable(name='b_h4', initial_value=tf.random_normal([n_nodes_hl4], stddev=0.05))}
                          
        output_layer = {'weights': tf.Variable(name='w_o', initial_value=tf.random_normal([n_nodes_hl4, n_classes], stddev=0.05)),
                        'biases': tf.Variable(name='b_o', initial_value=tf.random_normal([n_classes], stddev=0.05)) }

        # (input_data * weights) + biases
        l1 = tf.add(tf.matmul(data, hidden_1_layer['weights']), hidden_1_layer['biases'] )
        l1 = tf.nn.relu(l1)

        l2 = tf.add(tf.matmul(l1, hidden_2_layer['weights']), hidden_2_layer['biases'] )
        l2 = tf.nn.relu(l2)

        l3 = tf.add(tf.matmul(l2, hidden_3_layer['weights']), hidden_3_layer['biases'] )
        l3 = tf.nn.relu(l3)

        l4 = tf.add(tf.matmul(l3, hidden_4_layer['weights']), hidden_4_layer['biases'] )
        l4 = tf.nn.relu(l4)

        output = tf.add(tf.matmul(l4, output_layer['weights']), output_layer['biases'], name='output')

        return output

    def get_batch(self, inputX, inputY, batch_size):
        duration = len(inputX)
        for i in range(0, duration//batch_size):
            idx = i*batch_size + np.random.randint(0,10,(1))[0]

            yield inputX[idx:idx+batch_size], inputY[idx:idx+batch_size]
