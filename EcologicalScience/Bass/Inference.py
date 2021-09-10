import matplotlib.pyplot as plt
import numpy as np
from abcpy.statistics import Identity
import logging
logging.basicConfig(level=logging.INFO)
from model import Bass, LogNormal

problem = 'bass'
sample = False
plot = True
predict = False
algorithm = 'sabc' #'sabc'
dircertory = 'Figures'
# Define backend
from abcpy.backends import BackendDummy, BackendMPI
backend = BackendDummy()
#backend = BackendMPI()
# file=dircertory + "/"+'rejection_one_pcent.csv'
# rejection_abc = np.loadtxt(open(file, "rb"), delimiter=",", skiprows=1)
# print(np.mean(rejection_abc, axis=0))

# Define Graphical Model
from abcpy.continuousmodels import Uniform
H = Uniform([[0.25], [0.75]], name='H')
Am = Uniform([[0.000280], [0.00059]], name='Am')
AE = Uniform([[0], [0.003]], name='AE')
PM = Uniform([[0.045], [0.135]], name='PM')
I = Uniform([[2.5], [7.5]], name='I')
# H = LogNormal([[np.log(0.5)], [1]], name='H')
# Am = LogNormal([[np.log(0.00043)], [1]], name='Am')
# AE = LogNormal([[np.log(0.0015)], [1]], name='AE')
# PM = LogNormal([[np.log(0.03712)], [1]], name='PM')
# I = LogNormal([[np.log(5.00)], [1]], name='I')



Bass = Bass([H, Am, AE, PM, I], name='Bass')

# Example to Generate Data to check it's correct
#resultfakeobs1 = Bass.forward_simulate([.5,0.00043,0.0015,0.06,5.00], 1)
obs_data = [np.genfromtxt('data/valid_data_Warwick.csv', delimiter=',', skip_header=1)[:,1:].transpose().flatten().reshape(-1,)]

# Define Statistics
statistics_calculator = Identity(degree=1, cross=False)
#print('# Check whether the statistis works')
#print(statistics_calculator.statistics(resultfakeobs1))
#print(statistics_calculator.statistics(obs_data))
# Define distance
from Distance import WeightedEuclidean
wt = list(1000*np.ones(11)) + list(0*np.ones(11))
for ind in range(60):
    wt = wt + list((1/30)*np.ones(11))
distance_calculator = WeightedEuclidean(statistics_calculator, weight=np.array(wt))
#print('# Check whether the distance works')
#print(distance_calculator.distance(obs_data, resultfakeobs1))

if sample:
   # Define kernel
   from abcpy.perturbationkernel import DefaultKernel
   kernel = DefaultKernel([H, Am, AE, PM, I])

   if algorithm == 'sabc':
        ## SABC ##
        from abcpy.inferences import SABC
        sampler = SABC([Bass], [distance_calculator], backend, kernel, seed=1)
        steps, epsilon, n_samples, n_samples_per_param, ar_cutoff, full_output, journal_file = 10, 10e20, 111, 1, 0.001, 1, None
        print('SABC Inferring')
        # We use resultfakeobs1 as our observed dataset
        journal_sabc = sampler.sample(observations=[obs_data], steps=steps, epsilon=epsilon, n_samples = n_samples,
                                        n_samples_per_param = n_samples_per_param, beta = 2, delta = 0.2, v = 0.3,
                                        ar_cutoff = 0.001, resample = None, n_update = None, full_output = 1,
                                      journal_file = journal_file)
        print(journal_sabc.posterior_mean())
        journal_sabc.save(dircertory + "/" + algorithm +'_'+problem+'_obs.jrnl')
   elif algorithm == 'apmcabc':
        #APMCABC
        from abcpy.inferences import APMCABC
        sampler = APMCABC([Bass], [distance_calculator], backend, kernel, seed = 1)
        steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file =10, 111, 1, 0.6, 0.001, 2, 1, None
        print('APMCABC Inferring')
        journal_apmcabc = sampler.sample([obs_data], steps, n_samples, n_samples_per_param, alpha, acceptance_cutoff, covFactor, full_output, journal_file)
        print(journal_apmcabc.posterior_mean())
        journal_apmcabc.save('apmcabc_'+problem+'_obs.jrnl')

if predict:
    # from model import DrawFromPosteriorFile
    # dp = DrawFromPosteriorFile([Bass], backend)
    # print('Prediction')
    # parameters, simulations = dp.sample(file=dircertory + "/"+'rejection_one_pcent.csv')
    # np.savez(dircertory + "/"+'rejection_predict.npz', parameters=parameters, simulations=simulations)

    from model import DrawFromPosterior
    dp = DrawFromPosterior([Bass], backend)
    print('Prediction')
    parameters, simulations = dp.sample(journal_file=dircertory + "/" + algorithm +'_'+problem+'_obs.jrnl')
    np.savez(dircertory + "/" + algorithm +'_'+'predict.npz', parameters=parameters, simulations=simulations)

if plot:
    from abcpy.output import Journal
    journal_sabc = Journal.fromFile(dircertory + "/" +algorithm+"_" + problem + '_obs.jrnl')

    ### Converting posterior samples to a matrix
    accepted_parameters = journal_sabc.get_accepted_parameters(-1)
    parameters = []
    for ind in range(111):
        parameters.append([x[0] for x in accepted_parameters[ind]])
    parameters = np.array(parameters)
    parameters[:,-1] = parameters[:,-1] * 1e+13
    np.savetxt("Figures/sabc_parameters.csv", np.array(parameters), delimiter=",")

    print(journal_sabc.configuration)
    print(journal_sabc.posterior_mean())
    print(journal_sabc.posterior_cov())
    distance = []
    for ind in range(10):
        distance.append(journal_sabc.distances[ind])
    plt.figure()
    plt.boxplot(np.array(distance).transpose())
    plt.savefig(dircertory + '/'+algorithm+'_'+'distance.pdf')
    journal_sabc.plot_posterior_distr(path_to_save=dircertory + "/"+algorithm+"_" + problem + '_param.pdf',
                                       show_samples=False, parameters_to_show=['H','Am', 'AE', 'PM', 'I'])
    fig, ax = journal_sabc.plot_ESS()
    fig.savefig(dircertory + '/'+algorithm+'_'+'ess.pdf')

    predobs_sabc = np.load(dircertory + "/" + 'sabc' +'_'+'predict.npz')['simulations']
    #predobs_rej = np.load(dircertory + "/" + 'rejection' + '_' + 'predict.npz')['simulations']
    #predobs_try = np.loadtxt(open('Figures/checkrito.csv', "rb"), delimiter=",", skiprows=1)[:,:11]

    #print(predobs_rej.shape)
    print(predobs_sabc.shape)
    #print(predobs_try.shape)

    yearly_mean_var = np.zeros(shape=(4,682))
    yearly_mean_var[0,:] = np.mean(predobs_sabc, axis=0)
    yearly_mean_var[1, :] = np.var(predobs_sabc, axis=0)
    yearly_mean_var[2, :] = np.median(predobs_sabc, axis=0)
    yearly_mean_var[3, :] = np.subtract(*np.percentile(predobs_sabc, [75, 25], axis=0))

    np.savetxt("Figures/yearly_mean_var.csv", yearly_mean_var, delimiter=",")


    x_val = 2004 + np.linspace(0,10,11)
    ### SSB
    plt.figure()
    plt.plot(x_val, obs_data[0][:11], 'k.-', label='Observed')
    plt.plot(x_val, np.median(predobs_sabc[:, :11], axis=0), 'r', label='(SABC) Median predicted')
    plt.fill_between(x_val,  np.min(predobs_sabc[:, :11], axis=0), np.max(predobs_sabc[:, :11], axis=0), color = 'r', alpha=.05)
    plt.fill_between(x_val,  np.quantile(predobs_sabc[:, :11], .025, axis=0), np.quantile(predobs_sabc[:, :11], .975, axis=0), color = 'r', alpha=.1)
    # plt.plot(x_val, np.median(predobs_rej[:, :11], axis=0), 'y', label='(Rejection-ABC) Median predicted')
    # plt.fill_between(x_val,  np.min(predobs_rej[:, :11], axis=0), np.max(predobs_rej[:, :11], axis=0), color = 'y', alpha=.05)
    # plt.fill_between(x_val,  np.quantile(predobs_rej[:, :11], .025, axis=0), np.quantile(predobs_rej[:, :11], .975, axis=0), color = 'y', alpha=.1)
    #plt.ylabel('SSB', fontsize=30)
    plt.ylim(ymin=0)
    # plt.xticks(fontsize=30)
    # plt.yticks(fontsize=30)
    plt.title('SSB', fontsize=30, loc='right')
    plt.grid()
    #plt.legend()
    plt.tight_layout()
    plt.savefig(dircertory + '/Pred_SSB.pdf')
    plt.close()

    ### Rec
    plt.figure()
    plt.plot(x_val,obs_data[0][11:22], 'k.-', label='Observed')
    plt.plot(x_val,np.median(predobs_sabc[:, 11:22], axis=0), 'r', label='(SABC) Median predicted')
    plt.fill_between(x_val, np.min(predobs_sabc[:, 11:22], axis=0), np.max(predobs_sabc[:, 11:22], axis=0), color = 'r', alpha=.05)
    plt.fill_between(x_val, np.quantile(predobs_sabc[:, 11:22], .025, axis=0), np.quantile(predobs_sabc[:, 11:22], .975, axis=0), color = 'r', alpha=.1)

    # plt.plot(x_val,np.median(predobs_rej[:, 11:22], axis=0), 'g', label='(Rejection) Median predicted')
    # plt.fill_between(x_val, np.min(predobs_rej[:, 11:22], axis=0), np.max(predobs_rej[:, 11:22], axis=0), color = 'g', alpha=.05)
    # plt.fill_between(x_val, np.quantile(predobs_rej[:, 11:22], .025, axis=0), np.quantile(predobs_rej[:, 11:22], .975, axis=0), color = 'g', alpha=.1)
    #plt.ylabel('Rec', fontsize=30)
    plt.ylim(ymin=0)
    plt.title('Rec', fontsize=30, loc='right')
    #plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.savefig(dircertory + '/Pred_Rec.pdf')
    plt.close()
    for ind_M in range(30):
        plt.figure()
        plt.plot(x_val,obs_data[0][2*11+(ind_M*11):2*11+((ind_M+1)*11)], 'k.-', label='Observed')
        plt.plot(x_val,np.median(predobs_sabc[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], axis=0), 'r', label='(SABC) Median predicted')
        plt.fill_between(x_val, np.min(predobs_sabc[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], axis=0), np.max(predobs_sabc[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], axis=0), color = 'r', alpha=.05)
        plt.fill_between(x_val, np.quantile(predobs_sabc[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], 0.025, axis=0), np.quantile(predobs_sabc[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], 0.975, axis=0), color = 'r', alpha=.1)

        # plt.plot(x_val,np.median(predobs_rej[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], axis=0), 'g', label='(Rejection-ABC) Median predicted')
        # plt.fill_between(x_val, np.min(predobs_rej[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], axis=0), np.max(predobs_rej[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], axis=0), color = 'g', alpha=.05)
        # plt.fill_between(x_val, np.quantile(predobs_rej[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], 0.025, axis=0), np.quantile(predobs_rej[:, 2*11+(ind_M*11):2*11+((ind_M+1)*11)], 0.975, axis=0), color = 'g', alpha=.1)
        plt.title('M' + str(ind_M), fontsize=30, loc='right')
        plt.ylim(ymin=0)
        #plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig(dircertory + '/pred_M_' + str(ind_M) + '.pdf')
        plt.close()
    for ind_N in range(30):
        plt.figure()
        plt.plot(x_val,obs_data[0][32*11+(ind_N*11):32*11+((ind_N+1)*11)], 'k.-', label='Observed')
        plt.plot(x_val,np.median(predobs_sabc[:, 32*11+(ind_N*11):32*11+((ind_N+1)*11)], axis=0), 'r', label='(SABC) Median predicted')
        plt.fill_between(x_val, np.min(predobs_sabc[:, 32*11+(ind_N*11):32*11+((ind_N+1)*11)], axis=0),
                         np.max(predobs_sabc[:, 32*11+(ind_N*11):32*11+((ind_N+1)*11)], axis=0), color = 'r', alpha=.05)
        plt.fill_between(x_val, np.quantile(predobs_sabc[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], .025, axis=0),
                         np.quantile(predobs_sabc[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], 0.975, axis=0), color='r',
                         alpha=.1)

        # plt.plot(x_val, np.median(predobs_rej[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], axis=0), 'g',
        #          label='(Rejection-ABC) Median predicted')
        # plt.fill_between(x_val, np.min(predobs_rej[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], axis=0),
        #                  np.max(predobs_rej[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], axis=0),
        #                  color='g', alpha=.05)
        # plt.fill_between(x_val, np.quantile(predobs_rej[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], .025,
        #                                     axis=0),
        #                  np.quantile(predobs_rej[:, 32 * 11 + (ind_N * 11):32 * 11 + ((ind_N + 1) * 11)], 0.975,
        #                              axis=0), color='g',
        #                  alpha=.1)
        plt.title('N'+str(ind_N), fontsize=30, loc='right')
        plt.ylim(ymin=0)
        #plt.legend()
        plt.grid()
        plt.tight_layout()
        plt.savefig(dircertory + '/pred_N_'+str(ind_N)+'.pdf')
        plt.close()



