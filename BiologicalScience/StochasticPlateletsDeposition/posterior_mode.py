import numpy as np
from scipy import stats
from scipy.optimize import minimize
from abcpy.output import Journal
import os

def compute_posterior_mode(whichobs, dim=None, type='multiply', data_type='_obs_' ):
    print(data_type)
    journal_sabc = Journal.fromFile('Journals/sabc'+ data_type + str(whichobs) + '_'+type+'.jrnl')
    # filenamere = '../APMCABC_Results/apmcabc_obs_' + str(whichobs) + '_re.jrnl'
    # if os.path.isfile(filenamere):
    #     filename = '../APMCABC_Results/apmcabc_obs_'+str(whichobs)+'_re.jrnl'
    # else:
    #     filename = '../APMCABC_Results/apmcabc_obs_'+str(whichobs)+'_reweighted.jrnl'
    # journal_sabc = Journal.fromFile(filename)
    ### Converting posterior samples to a matrix
    weights = np.concatenate(journal_sabc.get_weights(-1))
    accepted_parameters = journal_sabc.get_accepted_parameters(-1)
    parameters = []
    for ind in range(511):
        parameters.append(np.array([x[0] for x in accepted_parameters[ind]]).reshape(-1, 1))
    parameters = np.array(parameters).squeeze()
    if dim is not None:
        parameters = parameters[:, dim]
    #print(weights.shape, parameters.shape)
    kernel = stats.gaussian_kde(parameters.transpose(), weights=weights, bw_method=0.45)
    def rosen(x):
        return -np.log(kernel(x))
    posterior_mean = np.mean(parameters, axis=0)
    res = minimize(rosen, posterior_mean,  method='nelder-mead', options={'xatol': 1e-8, 'disp': False})
    return(np.array(res.x).reshape(1,-1))

compute = False
plot = True
test = False
pvalue_corrected = False
clustering = False
renewed = list(range(6)) + list(range(23,65))

if compute:
    ### Multiply ###
    print(renewed)
    joint_posterior_mode = np.zeros(shape=(65, 7))
    marginal_posterior_mode = np.zeros(shape=(65, 7))
    for whichobs in renewed:
        print(whichobs)
        joint_posterior_mode[whichobs, :] = compute_posterior_mode(whichobs, dim=None)
        for dim_ind in range(7):
            marginal_posterior_mode[whichobs,dim_ind] = compute_posterior_mode(whichobs, dim=dim_ind)
        print(joint_posterior_mode[whichobs, :], marginal_posterior_mode[whichobs, :])
        np.savez('posterior_mode_multiply', joint_posterior_mode=joint_posterior_mode, marginal_posterior_mode=marginal_posterior_mode)

if plot:
    #### Post-processing (Boxplot) ################
    import pylab as plt
    mode = np.load('posterior_mode_multiply.npz')['joint_posterior_mode']
    mode = mode[renewed, :]
    print(mode.shape)
    #namestitle = [r'$p_{Ad}$',r'$p_{Ag}$',r'$p_T$',r'$p_F$',r'$a_T$','v_z_AP','v_z_NAP']
    namestitle = [r'$p_{Ag}$',r'$a_T$',r'$v^{NAP}_{z}$']
    #limit = [[40, 120], ]
    fig, axs = plt.subplots(1, 3, figsize=(30,10))
    for ind, ax in enumerate(axs.ravel()):
        ax.set_title(namestitle[ind], fontsize=30)
        ax.boxplot([ mode[32:,ind],mode[16:32,ind], mode[:16,ind]], labels = ['Healthy', 'Dialysis', 'COPD'])
        ax.tick_params(axis='y', which='major', labelsize=30)
        ax.tick_params(axis='x', which='major', labelsize=30)
    plt.savefig('FinalFigures/Fig2.pdf')

if test:
    #### Joint: Post-processing (KW-test, Boxplot) ################
    print('Joint')
    mode = np.load('posterior_mode_multiply.npz')['joint_posterior_mode']
    mode = mode[renewed, :]
    print(mode.shape)
    ### KW-Test ##########
    from scipy import stats
    names = ['pAD','pAg','pT','pF','aT','v_z_AP','v_z_NAP']
    # # Analysis of MAP
    # print('Testing using MAP estimate')
    def formating(x):
        return "{:.2e}".format(x)
    namestitle = ['$\Pad$', '$\Pg$', '$\Pt$', '$\Pf$', '$\Ra$', r'$\vap$', r'$\vnap$']
    for ind in range(7):
        A = namestitle[ind]
        #print('Kruskal-Wallis test between healty, diabetic and COPD using parameter '+names[ind])
        copd, dialysis, healthy = mode[:16,ind], mode[16:32,ind], mode[32:,ind]
        if stats.kruskal(healthy, dialysis, copd).pvalue > 0.05:
            A = A + ' & ' +formating(stats.kruskal(healthy, dialysis, copd).statistic) + ' (' +  formating(stats.kruskal(healthy, dialysis, copd).pvalue) + ')'
        else:
            A = A + r' & \textbf{' + formating(stats.kruskal(healthy, dialysis, copd).statistic) + ' (' + formating(stats.kruskal(healthy, dialysis, copd).pvalue) + ')}'
        if stats.kruskal(healthy, dialysis).pvalue > 0.05:
            A = A + ' & ' + formating(stats.kruskal(healthy, dialysis).statistic) + ' (' + formating(stats.kruskal(healthy, dialysis).pvalue) + ')'
        else:
            A = A + r' & \textbf{' + formating(stats.kruskal(healthy, dialysis).statistic) + ' (' + formating(stats.kruskal(healthy, dialysis).pvalue) + ')}'
        if stats.kruskal(healthy, copd).pvalue > 0.05:
            A = A + ' & ' + formating(stats.kruskal(healthy, copd).statistic) + ' (' + formating(stats.kruskal(healthy, copd).pvalue) + ')'
        else:
            A = A + r' & \textbf{' + formating(stats.kruskal(healthy, copd).statistic) + ' (' + formating(stats.kruskal(healthy, copd).pvalue) + ')}'
        if stats.kruskal(dialysis, copd).pvalue > 0.05:
            A = A + ' & ' + formating(stats.kruskal(dialysis, copd).statistic) + ' (' + formating(stats.kruskal(dialysis, copd).pvalue) + r')\\'
        else:
            A = A + r' & \textbf{' + formating(stats.kruskal(dialysis, copd).statistic) + ' (' + formating(stats.kruskal(dialysis, copd).pvalue) + r')}\\'
        print(A)

if pvalue_corrected:
    correction_method = 'fdr_bh'
    alpha_value = 0.05
    #All three classes
    p = [2.56e-01, 9.51e-03, 9.88e-01, 4.9e-01, 8.26e-02, 1.56e-01, 1.42e-02]
    import statsmodels.stats.multitest as multitest
    p_corrected = multitest.multipletests(p, alpha=alpha_value, method=correction_method)
    print(p_corrected)
    #Healthy vs Dialysis
    p = [1.32e-01, 6.24e-01, 9.7e-01, 6.78e-01, 5.00e-02, 9.73e-02, 8.99e-02]
    import statsmodels.stats.multitest as multitest
    p_corrected = multitest.multipletests(p, alpha=alpha_value, method=correction_method)
    print(p_corrected)
    #Healthy vs COPD
    p = [8.51e-01, 6.66e-03, 8.8e-01, 5.22e-01, 9.1e-01, 8.51e-01, 1.13e-01]
    import statsmodels.stats.multitest as multitest
    p_corrected = multitest.multipletests(p, alpha=alpha_value, method=correction_method)
    print(p_corrected)
    #Dialysis vs COPD
    p = [1.87e-01, 1.16e-02, 9.1e-01, 2.14e-01, 5.95e-02, 9.73e-02, 6.66e-03]
    import statsmodels.stats.multitest as multitest
    p_corrected = multitest.multipletests(p, alpha=alpha_value, method=correction_method)
    print(p_corrected)


if clustering:
    # #### Joint: Post-processing (Clustering) ################
    mode = np.load('posterior_mode_multiply.npz')['joint_posterior_mode']
    mode = mode[renewed, :]
    print(mode.shape)
    from abcpy.statistics import Identity
    stat = Identity(degree=4, cross=True)
    mode = stat.statistics([[mode[i,:]] for i in range(mode.shape[0])])
    print(mode.shape)
    label = [1 for i in range(16)] + [2 for i in range(16)] + [3 for i in range(16)]
    print(label)
    from metric_learn import LMNN
    metric = LMNN(init='auto', k=6, min_iter=10000, max_iter=50000, convergence_tol=1e-6, learn_rate=1e-10,
                  regularization=.5, n_components=2)
    metric.fit(mode, label)
    L = metric.components_
    np.savez('L_all_3_cross_parameters.npz', L=L)

    L = np.load('L_all_3_cross_parameters.npz')['L']
    mode_lmnn = mode.dot(L.T)
    print(mode_lmnn.shape)
    import pylab as plt

    plt.figure()
    plt.plot(mode_lmnn[:16, 0], mode_lmnn[:16, 1], 'k*', label='Patients with COPD')
    plt.plot(mode_lmnn[16:32, 0], mode_lmnn[16:32, 1], 'r*', label='Patients having dialysis')
    plt.plot(mode_lmnn[32:, 0], mode_lmnn[32:, 1], 'b*', label='Healthy volunteers')
    plt.xlabel(r'$\theta_1$', fontsize=20)
    plt.ylabel(r'$\theta_2$', fontsize=20)
    # plt.title('2-dim projected space', fontsize=20)
    plt.legend(fontsize=10)# import pylab as plt
    plt.savefig('FinalFigures/Fig8.pdf')

    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics.cluster import adjusted_rand_score

    cluster_lmnn = AgglomerativeClustering(n_clusters=3, affinity='euclidean', linkage='ward')
    cluster_lmnn.fit_predict(mode_lmnn)
    print('AR score of Euclidean distance btwn LMNN: ' + str(adjusted_rand_score(label, cluster_lmnn.labels_)))

