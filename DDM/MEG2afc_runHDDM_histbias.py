#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 21:44:05 2016

@author: kloosterman
"""
#reset

import pandas as pd
import matplotlib.pyplot as plt
import hddm
import os
import numpy

# # to get plots on screen
# %matplotlib

path = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM'
os.chdir(path)
data = hddm.load_csv('./ddmdat_histbias.csv')

# rename column values
data.loc[data.drug==0, 'drug'] = 'plac'
data.loc[data.drug==1, 'drug'] = 'atx'
data.loc[data.simon==0, 'simon'] = 'ipsi'
data.loc[data.simon==1, 'simon'] = 'contra'
data.loc[data.difficulty==0, 'difficulty'] = 'easy'
data.loc[data.difficulty==1, 'difficulty'] = 'hard'
data.loc[data.prevresp==0, 'prevresp'] = 'Lprev'
data.loc[data.prevresp==1, 'prevresp'] = 'Rprev'

data = data[data.rt > 0.2] # drop too fast RT's
# data = data[data.subj_idx != 1] # drop NK2

# # for historybias analysis
# data_stimcoded = hddm.utils.flip_errors(data)

# move correct col to response col, keep actual response in button, flip data to make errors neg
data_acc_coded = hddm.utils.flip_errors( data.rename(columns={'response': 'button', 'correct': 'response'}) )

# for historybias analysis
# data_stimcoded = hddm.utils.flip_errors( data.rename(columns={'response': 'button', 'correct': 'response'}) )
# no error RT flipping, needed??
data_stimcoded =  data.rename(columns={'response': 'button', 'correct': 'response'})
# # % ompolen response contra session: choice
# data_stimcoded.response = data_stimcoded.button
data_stimcoded = data_stimcoded.rename(columns={'response':'correct', 'button': 'response' })
data_stimcoded.loc[data_stimcoded.simon=='contra', 'response'] = data_stimcoded[data_stimcoded.simon=='contra'].response.replace({0:1, 1:0})



fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data_stimcoded.groupby('subj_idx'):
    #print(subj_data.subj_idx)
    subj_data.rt.hist(bins=100, histtype='step', ax=ax)
#plt.savefig('RTdistribution_{}.pdf'.format(phase))
plt.savefig('RTdistribution_stimcoded.pdf')

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
for i, subj_data in data_acc_coded.groupby('subj_idx'):
    #print(subj_data.subj_idx)
    subj_data.rt.hist(bins=100, histtype='step', ax=ax)
#plt.savefig('RTdistribution_{}.pdf'.format(phase))
plt.savefig('RTdistribution_acc_coded.pdf')


model_name = 'ddm_histbias_perrun_ol' # ddm_histbias_perrun ddm_acc_perrun
print(model_name)


# Regular Chi square DDMs
if model_name == 'chi_prevresp_z_dc':
    # TODO plot fits etc
    subj_params = []
    for subj_idx, subj_data in data_stimcoded.groupby('subj_idx'):
        m_subj = hddm.HDDMStimCoding(subj_data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
        depends_on={'a': ['drug','simon'], 
        't': ['drug','simon'], 
        'v': ['drug','simon', 'difficulty'], 
        'dc':['drug','simon', 'prevresp'], 
        'z': ['drug','simon', 'prevresp']})       
        subj_params.append(m_subj.optimize(method='chisquare'))    
    params = pd.DataFrame(subj_params)
    params.to_csv('{}.csv'.format(model_name))        

if model_name == 'chi_prevresp_z_dc_nomotor':
    # TODO plot fits etc
    subj_params = []
    for subj_idx, subj_data in data_stimcoded.groupby('subj_idx'):
        m_subj = hddm.HDDMStimCoding(subj_data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
        depends_on={'a': ['drug'], 
        't': ['drug'], 
        'v': ['drug', 'difficulty'], 
        'dc':['drug', 'prevresp'], 
        'z': ['drug', 'prevresp']})       
        subj_params.append(m_subj.optimize(method='chisquare'))    
    params = pd.DataFrame(subj_params)
    params.to_csv('{}.csv'.format(model_name))        

elif model_name == 'chi_prevresp_z_dc_runs': 
    for subj_idx, subj_data in data_stimcoded.groupby('subj_idx'):
        subj_params = []
        m_subj = hddm.HDDMStimCoding(subj_data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True,
        depends_on={'a': ['drug','simon', 'run_nr'],
        't': ['drug','simon', 'run_nr'],
        'v': ['drug','simon', 'difficulty', 'run_nr'],
        'dc':['drug','simon', 'prevresp', 'run_nr'],
        'z': ['drug','simon', 'prevresp', 'run_nr']})       
        subj_params.append(m_subj.optimize(method='chisquare'))    
        params = pd.DataFrame(subj_params)
        params.to_csv('{}_subj{}.csv'.format(model_name, subj_idx))        

elif model_name == 'chi_accuracy_basic':
    subj_params = []
    for subj_idx, subj_data in data_acc_coded.groupby('subj_idx'):
        m_subj = hddm.HDDM(subj_data, 
        depends_on={'a': ['drug','simon',], 
        't': ['drug','simon',], 
        'v': ['drug','simon', 'difficulty']})     
        subj_params.append(m_subj.optimize(method='chisquare'))
    params = pd.DataFrame(subj_params)
    params.to_csv('{}.csv'.format(model_name))      
      

elif model_name == 'chi_accuracy_basic_nomotor':
    subj_params = []
    for subj_idx, subj_data in data_acc_coded.groupby('subj_idx'):
        m_subj = hddm.HDDM(subj_data, 
        depends_on={'a': ['drug',], 
        't': ['drug',], 
        'v': ['drug', 'difficulty']})       
        subj_params.append(m_subj.optimize(method='chisquare'))    
    params = pd.DataFrame(subj_params)
    params.to_csv('{}.csv'.format(model_name))        

# ori, let HDDM do the looping
elif model_name == 'chi_accuracy_basic_runs':
    for subj_idx, subj_data in data_acc_coded.groupby('subj_idx'):
        subj_params = []
        m_subj = hddm.HDDM(subj_data,
        depends_on={'a': ['drug','simon', 'run_nr'],
        't': ['drug','simon', 'run_nr'],
        'v': ['drug','simon', 'difficulty', 'run_nr']})
        subj_params.append(m_subj.optimize(method='chisquare'))
        params = pd.DataFrame(subj_params)
        params.to_csv('{}_subj{}.csv'.format(model_name, subj_idx))

elif model_name == 'ddm_acc_perrun' or model_name == 'ddm_acc_perrun_ol':
    fig = plt.figure(figsize=(8.27,11.69)); 
    plt.rcParams.update({'font.size': 7})   
    iplot=0    
    for subj_idx, subj_data in data_acc_coded.groupby('subj_idx'):
    # for subj_idx, subj_data in data_acc_coded.groupby('subj_idx').filter(lambda x: (x.subj_idx == 1).any()):
        # todo start run_idx from 0?
        for drug_idx, drug_data in subj_data.groupby('drug'):
            for ses_idx, ses_data in drug_data.groupby('simon'):
                subj_params = []
                sim_data = []                                     
                for run_idx, run_data in ses_data.groupby('run_nr'):
                    if model_name == 'ddm_histbias_perrun_ol':                        
                         rt = run_data['rt'].abs()
                         rtz = (rt - rt.mean()) / rt.std()
                         run_data = run_data[rtz < 3]                    
                    m_subj = hddm.HDDM(run_data, depends_on={'v': 'difficulty'})     
                    run = m_subj.optimize(method='chisquare')
                    subj_params.append(run)   
         
                    sim_params = {'easy': {'v': run['v(easy)'], 'a': run['a'], 't': run['t']},
                                  'hard': {'v': run['v(hard)'], 'a': run['a'], 't': run['t']}} 
                    sim_data_run, sim_params = hddm.generate.gen_rand_data(sim_params, size=800)
                    sim_data_run['drug'] = drug_idx
                    sim_data_run['simon'] = ses_idx
                    sim_data_run['run_nr'] = run_idx
                    sim_data.append(sim_data_run)

                sim_data = pd.concat(sim_data)
                
                sim_data = hddm.utils.flip_errors(sim_data)                
                iplot+=1
                # curplot = numpy.ravel_multi_index((subj_idx,iplot), (19,4))
                ax = fig.add_subplot(19,4,iplot, xlabel='RT', ylabel='count', title= ' Subj{} {} {}'.format(subj_idx, drug_idx, ses_idx) )
                # ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title= ' Subj{} {} {}'.format(subj_idx, drug_idx, ses_idx) )
                ses_data.rt.hist(bins=numpy.linspace(-3,3,50), histtype='step', normed=True, ax=ax)
                sim_data.rt.hist(bins=numpy.linspace(-3,3,50), histtype='step', normed=True, ax=ax)
                ax.grid(False)
                
                # #save params and sim data to csv:
                subj_params = pd.DataFrame(subj_params)
                subj_params.to_csv('params_{}_subj{}_{}_{}.csv'.format(model_name, subj_idx, drug_idx, ses_idx))
                # sim_data.to_csv('sim_data_{}_subj{}.csv'.format(model_name, subj_idx))
        plt.savefig('RThist_{}.pdf'.format(model_name))
        if subj_idx == 1:
            iplot+=1
            
                                
elif model_name == 'ddm_histbias_perrun' or model_name == 'ddm_histbias_perrun_ol':
    fig = plt.figure(figsize=(8.27,11.69)); 
    plt.rcParams.update({'font.size': 7})   
    iplot=0    
    for subj_idx, subj_data in data_stimcoded.groupby('subj_idx'):
        for drug_idx, drug_data in subj_data.groupby('drug'):
            for ses_idx, ses_data in drug_data.groupby('simon'):
                subj_params = []
                sim_data = []                                     
                for run_idx, run_data in ses_data.groupby('run_nr'):
                    if model_name == 'ddm_histbias_perrun_ol':                        
                         rt = run_data['rt'].abs()
                         rtz = (rt - rt.mean()) / rt.std()
                         run_data = run_data[rtz < 3]
                    # m_subj = hddm.HDDMStimCoding(run_data,  stim_col='stimulus', split_param='v', depends_on={'v': 'prevresp', 'z': 'prevresp'}, bias=True)
                    m_subj = hddm.HDDMStimCoding(run_data,  stim_col='stimulus', split_param='v', depends_on={'v': 'difficulty', 'dc': 'prevresp', 'z': 'prevresp'}, drift_criterion=True, bias=True)
                                            
                    run = m_subj.optimize(method='chisquare')
                    subj_params.append(run)   
         
                #     # TODO FIX simulations with v: difficulty in there
                #     sim_params = {'Lprev': {'v': run['v'], 'z': run['z(Lprev)'], 'dc': run['dc(Lprev)'], 'a': run['a'], 't': run['t']},
                #                   'Rprev': {'v': run['v'], 'z': run['z(Rprev)'], 'dc': run['dc(Rprev)'], 'a': run['a'], 't': run['t']}
                #               }
                #     # sim_params = {'Lprev': {'v': run['dc(Lprev)'],'z': run['z(Lprev)'], 'a': run['a'], 't': run['t']},
                #     #               'Rprev': {'v': run['dc(Rprev)'],'z': run['z(Rprev)'], 'a': run['a'], 't': run['t']}
                #     #           }
                #     sim_data_run, sim_params = hddm.generate.gen_rand_data(sim_params, size=800)
                #     sim_data_run['drug'] = drug_idx
                #     sim_data_run['simon'] = ses_idx
                #     sim_data_run['run_nr'] = run_idx
                #     sim_data.append(sim_data_run)

                # sim_data = pd.concat(sim_data)
                
                # sim_data = hddm.utils.flip_errors(sim_data)                
                # # ses_data = hddm.utils.flip_errors(ses_data)
                # ses_data =  data_acc_coded[data_acc_coded.subj_idx == subj_idx]
                # iplot+=1
                # # curplot = numpy.ravel_multi_index((subj_idx,iplot), (19,4))
                # ax = fig.add_subplot(19,4,iplot, xlabel='RT', ylabel='count', title= ' Subj{} {} {}'.format(subj_idx, drug_idx, ses_idx) )
                # # ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title= ' Subj{} {} {}'.format(subj_idx, drug_idx, ses_idx) )
                # ses_data.rt.hist(bins=numpy.linspace(-3,3,50), histtype='step', normed=True, ax=ax)
                # sim_data.rt.hist(bins=numpy.linspace(-3,3,50), histtype='step', normed=True, ax=ax)
                # ax.grid(False)
                
                # #save params and sim data to csv:
                subj_params = pd.DataFrame(subj_params)
                subj_params.to_csv('params_{}_subj{}_{}_{}.csv'.format(model_name, subj_idx, drug_idx, ses_idx))
                # sim_data.to_csv('sim_data_{}_subj{}.csv'.format(model_name, subj_idx))
        plt.savefig('RThist_{}.pdf'.format(model_name))
        if subj_idx == 1:
            iplot+=1


elif model_name == 'ddm_histbias_perses':
    fig = plt.figure(figsize=(8.27,11.69)); 
    plt.rcParams.update({'font.size': 7})   
    iplot=0    
    for subj_idx, subj_data in data_stimcoded.groupby('subj_idx'):
        for drug_idx, drug_data in subj_data.groupby('drug'):
            for ses_idx, ses_data in drug_data.groupby('simon'):
                subj_params = []
                sim_data = []                                     
                # for run_idx, run_data in ses_data.groupby('run_nr'):
                # m_subj = hddm.HDDMStimCoding(run_data,  stim_col='stimulus', split_param='v', depends_on={'v': 'prevresp', 'z': 'prevresp'}, bias=True)
                m_subj = hddm.HDDMStimCoding(ses_data,  stim_col='stimulus', split_param='v', depends_on={'dc': 'prevresp', 'z': 'prevresp'}, drift_criterion=True, bias=True)     
                
                run = m_subj.optimize(method='chisquare')
                subj_params.append(run)   
     
                # sim_params = {'Lprev': {'v': run['v(Lprev)'],'z': run['z(Lprev)'], 'a': run['a'], 't': run['t']},
                #               'Rprev': {'v': run['v(Rprev)'],'z': run['z(Rprev)'], 'a': run['a'], 't': run['t']}
                #           }
                # sim_params = {'Lprev': {'v': run['dc(Lprev)'],'z': run['z(Lprev)'], 'a': run['a'], 't': run['t']},
                #               'Rprev': {'v': run['dc(Rprev)'],'z': run['z(Rprev)'], 'a': run['a'], 't': run['t']}
                #           }
                sim_params = {'Lprev': {'v': run['v'], 'z': run['z(Lprev)'], 'dc': run['dc(Lprev)'], 'a': run['a'], 't': run['t']},
                              'Rprev': {'v': run['v'], 'z': run['z(Rprev)'], 'dc': run['dc(Rprev)'], 'a': run['a'], 't': run['t']}
                        }
                sim_data_run, sim_params = hddm.generate.gen_rand_data(sim_params, size=800)
                sim_data_run['drug'] = drug_idx
                sim_data_run['simon'] = ses_idx
                # sim_data_run['run_nr'] = run_idx
                sim_data.append(sim_data_run)

                sim_data = pd.concat(sim_data)
                
                sim_data = hddm.utils.flip_errors(sim_data)
                ses_data =  data_acc_coded[data_acc_coded.subj_idx == subj_idx]
                                
                iplot+=1
                # curplot = numpy.ravel_multi_index((subj_idx,iplot), (19,4))
                ax = fig.add_subplot(19,4,iplot, xlabel='RT', ylabel='count', title= ' Subj{} {} {}'.format(subj_idx, drug_idx, ses_idx) )
                # ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title= ' Subj{} {} {}'.format(subj_idx, drug_idx, ses_idx) )
                ses_data.rt.hist(bins=numpy.linspace(-3,3,50), histtype='step', normed=True, ax=ax)
                sim_data.rt.hist(bins=numpy.linspace(-3,3,50), histtype='step', normed=True, ax=ax)
                ax.grid(False)
                
                # #save params and sim data to csv:
                subj_params = pd.DataFrame(subj_params)
                subj_params.to_csv('params_{}_subj{}_{}_{}.csv'.format(model_name, subj_idx, drug_idx, ses_idx))
                # sim_data.to_csv('sim_data_{}_subj{}.csv'.format(model_name, subj_idx))
            
        plt.savefig('RThist_{}.pdf'.format(model_name))










# TODO hddm versions
elif model_name == 'hddm_prevresp_z_dc':
    #%%
    # Instantiate model object passing it our data (no need to call flip_errors() before passing it).
    # This will tailor an individual hierarchical DDM around your dataset.
    #m = hddm.HDDM(data)
    # find a good starting point which helps with the convergence.
    #m.find_starting_values()
    ## start drawing 7000 samples and discarding 5000 as burn-in
    #m.sample(2000, burn=20)

    # m = hddm.HDDM(data, depends_on={'v': 'age', 'a': 'age', 't': 'age'}, bias=False) # , include='all'
    # let drift depend on both drug and pupil bin


    # Regular Chi square DDM
    subj_params = []
    for subj_idx, subj_data in data.groupby('subj_idx'):
        # m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'occ_alphabin'], 'a':['drug', 'occ_alphabin'], 't':['drug', 'occ_alphabin'] }, p_outlier=0.05,)
        # m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'occ_alphabin', 'diff'], 'a':['drug', 'occ_alphabin'], 't':['drug', 'occ_alphabin'] }, p_outlier=0.05,)
        # m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'diff'], 'a':['drug'], 't':['drug'] }, p_outlier=0.05,)
        # m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'diff'], 'a':['drug', 'diff'], 't':['drug', 'diff'] }, p_outlier=0.05,)
        # m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'pupilbin0_4'], 'a':['drug', 'pupilbin0_4'], 't':['drug', 'pupilbin0_4'] }, p_outlier=0.05,)
        # subj_params.append(m_subj.optimize(method='chisquare', quantiles=(.1, .2, .3, .4, .5, .6, .7, .8, .9 )))

        # own code
        m_subj = hddm.HDDMStimCoding(subj_data, stim_col='stimulus', split_param='v',
        drift_criterion=True, bias=True,
        depends_on={'a': ['drug'], 't': ['drug'], 'v': ['drug', 'difficulty'], 'dc':['drug', 'prevresp'], 'z':['drug', 'prevresp']})
        # Urai code
        # m_subj = hddm.HDDMStimCoding(subj_data, stim_col='stimulus', split_param='v',
        #    drift_criterion=True, bias=True, p_outlier=0.05,
        #   include=('sv', 'sz'), group_only_nodes=['sv', 'sz'],
        #  depends_on={'v': ['difficulty'], 'dc':['prevresp'], 'z':['prevresp']})      
    
        subj_params.append(m_subj.optimize(method='chisquare'))

        params = pd.DataFrame(subj_params)

        params.to_csv('params_histbias.csv' )


## m = hddm.HDDM(data, depends_on={'v':['drug', 'occ_alphabin'], 'a':['drug', 'occ_alphabin'], 't':['drug', 'occ_alphabin'] }, bias=False) # , include='all'
#m = hddm.HDDM(data, depends_on={'v':['drug', 'diff'], 'a':['drug'], 't':['drug'] } , bias=False) # , include='all'
##m = hddm.HDDM(data, depends_on={'v':['drug', 'diff'], 'a':['drug'], 't':['drug'], 'sv':['drug'] } , bias=False) # , include='all'
#m.find_starting_values()
##m.sample(10000, burn=1000)
#m.sample(10000, burn=10000/10)
##m.sample(100, burn=100/10)
#m.plot_posterior_predictive(figsize=(20, 10))
## plt.show()
##
#test = m.gen_stats()
#test.to_csv('HDDM_drift_drugdiff.csv' )

# TODO sort out regression model, st like this:
# m = hddm.models.HDDMRegressor(data, 'v ~ alpha:C(drug)')

#v_easyatx, v_easyplac= m.nodes_db.node[['v(easy.atx)', 'v(easy.plac)']]
#print ("P_v(v_easyatx > v_easyplac) = ", (v_easyatx.trace() > v_easyplac.trace()).mean())
#hddm.analyze.plot_posterior_nodes([v_easyatx, v_easyplac])
#
#v_hardatx, v_hardplac= m.nodes_db.node[['v(hard.atx)', 'v(hard.plac)']]
#print ("P_v(v_hardatx > v_hardplac) = ", (v_hardatx.trace() > v_hardplac.trace()).mean())
#hddm.analyze.plot_posterior_nodes([v_hardatx, v_hardplac])


##TODO make loop across ddm para's
##paralist=["v" "a" "t"]
#paralist=[ 'v', 'a', 't' ]
#paranamelist = [ 'drift-rate', 'boundary separation', 'non-decision time' ]
## TODO f, ax = plt.subplots(2,2)
#for idx, ipara, in enumerate(paralist):    
#    v_0, v_1 = m.nodes_db.node[['{}(0)'.format(ipara), '{}(1)'.format(ipara)]]    
#    hddm.analyze.plot_posterior_nodes([v_0, v_1])
#    plt.xlabel(ipara)
#    plt.ylabel('Posterior probability')
#    plt.title('Posterior of {} group means, p = {}'.format(paranamelist[idx], (v_0.trace() < v_1.trace()).mean()))
#    plt.savefig('hddm_{}_{}.pdf'.format(phase, paranamelist[idx]))
#   
                                   
                                   
#%% OLD

# # TODO find out how to plot young and old in different figures
# m.plot_posterior_predictive(figsize=(20, 20))
# plt.savefig('hddm_subject_fits {}.pdf'.format(phase))
#
#
# #%%
#
# #TODO make loop across ddm para's
# #paralist=["v" "a" "t"]
# paralist=[ 'v', 'a', 't' ]
# paranamelist = [ 'drift-rate', 'boundary separation', 'non-decision time' ]
# # TODO f, ax = plt.subplots(2,2)
# for idx, ipara, in enumerate(paralist):
#     v_0, v_1 = m.nodes_db.node[['{}(0)'.format(ipara), '{}(1)'.format(ipara)]]
#     hddm.analyze.plot_posterior_nodes([v_0, v_1])
#     plt.xlabel(ipara)
#     plt.ylabel('Posterior probability')
#     plt.title('Posterior of {} group means, p = {}'.format(paranamelist[idx], (v_0.trace() < v_1.trace()).mean()))
#     plt.savefig('hddm_{}_{}.pdf'.format(phase, paranamelist[idx]))
    

