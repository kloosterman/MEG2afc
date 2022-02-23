#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 21:44:05 2016
conda activate hddm_env
ipcluster start
spyder

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

path = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/behav/HDDM'
os.chdir(path)
# data = hddm.load_csv('./ddmdat_histbias.csv')

model_name = 'chi_accuracy_basic_runs' # ddm_histbias_perrun ddm_acc_perrun
print(model_name)

subj = range(19)
for isub in subj:
    data = hddm.load_csv('./subj%d_ddmdat.csv' % isub)
    data = data.dropna()

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
    plt.savefig('RTs_%s_stimcoded.pdf' % isub)

    fig = plt.figure()
    ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions')
    for i, subj_data in data_acc_coded.groupby('subj_idx'):
        #print(subj_data.subj_idx)
        subj_data.rt.hist(bins=100, histtype='step', ax=ax)
    #plt.savefig('RTdistribution_{}.pdf'.format(phase))
    plt.savefig('RTs_%s_acc_coded.pdf' % isub)


    if model_name == 'chi_accuracy_basic_runs':
        m = hddm.HDDM(data_acc_coded,
        depends_on={'a': ['drug', 'simon', 'difficulty'],
        't': ['drug', 'simon', 'difficulty'],
        'v': ['drug', 'simon', 'difficulty']})
        
        m.find_starting_values()
        m.sample(1000,500)
        
        test = m.gen_stats()
        test.to_csv('subj%d_ddmparams.csv' % isub )  # './subj%d_ddmdat.csv' % isub)




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
    

