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

path = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM'
os.chdir(path)

#phase = 'test'
phase = 'study'

#data = hddm.load_csv('./EyeMem_hddm_{}.csv'.format(phase))
#data = hddm.load_csv('./EyeMem_hddm_study.csv')
data = hddm.load_csv('./data_alltrials_3bins_outliersout2_char.csv')
#data.head(10)

data = data[data.rt > 0.2] # drop too fast RT's

#data = data[data.subj_idx != 5] # drop subj_idx 5, only one going in opposite dir

data = hddm.utils.flip_errors(data)

#fig = plt.figure()
#ax = fig.add_subplot(111, xlabel='RT', ylabel='count', title='RT distributions {}'.format(phase))
#for i, subj_data in data.groupby('subj_idx'):
##    print(subj_data.subj_idx)
#    subj_data.rt.hist(bins=20, histtype='step', ax=ax)


#plt.savefig('RTdistribution_{}.pdf'.format(phase))
# plt.savefig('RTdistribution.pdf')

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
    m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'diff'], 'a':['drug', 'diff'], 't':['drug', 'diff'] }, p_outlier=0.05,)
    # m_subj = hddm.HDDM(subj_data, depends_on={'v':['drug', 'pupilbin0_4'], 'a':['drug', 'pupilbin0_4'], 't':['drug', 'pupilbin0_4'] }, p_outlier=0.05,)
    # subj_params.append(m_subj.optimize(method='chisquare', quantiles=(.1, .2, .3, .4, .5, .6, .7, .8, .9 )))
    subj_params.append(m_subj.optimize(method='chisquare'))
params = pd.DataFrame(subj_params)

params.to_csv('params_drug_diff2.csv' )



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
    

