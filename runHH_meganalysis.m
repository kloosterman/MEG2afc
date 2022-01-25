% function runHH_meganalysis

% restoredefaultpath
if ismac
    basepath = '/Users/kloosterman/Documents/GitHub/'; % local
    backend = 'none'; % local torque2    addpath(fullfile(basepath, 'MEG2afc'))
    addpath(genpath(fullfile(basepath, 'plotting_tools/')))
    addpath(genpath(fullfile(basepath, 'stats_tools/')))
else
    basepath = '/mnt/beegfs/home/kloosterman/GitHub'; % on the cluster
%     addpath(fullfile(basepath, 'tools'))
    backend = 'slurm'; % local torque slurm
end
% addpath(genpath(fullfile(basepath, 'MEG_HH_analysis')))
% rmpath(genpath(fullfile(basepath, 'MEG_HH_analysis/.git/')))
addpath(genpath_exclude(fullfile(basepath, 'MEG2afc'), {'\.git'} ))

% addpath(fullfile(basepath, 'tools', 'fieldtrip-20170611')) %inc JJ edit ft_artifact_zvalue
addpath(fullfile(basepath, 'fieldtrip')) % cloned on 13 09 19
addpath(fullfile(basepath, 'zapline-plus')) 
addpath(fullfile(basepath, 'qsub-tardis')) %inc JJ edit ft_artifact_zvalue
% addpath(fullfile(basepath, 'tools/qsub_tardis_slurmpreview'))

ft_defaults
% addpath(fullfile(basepath, 'tools', 'mmse')) %
% addpath(fullfile(basepath, 'tools', 'dva')) 
% addpath(fullfile(basepath, 'tools', 'statfun')) 
% addpath(fullfile(basepath, 'tools/NoiseTools')) % robust detrend (MSE)
% addpath(genpath(fullfile(basepath, 'tools', 'pls')))
% addpath(fullfile(basepath,  'critEEG_analysis/preprocessing/')) % for remove_ERP_fromdata function

% addpath(fullfile(basepath, 'tools/LCMV_pipelines')) % var based artf rejection oostenveld paper

%% FOR spectral latr history bias paper
%%
%% new smooth short functions 
%%
%%
%% preprocessing
MEG2afc_preproc_setup()

%% behavioral analysis
% behavior = MEG2afc_readbehavior_setup()

%% load behavior
% load '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/behav/behavstruct.mat'; 

%% plot behavior TODO this
% MEG2afc_behaviorplot(behavior)

%% freqanalysis
% MEG2afc_freq_setup()

%% load and combine freqs of sessions
% [megdat] = MEG2afc_mergefreq( 'BLC' ); % 'raw' or BLC, TODO get both, will fit

%% corr pow with behavior
% megdat = MEG2afc_mergefreq_corrstats(megdat);

%% plot corr stats, script for now
% MEG2afc_mergefreq_corrstats_plot %(megdat)

%% run drug vs plac stats
% [megdat] = MEG2afc_mergefreq_stats(megdat);

%% plot 3D stats
% MEG2afc_mergefreq_plotstats(megdat)

%%
%%
%%

