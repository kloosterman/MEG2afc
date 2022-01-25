% function MEG2afc_setup_5D_freqstats_corr()
% Make jobs for each condition to run freqanalysis

% restoredefaultpath
if ismac
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman'; % on the cluster
    backend = 'none'; % local torque
else
    basepath = '/home/mpib/kloosterman'; % on the cluster
    addpath(fullfile(basepath, 'MATLAB',  'tools/qsub_tardis'))
    backend = 'torque'; % local torque
end
addpath(genpath(fullfile(basepath, 'MATLAB',   'MEG_HH_analysis')))
rmpath(genpath(fullfile(basepath, 'MATLAB', 'MEG_HH_analysis/.git/')))
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip-20161220')) %inc JJ edit ft_artifact_zvalue
ft_defaults
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip-20161220/external/spm8/')) %for spm_bwlabel


%% get blink stamps from EOG
preproc = load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/preproc/NK12/plac_ipsi/NK12_ses1_20141213_plac_ipsi_run1_stim_data.mat');
EOGart = ft_findcfg(preproc.data.cfg, 'artfctdef');

%%
datapath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/data/NK12/A/meg';

cfg=[];
cfg.channel = {'EEG057' 'EEG058'};
cfg.dataset = fullfile( datapath, 'NK12_NielsDet_20141213_01.ds');
cfg.checkmaxfilter=false;
cfg.continuous = 'yes';


event = ft_read_event(cfg.dataset);

trgval1 = [event(find(strcmp('UPPT001',{event.type}))).value];
trgsmp1 = [event(find(strcmp('UPPT001',{event.type}))).sample];

start = find(trgval1 == 127);
start_smp = trgsmp1(start);

data = ft_preprocessing(cfg);

[~,~,start_timestamp,~,end_timestamp] = import_blink_file('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/MEG_HH_analysis/eye_counts/nk12_1_blink.csv');

blinks_EL = [start_timestamp(2:end) end_timestamp(2:end)]   / 1000 * 1200;
% shift to match Eyelink to eye data
blinks_EL = round(blinks_EL + start_smp);


%%
close all
cfg2=[];
% cfg2.trl = blinks_EL;
cfg2.artfctdef.blinks_EyeLink.artifact = blinks_EL;
cfg2.artfctdef.EOG_ver.artifact = EOGart.eog_ver.artifact;
cfg2.artfctdef.EOG_hor.artifact = EOGart.eog_hor.artifact;
cfg2.demean = 'yes';
cfg2.blocksize = 100;
cfg3 = ft_databrowser(cfg2, data)
