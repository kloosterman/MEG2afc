function MEG2afc_mse_analysis(cfg)

% load('/Volumes/LNDG/user/Niels/MEG2afc/preproc/NK7/drug_ipsi/NK7_ses2_20141128_drug_ipsi_run3_stim_data.mat')
% load('/home/mpib/kloosterman/projectdata/MEG2afc/preproc/NK7/drug_ipsi/NK7_ses2_20141128_drug_ipsi_run3_stim_data.mat')
% data_ori = data;

% cfg = [];
% cfg.channel = 'MEG';
% ft_databrowser(cfg, data)

%% detrend
for itrial=1:length(data.trial)
    data.trial{itrial} = transpose(nt_detrend(data.trial{itrial}', 2));
end

%% downsample
cfg=[];
cfg.resamplefs = 250;
data = ft_resampledata(cfg, data);

%% MSE

cfg = [];
cfg.timwin  = 0.5; 
cfg.m        = 2; % pattern length
cfg.r        = 0.5; % similarity criterion
cfg.timescales   = 1:40; %1:40; % scale list
cfg.normalized_r = 1;
cfg.coarsegrainmethod = 'filtskip';
cfg.toi  = -0.75:0.05:1;
cfg.channel = 'MEG';

mse = ft_entropyanalysis(cfg, data)
