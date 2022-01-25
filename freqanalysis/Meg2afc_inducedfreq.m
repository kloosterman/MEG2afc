close all
if ismac
  basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
end

cd(fullfile(basepath, 'preproc'))

load NK11_drug_contra_run3.mat % ssvep @ 10 Hz?
% load NK21_drug_ipsi_run2.mat

% cfg=[];
% cfg.offset = -round((data.trialinfo(:,5) / ft_findcfg(data.cfg, 'origfs')) * ft_findcfg(data.cfg, 'resamplefs'));
% data = ft_redefinetrial(cfg,data);

disp(' synthetic planar computation')
cfg              = [];
cfg.feedback     = 'no';
cfg.method       = 'template';
if ismac
  cfg.template     = 'ctf275_neighb.mat';
else
  cfg.template     = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
end
cfg.planarmethod = 'sincos';
cfg.channel      = {'MEG'};
cfg.trials       = 'all';
cfg.neighbours   = ft_prepare_neighbours(cfg, data);
data =  ft_megplanar(cfg, data);

data_BLC = remove_ERP_fromdata(data, 'subtract'); % on planar grads?

timelock_BLC = ft_timelockanalysis([], data_BLC);
timelock_BLC = ft_combineplanar([], timelock_BLC);

cfg=[];
cfg.output= 'pow';
cfg.channel= 'MEG';
cfg.keeptapers= 'no';
cfg.pad = 7;
cfg.method= 'mtmconvol';
cfg.trials= 'all';
cfg.toi = -0.5:0.05:0.75;

% % % low frequency-optimized analysis
cfg.taper = 'hanning';
cfg.keeptrials  = 'no'; % needed for fourier-output
cfg.foi = 1:35;
cfg.t_ftimwin = ones(length(cfg.foi),1) .* 1;
cfg.tapsmofrq = ones(length(cfg.foi),1) .* 0;

% % % high frequency-optimized analysis (smooth)
%         cfg.taper = 'dpss';
%         cfg.keeptrials  = 'yes';
%         cfg.foi = 36:2:120;
%         cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
%         cfg.tapsmofrq = ones(length(cfg.foi),1) .* 8;

freq = ft_freqanalysis(cfg, data); % run it
freq = ft_combineplanar([], freq);

freq_erp = ft_freqanalysis(cfg, timelock); % run it
freq_erp = ft_combineplanar([], freq_erp);

cfg=[];
cfg.layout = 'CTF275_helmet.lay';
cfg.baseline = [-0.25 0];
cfg.colorbar = 'yes';
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.zlim = 'maxabs';
cfg.baselinetype = 'relchange';
cfg.hotkeys = 'yes';
figure
ft_multiplotTFR(cfg, freq)

figure
ft_multiplotTFR(cfg, freq_erp)

timelock = ft_combineplanar([], timelock);
cfg=[];
cfg.layout = 'CTF275_helmet.lay';
figure
% ft_multiplotER(cfg, timelock)
ft_multiplotER(cfg, timelock_noBLC, timelock_BLC)


