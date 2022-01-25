function MEG2afc_mse_analysis(cfg)

fileload = cfg.fileload;
outputfile =     cfg.outputfile;
trigger = cfg.trigger;
removeERP = cfg.removeERP;
removeERPmethod = cfg.removeERPmethod;
SUBJ = cfg.SUBJ;

clear cfg

disp(fileload);
load(fileload); % data comes out

disp(outputfile);

% cfg2 = [];
% cfg2.channel = 'MEG';
% ft_databrowser(cfg2, data)

% downsample
cfg=[];
cfg.resamplefs = 256;
cfg.channel = 'MEG';
data = ft_resampledata(cfg, data);

disp('discard trials with RT < 0.2 and > 2.5 and crop a bit')
cfg=[];
cfg.trials = data.trialinfo(:,5)/1200 > 0.2 & data.trialinfo(:,5)/1200 < 2.5;
data = ft_selectdata(cfg, data);

cfg=[];
cfg.toilim = [-1 3];
data = ft_redefinetrial(cfg,data);

close all
%%
disp('realigning MEG')
cfg=[];
if ismac
  cfg.template       = {'/Users/kloosterman/Dropbox/tardis_code/MATLAB/MEG_HH_analysis/megrealign/ctf275.mat'};
  cfg.headmodel   = fullfile('/Volumes/LNDG/user/Niels/MEG2afc/MRI/NKdet', SUBJ, [SUBJ '_hdm.mat']);
else
  cfg.template       = {'/home/mpib/kloosterman/MATLAB/MEG_HH_analysis/megrealign/ctf275.mat'};
  cfg.headmodel   = fullfile('/home/beegfs/kloosterman/projectdata/MEG2afc/MRI', SUBJ, [SUBJ '_hdm.mat']);
end
cfg.inwardshift = 1;
cfg.feedback = 'no';
data = ft_megrealign(cfg, data);

if strcmp(removeERP, 'yes')
  data = remove_ERP_fromdata(data, removeERPmethod); % TODO add remove type to cfg
  data.cfg.previous = data.cfg; % log detrending done in cfg to keep overview
  data.cfg = keepfields(data.cfg, 'previous');
  data.cfg.funthatwasrun = 'remove_ERP_fromdata';
  data.cfg.method = [removeERPmethod '_percondition'];
  clear data_ERP1
end

if strcmp(trigger, 'resp')
  cfg=[];
  cfg.offset = -data.trialinfo(:,4);
  cfg.trials = find(cfg.offset < 1);
  data = ft_redefinetrial(cfg, data);
end

% convert timelock back to raw
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');

% make sure there are no nans in raw data (e.g. if coming from timelock with var trl lengths)
ntrials = length(data.trial);
cfg = [];
cfg.begsample = nan(ntrials,1);
cfg.endsample = nan(ntrials,1);
for itrial = 1:ntrials
  nonnans = find(~isnan(data.trial{itrial}(1,:)));
  cfg.begsample(itrial,:) = nonnans(1);
  cfg.endsample(itrial,:) = nonnans(end);
end
data = ft_redefinetrial(cfg, data);

% %% synthetic planar computation
% cfg              = [];
% cfg.feedback     = 'no';
% cfg.method       = 'template';
% cfg.planarmethod = 'sincos';
% cfg.channel      = {'MEG'};
% cfg.trials       = 'all';
% cfg.neighbours   = ft_prepare_neighbours(cfg, data);
% data             = ft_megplanar(cfg,data);

%% entropy analysis
cfg = [];
cfg.timwin        = 0.5; % sliding window size
cfg.m             = 2; % pattern length
cfg.r             = 0.5; % similarity criterion
cfg.timescales    = 1:42; %1:40; % scale list
cfg.recompute_r  = 'perscale_toi_sp'; 
cfg.coarsegrainmethod = 'filtskip';  % pointavg
cfg.mem_available = 16e9; % in bytes, 8e9 default
cfg.allowgpu = true;
cfg.toi  = -0.75:0.05:1;
cfg.filtmethod = 'lp';
cfg.outputfile = outputfile;

% cfg.hpfilter      = 'no';
% cfg.match_trlcounts = 'no';
% cfg.timescales    = 42; cfg.toi  = 1;  % for testing
% mse = ft_entropyanalysis_reverseFilterChunking(cfg, data);
mse = ft_entropyanalysis(cfg, data);

freq=[];
freq.label = mse.label;
freq.freq = mse.timescales;
freq.time = mse.time;
freq.powspctrm = mse.sampen;
freq.dimord = 'chan_freq_time';
freq.trialinfo = mse.trialinfo;
freq.cfg = mse.cfg;

cfg               = [];
cfg.method = 'sum';
mse     = ft_combineplanar(cfg, freq);


freq=[];
freq.label = mse_ax.label;
freq.freq = mse_ax.timescales;
freq.time = mse_ax.time;
freq.powspctrm = mse_ax.sampen;
freq.dimord = 'chan_freq_time';
freq.trialinfo = mse_ax.trialinfo;
freq.cfg = mse_ax.cfg;
mse_ax = freq

%%
load colormap_jetlightgray.mat
cfg=[];
cfg.layout = 'CTF275_helmet.lay';
cfg.colormap = cmap;
cfg.zlim = 'maxabs';
cfg.baseline = [-0.5 -0.25];
cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.hotkeys = 'yes';
figure
ft_multiplotTFR(cfg,mse_plan);
figure
ft_multiplotTFR(cfg,mse_ax);
% 


%% old
% %% MSE
% cfg = [];
% cfg.timwin  = 0.5; 
% cfg.m        = 2; % pattern length
% cfg.r        = 0.5; % similarity criterion
% cfg.timescales   = 1:42; %1:40; % scale list
% cfg.normalized_r = 1;
% cfg.coarsegrainmethod = 'filtskip';
% cfg.toi  = -0.75:0.05:1;
% cfg.outputfile = cfg.outputfile;
% cfg.channel = 'MEG';
% 
% mkdir(fileparts(cfg.outputfile))
% mse = ft_entropyanalysis(cfg, data)


