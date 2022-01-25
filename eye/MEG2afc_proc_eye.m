function MEG2afc_proc_eye(cfg)
% preproc edf data
if ismac
  edf2asc = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/eyelink/mac/edf2asc';
  plot = 0;
else
  edf2asc = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/eyelink/linux/edf2asc';
  plot = 0;
end

isub = cfg.isub;
PREIN = cfg.PREIN;
PREOUT = cfg.PREOUT;

cd(PREIN);
runlist = dir(sprintf('nk%d_*', isub));
[~,idx] = sort([runlist.datenum]);
runlist = runlist(idx); % sort by date, hope it is not changed

alldata = {};
for irun = 1:length(runlist)
  [~,runname] = fileparts( runlist(irun).name);
  filename_eye = sprintf('%s.asc', runname);
  if ~exist(filename_eye)
    system(sprintf('%s %s', edf2asc, runlist(irun).name )); %% convert edf to asc
  end
  
  cfg = [];
  cfg.dataset          = filename_eye;
  cfg.montage.tra      = eye(4);
  cfg.montage.labelorg = {'1', '2', '3', '4'};
  cfg.montage.labelnew = {'EYE_TIMESTAMP', 'EYE_HORIZONTAL', 'EYE_VERTICAL', 'EYE_DIAMETER'};
  data_eye = ft_preprocessing(cfg);
  
  %% interpolate blinks
  hdr = ft_read_header(filename_eye); %, 'headerformat', 'eyelink_asc');
  data_eye = interpolate_blinks(hdr, data_eye);
  
  %% low pass filter
  cfg = [];
  cfg.lpfilter = 'yes';
  cfg.lpfreq = 6; % cf de gee 2019 biorxiv
  cfg.lpfiltord = 3;
  cfg.channel = 'EYE_DIAMETER';
  data_eye = ft_preprocessing(cfg, data_eye);

%   %% high pass filter
%   cfg = [];
%   cfg.hpfilter = 'yes';
%   cfg.hpfreq = 0.5;
%   cfg.hpfiltord = 4;
%   cfg.channel = 'EYE_DIAMETER';
%   data_eye = ft_preprocessing(cfg, data_eye);
  
  if plot
    cfg = [];
    cfg.preproc.demean = 'yes';
    cfg.viewmode = 'vertical';
    cfg.channel = 'EYE_DIAMETER';
    ft_databrowser(cfg, data_eye)
  end
  
  %% make trials
  cfg = [];
  cfg.dataset = filename_eye;
  cfg.headerfile = cfg.dataset;
  cfg.fsample = 1000;         % in Hz
  cfg.trialdef.trg = 'stim';    % 'stim' or 'resp' or 'baseline'
  cfg.trialdef.begtim = -1; % in seconds
  cfg.trialdef.endtim = 1; % in seconds
  cfg.datatype = 'eye'; % MEG or eye
  cfg.irun = irun;
  cfg.trialfun = 'sortTrials_MEGhh_2afc';
  
  [cfg] = ft_definetrial(cfg);
  trl = cfg.trl;
  event = cfg.event;
  
  if plot
    cfg=[];
    cfg.event = event;
    cfg.preproc.demean = 'yes';
    cfg.viewmode = 'vertical';
    % cfg.channel = 'EYE_DIAMETER';
    ft_databrowser(cfg, data_eye)
  end
  
  cfg=[];
  cfg.trl = trl;
  data_eye = ft_redefinetrial(cfg, data_eye);
  
%   delete(filename_eye) % remove asc file
  
  alldata = [alldata data_eye];
  clear data_eye
  
  
end

cfg = [];
cfg.appenddim = 'rpt';
cfg.keepsampleinfo = 'no';
data = ft_appenddata(cfg, alldata{:});

%%
cfg = [];
cfg.keeptrials = 'yes';
timelock = ft_timelockanalysis(cfg, data_eye);

cfg=[];
cfg.baseline     = [-0.2 0];
timelock = ft_timelockbaseline(cfg, timelock);

cfg = [];
cfg.keeptrials = 'no';
timelock = ft_timelockanalysis(cfg, timelock);

if plot
cfg = [];
cfg.preproc.demean = 'yes';
cfg.viewmode = 'vertical';
cfg.channel = 'EYE_DIAMETER';
ft_databrowser(cfg, timelock)
end




