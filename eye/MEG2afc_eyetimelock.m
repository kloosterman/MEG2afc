function timelock = MEG2afc_eyetimelock(cfg)
% check resplocked timelock, better do fixed effects with appenddata?, no good for
% removing ERP per run
PREIN = cfg.PREIN;
SUBJ = cfg.SUBJ;
PREOUT = cfg.PREOUT;
sesname = cfg.sesname;

%%
cd(PREIN)

timelock = {};
runlist = dir(sprintf('%s_%s_run*.mat', SUBJ, sesname)); % put simon together for now
for irun = 1:length(runlist)
  disp('Loading');    disp(runlist(irun).name)
  load(runlist(irun).name, 'data_eye');
  data = data_eye;
  
  disp('compute resplocked')
  temp = data;    data = {};
  data{1} = temp;   clear temp
  cfg=[];
  cfg.offset = -round((data{1}.trialinfo(:,5) / ft_findcfg(data{1}.cfg, 'origfs')) * data{1}.fsample);
  data{2} = ft_redefinetrial(cfg,data{1});
  
  disp 'resample on time axis'
  cfg=[];
  cfg.time = -0.5:0.025:1; % 40 Hz
  cfg.time  = cellfun(@double,repmat({cfg.time},length(data{1}.trial),1),'Un',0);
  data{1} = ft_resampledata(cfg, data{1});
  cfg.time = -0.75:0.025:0.5; % 40 Hz
  cfg.time  = cellfun(@double,repmat({cfg.time},length(data{2}.trial),1),'Un',0);
  data{2} = ft_resampledata(cfg, data{2});
  
  cfg=[];
  cfg.method = 'within';
  for itrig = 1:2
    for idiff = 1:2
      cfg.trials = data{itrig}.trialinfo(:,1) == idiff;
      timelock{itrig, idiff}{irun} = ft_timelockanalysis(cfg, data{itrig});
    end
  end
end

disp('average over runs')
cfg=[];
% cfg.latency = [-0.5 1]
timelock = cellfun(@(x) ft_timelockgrandaverage(cfg, x{:}), timelock); % makes struct array

%%
if ismac
  cfg=[];
  cfg.preproc.demean = 'yes';
  cfg.channel = 'EYE_DIAMETER'; % EYE_BLINKS EYE_DIAMETER EYE_SACCADES EYE_HORIZONTAL EYE_VERTICAL
  ft_databrowser(cfg, timelock(1,1))
end
%%
outfile = fullfile(PREOUT, sprintf('%s_%s_freq.mat', SUBJ, sesname));
fprintf('Saving %s\n', outfile)
save(outfile, 'timelock');

% 
%% tryout normalization and plotting
if ismac
  disp 'stimlocked baseline correction'
  cfg=[];
  cfg.baseline = [-0.25 0];
  cfg.baselinetype = 'relchange';
  freq_blc = arrayfun(@(x) ft_freqbaseline(cfg, x), timelock(1,:,:,:));
  
  disp 'get baseline from stim for resp'
  cfg=[];
  cfg.latency = [-0.25 0];
  cfg.avgovertime = 'yes';
  freq_stimbl = arrayfun(@(x) ft_selectdata(cfg, x), timelock(1,:,:,:));
  
  disp 'repmat baseline matrix over time'
  ntim = length(timelock(2,1,1,1).time);
  powspctrm = arrayfun(@(x) repmat(x.powspctrm, 1, 1, ntim), freq_stimbl, 'uni', false);
  [freq_stimbl.powspctrm] = powspctrm{:};
  [freq_stimbl.time] = timelock(2,:,:,:).time;
  
  disp 'baseline correction resp with stimbaseline'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = '((x1-x2)/x2)*100'; % relchange
  freq_blc(2,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), timelock(2,:,:,:), freq_stimbl); % raw
  
    disp 'drug?placebo contrast'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract';
  freq_blc(:,:,3,:) = arrayfun(@(x,y) ft_math(cfg, x,y), freq_blc(:,:,1,:), freq_blc(:,:,2,:));
  
  load colormap_jetlightgray.mat
  cfg=[];
  cfg.layout = 'CTF275_helmet.mat';
  %   cfg.baseline = [-0.25 0];
  %   cfg.baselinetype = 'relchange';
  cfg.colorbar = 'yes';
  cfg.colormap = cmap;
  cfg.zlim = 'maxabs';
  cfg.hotkeys = 'yes';
  figure
  ft_multiplotTFR(cfg, freq_blc(1)) % timelock{itrig, ifreq, idrug, idiff}
end