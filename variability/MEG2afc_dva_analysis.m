function dvalock = MEG2afc_dva_analysis( cfg )
%Run Schurger dva analysis
% 18 drugipsi has high norm towards end!
%  To avoid confusion, I would only calculate the DVA time course up to the shortest stimulus duration (i.e., shortest RT)
% minus 100 ms (motor latency) within subject and then calculate a NaN mean across subjects.

PREIN = cfg.PREIN;
SUBJ = cfg.SUBJ;
PREOUT = cfg.PREOUT;
sesname = cfg.sesname;

cd(PREIN)

dvalock = {};
runlist = dir(sprintf('%s_%s_run*.mat', SUBJ, sesname));
alldata = {};
disp 'append data'
for irun = 1:length(runlist)
  disp('Loading');    disp(runlist(irun).name)
  load(runlist(irun).name);
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
  
  alldata{irun} = data;
end
data = ft_appenddata([], alldata{:});
clear alldata

% disp 'remove trls with RT < 10th percentile'
% cfg=[];
% [~, trlind] = sort(data.trialinfo(:,5));
% cfg.trials = trlind(round(prctile(trlind, 10)):end); % figure; histogram( data.trialinfo(:,5))
% data = ft_selectdata(cfg, data);
disp 'remove trls with RT < 0.6 s'
cfg=[];
cfg.trials = data.trialinfo(:,5)/1200 > 0.6;
data = ft_selectdata(cfg, data);

disp 'cut off trials at shortest RT'
cfg=[];
cfg.toilim = [-0.3 min(data.trialinfo(:,5))/1200];
data = ft_redefinetrial(cfg, data);

disp('compute resplocked') % not used
temp = data;    data = {};
data{1} = temp;   clear temp
% cfg=[];
% cfg.offset = -round((data{1}.trialinfo(:,5) / ft_findcfg(data{1}.cfg, 'origfs')) * ft_findcfg(data{1}.cfg, 'resamplefs'));
% data{2} = ft_redefinetrial(cfg,data{1});

% if ismac
%   cfg=[];
%   cfg.layout = 'CTF275.lay';
%   cfg.method = 'trial';
%   ft_rejectvisual(cfg, data{1})
% end

dvatypes = { 'withindva', 'acrossdva'};
for itrig = 1%:2
  for idva = 1:2 %:2 % withintrial or acrosstrials
    for idiff = 1:3 %1:2
      cfg = [];
      cfg.dvatype = dvatypes{idva};
      cfg.timwin = 0.1; % in s for within
      cfg.resamplefs = 250;
      
      if itrig == 1 % stim
        %           cfg.toi = -0.2:1/250:1; %    -0.3:0.05:1;
%         cfg.toilim = [-0.2 1.5]; %[-0.2 0.9]; %    -0.3:0.05:1;
        cfg.toilim = [-0.2 data{itrig}.time{1}(end)-0.1];
        cfg.toilim = [-0.2 0.5];
      elseif itrig == 2 % resp
        cfg.toi = -0.75:0.05:0.3; %TODO
      end
      if idiff < 3
        cfg.trials = data{itrig}.trialinfo(:,1) == idiff;
      end
      
      dvalock{idva, idiff} = ft_dva_analysis(cfg, data{itrig}); % TODO make function
      
      temp = dvalock{idva, idiff};
      if ismac
        figure; subplot(211)
        plot(temp.time, squeeze(temp.avg(1,:)))
        subplot(212)
        plot(temp.time, squeeze(temp.avg(2,:)))
      end
      
    end
  end
end
% end

% disp('average over runs')
% dvalock = cellfun(@(x) ft_timelockgrandaverage([], x{:}), dvalock); % makes struct array

outfile = fullfile(PREOUT, sprintf('%s_%s_dvalock.mat', SUBJ, sesname));
fprintf('Saving %s\n', outfile)
save(outfile, 'dvalock');

%
%% tryout normalization and plotting
if ismac
  disp 'stimlocked baseline correction'
  cfg=[];
  cfg.baseline = [-0.25 0];
  cfg.baselinetype = 'relchange';
  freq_blc = arrayfun(@(x) ft_freqbaseline(cfg, x), freq(1,:,:,:));
  
  disp 'get baseline from stim for resp'
  cfg=[];
  cfg.latency = [-0.25 0];
  cfg.avgovertime = 'yes';
  freq_stimbl = arrayfun(@(x) ft_selectdata(cfg, x), freq(1,:,:,:));
  
  disp 'repmat baseline matrix over time'
  ntim = length(freq(2,1,1,1).time);
  powspctrm = arrayfun(@(x) repmat(x.powspctrm, 1, 1, ntim), freq_stimbl, 'uni', false);
  [freq_stimbl.powspctrm] = powspctrm{:};
  [freq_stimbl.time] = freq(2,:,:,:).time;
  
  disp 'baseline correction resp with stimbaseline'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = '((x1-x2)/x2)*100'; % relchange
  freq_blc(2,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), freq(2,:,:,:), freq_stimbl); % raw
  
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
  ft_multiplotTFR(cfg, freq_blc(1)) % freq{itrig, ifreq, idrug, idiff}
end
% end freq function to edit






% %% OLD OLD OLD
% memtic
% cd('~/qsub');
% 
% [~,SUBJ] = fileparts(subjdir);
% 
% sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
% 
% % specify a time axis on%   which you want the data to be resampled.
% fsample_desired = 250;
% % fsample_desired = 50;
% if strcmp(trigger, 'stim')
%   taxis = -1.0:1/fsample_desired:3;
%   basetind = find(taxis >= -0.25 & taxis <= 0);
% else
%   taxis = -1.5:1/fsample_desired:2;
% end
% 
% % 1     2    3      4        5          6       7       8       9       10           11     12
% subjrespavg = nan( 4, 3, length(taxis),2,2,3,3,3, 'single' ); % across var, across norm, within var, within norm
% 
% subjrespavg_accum = nan(2, 3, 2, 2, 3, 3, 3, 'single' );
% 
% for ises = 1:4
%   fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
%   fprintf('Concatenating runs Subject %s Session %d: . . .\n', SUBJ, ises)
%   
%   if ises==1, ipharm = 1;   imotor=1;
%   elseif ises==2, ipharm = 1;   imotor=2;
%   elseif ises==3, ipharm = 2;   imotor=1;
%   elseif ises==4, ipharm = 2;   imotor=2;
%   end
%   
%   sesdir = fullfile(subjdir, sesdirs{ises});
%   runlist = dir(fullfile(sesdir, '*_data.mat'));
%   if isempty(runlist);     continue;            end
%   cd(sesdir)
%   
%   fprintf('Loading data and concatenating . . .\n')
%   
%   concat_data = cell(size(runlist,1),1);
%   for irun = 1:length(runlist)
%     fprintf('Loading run %d\n', irun)
%     load(runlist(irun).name);
%     
%     % TODO interpolate missing channels
%     % compute planar gradient
%     cfg                 = [];
%     cfg.feedback        = 'no';
%     cfg.method          = 'template';
%     if isdeployed
%       cfg.template = fullfile(getenv('HOME'), 'MATLAB/toolbox/fieldtrip-20150803/template/neighbours/ctf275_neighb.mat');
%     else
%       cfg.template = 'ctf275_neighb.mat'; % on 'gridmaster2012',
%     end
%     cfg.neighbours      = ft_prepare_neighbours(cfg, data);
%     
%     cfg.planarmethod    = 'sincos';
%     dataplanar        = ft_megplanar(cfg, data);
%     planar_comp = 1;
%     %         cfg = [];
%     %         dataplanarComb = ft_combineplanar(cfg,dataplanar);
%     
%     %         concat_data{irun} = dataplanarComb;
%     concat_data{irun} = dataplanar;
%   end
%   data = ft_appenddata([], concat_data{:});
%   clear concat_data dataplanarComb dataplanar
%   
%   % select trials with RT between x and y
%   %     minRT = 0.4;
%   minRT = 0.7; % equal amount of samples at each timepoint
%   %     maxRT = 1% 3 sd's
%   cfg=[];
%   cfg.trials = find(data.trialinfo(:,5) > minRT*1200 & ... %bigger than minRT
%     data.trialinfo(:,5) < mean(data.trialinfo(:,5)) + std(data.trialinfo(:,5))*3); % % less than 3 sd's from mean
%   data = ft_selectdata(cfg,data);
%   
%   if strcmp(trigger, 'resp')
%     cfg=[];
%     cfg.offset = -round((data.trialinfo(:,5) / ft_findcfg(data.cfg, 'origfs')) * ft_findcfg(data.cfg, 'resamplefs'));
%     cfg.trials = find(cfg.offset < 1);
%     data = ft_redefinetrial(cfg, data);
%   end
%   
%   time = data.time;
%   [time{:}] = deal(taxis); % defined above
%   cfg=[];
%   cfg.resample = 'yes';
%   cfg.fsample = 500;
%   cfg.time = time; % instead of % cfg.resamplefs = 250;
%   data = ft_resampledataNK(cfg,data);
%   
%   % then, take only trials containing data at all tp
%   cfg=[];
%   cfg.channel = 'MEG';
%   cfg.vartrllength       =  0;
%   cfg.keeptrials         = 'yes';
%   timelock = ft_timelockanalysis(cfg, data);
%   clear data
%   
%   trialinfo = timelock.trialinfo;
%   
%   % make tp's after RT-0.1 go away, to remove resp from across dva
%   % analysis TODO remove trials RT < 1 s
%   timelock_crop = timelock;
%   for itrial = 1:size(trialinfo,1)
%     rt_ind = trialinfo(itrial,5)/1200 - 0.1;
%     if rt_ind > 0
%       if strcmp(trigger, 'stim')
%         timelock_crop.trial(itrial,:, taxis >= rt_ind) = nan; % gum tinds na RT uit
%       else
%         timelock_crop.trial(itrial,:, taxis <= -rt_ind) = nan; % gum tinds na RT uit
%       end
%     end
%     %         timelock_crop.trial(itrial,:, taxis >= rt_ind) = nan; % gum tinds na RT uit
%   end
%   %%
%   % Split up by condition and average over trials
%   % freq.trialinfo   % 1 difficulty condition% 2 signal location% 3 button% 4 correct % 5 RT% 6 empty% 7 trial index% 8 N button presses in trial (omission = 0)
%   
%   NKsensorselection
%   if planar_comp
%     for isoi=1:3 % take both h and v sensors
%       sens.ind{isoi} = [sens.ind{isoi}; sens.ind{isoi} + 268];
%     end
%   end
%   
%   dvaleg = {'across dva', 'across norm', 'within dva', 'within norm'};
%   % 1     2    3      4        5          6       7       8       9       10           11        12
%   % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) RT bin(4) correct (3)
%   fprintf('Split up by condition and average over trials\n')
%   ctr=0;
%   diffs = 1:3;    stims = 1:3;    resps = 1:3;    isois = 1:3;
%   ft_progress('init', 'etf',     'Please wait...');
%   nconds=length(diffs) * length(stims) * length(resps) * length(isois);
%   for isoi=isois % occipital, motor, all
%     for idiff = diffs % easy, hard
%       for istim = stims % 1 = left, 2 = right
%         for iresp = resps % left or right
%           ctr=ctr+1;
%           diff_ind = trialinfo(:,1) == idiff;
%           if ~any(diff_ind), diff_ind(:)=true; end
%           stim_ind = trialinfo(:,2) == istim;
%           if ~any(stim_ind), stim_ind(:)=true; end
%           resp_ind = trialinfo(:,3) == iresp;
%           if ~any(resp_ind), resp_ind(:)=true; end
%           
%           trial_inds = diff_ind & stim_ind & resp_ind;
%           
%           ft_progress(ctr/nconds, 'Processing diff %d stim %d resp %d soi %d: %d trials', idiff, istim, iresp, isoi, ...
%             length(find(trial_inds)));
%           
%           if length(find(trial_inds)) < 1
%             warning('No trials remain for this condition')
%             continue
%           end
%           
%           %select channels
%           chans = sens.ind{isoi};
%           
%           %dva analysis per condition
%           % across-trial variability
%           %                     RT_inds = find(timelock.trialinfo(:,5) > 300); %300 samples = 0.25 s minimum
%           %                     shortest_RT = min(timelock.trialinfo(RT_inds,5))/ 1200;  %  ft_findcfg(data.cfg.previous, 'fsample')                    RT_inds = find(timelock.trialinfo(:,5) > 300); %300 samples = 0.25 s minimum
%           %                     shortest_RT = min(timelock.trialinfo(:,5)) / 1200;  %  ft_findcfg(data.cfg.previous, 'fsample')
%           
%           
%           dvadat = []; mn = [];
%           for it = 1:length(taxis)
%             %                     for it = find(taxis <= shortest_RT - 0.1) % only compute until shortest RT -0.1 (motor latency)
%             vecSet = squeeze(timelock_crop.trial(trial_inds,chans,it))'; %timelock2: RT cropped
%             vecSet = vecSet(:, ~isnan(vecSet(1,:))); % remove nan trials
%             
%             [dvadat(it), mn(it)] = dva(vecSet);
%             
%           end
%           %                             figure; hold on
%           %                             plot(taxis, dvadat(:,:))
%           subjrespavg(1,isoi, 1:it, ipharm, imotor, idiff, istim, iresp) = dvadat;
%           subjrespavg(2,isoi, 1:it, ipharm, imotor, idiff, istim, iresp) = mn;
%           
%           
%           % within-trial variability
%           %                     winSize = 25; % 100 ms @ 250 Hz
%           winSize = fsample_desired  / 10; % 100 ms
%           ntrials = size(timelock.trial,1);
%           dvadat = nan(ntrials, length(timelock.time));
%           L2norm = nan(ntrials, length(timelock.time));
%           %                     dvadat = nan(ntrials, it); % % only compute until shortest RT -0.1 (motor latency)
%           %                     L2norm = nan(ntrials, it);
%           
%           dvadat_accum = nan(ntrials, 1);
%           for itrial= find(trial_inds)'
%             
%             vecSet = squeeze(timelock_crop.trial(itrial,chans,:));
%             % trim nan timepoints in trial
%             halfWin = floor(winSize/2);
%             non_nan_ind = find(~isnan(vecSet(1,:))); % remember non_nan timepoints
%             trl_start = non_nan_ind(1) + halfWin; % remove tp's where sliding window was hanging over the edge
%             trl_end = non_nan_ind(end) - halfWin;
%             vecSet = vecSet(:, ~isnan(vecSet(1,:))); % remove nan time inds in trial
%             
%             [swv, swn] = fast_sw_dva(vecSet,winSize);
%             
%             %                         % take trial RT-0.1 as trial end
%             %                         it = find(taxis <= trialinfo(itrial,5)/1200 - 0.1, 1, 'last');
%             try
%               dvadat(itrial, trl_start:trl_end) = swv(1+halfWin:length(swv)-halfWin); % remove tp's where sliding window was hanging over the edge, go until minimum RT
%             catch
%               disp(itrial)
%             end
%             L2norm(itrial, trl_start:trl_end) = swn(1+halfWin:length(swv)-halfWin);
%             
%             %                         % take shortest RT per session as trial end
%             %                         dvadat(itrial, trl_start:it) = swv(1+halfWin:it); % remove tp's where sliding window was hanging over the edge, go until minimum RT
%             %                         L2norm(itrial, trl_start:it) = swn(1+halfWin:it);
%             
%             %                         dvadat(itrial, trl_start:trl_end) = swv(1+halfWin:length(swv)-halfWin); % remove tp's where sliding window was hanging over the edge
%             %                         L2norm(itrial, trl_start:trl_end) = swn(1+halfWin:length(swv)-halfWin);
%             
%             % for each trial average dva between t = 0 until RT
%             % Update: for each trial average dva between t = 0.25 until 0.75
%             
%             % TOBI: * It is important for us to see these time courses, but ultimately, we want to look at 3 bar graphs:
%             % (i) one for ?raw? baseline DVA; (ii) another for ?raw? trial DVA, (iii) a third for difference between
%             % baseline and raw (your middle row). For trial DVA I would simply collapse across the time window around the dip, i.e., 500 ms after decision onset. This will also give you more robust statistics.
%             %                         if strcmp(trigger, 'stim')
%             %                             RT_in_s = trialinfo(itrial,5)/1250;
%             % %                             ev_accum_inds =   taxis >= 0 & taxis <= RT_in_s;
%             %                             ev_accum_inds =   taxis >= 0.25 & taxis <= 0.75;
%             % %                             dvadat_accum(itrial) = mean(dvadat(itrial, ev_accum_inds), 2);
%             %                             dvadat_accum(itrial) = mean(swv(ev_accum_inds));
%             %                         end
%           end
%           %                             figure; hold on
%           %                             plot(timelock.time, nanmean(dvadat(:,:)))
%           %                     basedat=nanmean(nanmean(dvadat(:,basetind),2));
%           if strcmp(trigger, 'stim')
%             basedat(isoi, ipharm, imotor, idiff, istim, iresp) = nanmean(nanmean(dvadat(:,basetind),2));
%             %                         subjrespavg_accum(1, isoi, ipharm, imotor, idiff, istim, iresp) = nanmean(dvadat_accum);
%             %                         subjrespavg_accum(2, isoi, ipharm, imotor, idiff, istim, iresp) = nanmean(dvadat_accum) -  basedat(isoi, ipharm, imotor, idiff, istim, iresp);
%           else
%             basedatfile = fullfile(fileparts(subjdir), 'dva', 'stim', [SUBJ '_dva.mat']);
%             load(basedatfile, 'basedat')
%           end
%           subjrespavg(3,isoi, 1:size(dvadat,2), ipharm, imotor, idiff, istim, iresp) = nanmean(dvadat); % within raw
%           subjrespavg(4,isoi, 1:size(dvadat,2), ipharm, imotor, idiff, istim, iresp) = nanmean(L2norm);
%           
%           subjrespavg(5,isoi, 1:size(dvadat,2), ipharm, imotor, idiff, istim, iresp) = nanmean(dvadat) - basedat(isoi, ipharm, imotor, idiff, istim, iresp);
%           
%           % TODO save single trial data for mean matching before
%           % computing norm and dva
%           
%         end
%       end
%     end
%   end
%   ft_progress('close')
%   clear data timelock L2norm dvadat vecset timelock_crop
% end % ises
% %
% % figure; hold on; plot(taxis, subjrespavg(3,:,1,1,3,3,3))
% % plot(taxis, subjrespavg(3,:,2,1,3,3,3))
% % figure; hold on; plot(taxis, subjrespavg(1,:,1,1,3,3,3))
% % plot(taxis, subjrespavg(1,:,2,1,3,3,3))
% 
% % TODO meanmatching, compute dva
% 
% respavgout = fullfile(fileparts(subjdir), 'dva', trigger);
% mkdir(respavgout)
% 
% % save respavg per subj
% filesave = fullfile(respavgout, [SUBJ '_dva.mat']);
% fprintf('Saving %s . . . \n', filesave)
% % save(filesave, 'subjrespavg')
% save(filesave) % save all vars
% 
% memtoc
% 
