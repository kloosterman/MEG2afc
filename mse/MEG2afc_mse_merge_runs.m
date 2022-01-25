function mseavg = MEG2afc_mse_merge_runs()

loadtrialinfo = 0; % to compute behav for each run

if ismac
%   basepath = '/Users/kloosterman/gridmaster2012/kloosterman';
  basepath = '/Users/kloosterman/beegfs/';
else
  basepath = '/home/mpib/kloosterman'; % on the cluster
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sesdirs = {'drug_contra'	'drug_ipsi'	'plac_contra'	'plac_ipsi'};

triggers = {'stim'};
% triggers = {'resp'};
PREIN = fullfile(basepath, 'projectdata/MEG2afc/mse');
cd(PREIN);

scalelim = [1 500];
if strcmp(triggers{1}, 'stim')
  timelim = [-0.5 1];
else
  timelim = [-0.75 2];
end

% SUBJ= {'NK17'};
% ALL Subjects:
SUBJ  = { 
    'NK1' ... has funny driftrate
    'NK2' ...
    'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
         'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

nsub = length(SUBJ);

% load example mse to get nchan etc
% subjdir = fullfile(PREIN, SUBJ(1).name);
cd(PREIN)
runlist = dir('*.mat');
fprintf('Loading Subject %s Session %d: . . .\n', runlist(1).name, 1)

% runname = sprintf('*istim%d_iresp%d_irun%d_mse_%s_%s.mat', 3, 3, 1, triggers{1}, fun2run);
% run = dir(runname);

% if isempty(run)
%     runname = sprintf('*icond%d_istim%d_iresp%d_mse_%s.mat', 1, 1, 1, fun2run);
%     run = dir(runname);
% end

load(runlist(1).name)

scind = mse.timescales >= scalelim(1) & mse.timescales <= scalelim(2);
tind = mse.time >= timelim(1) & mse.time <= timelim(2);

nchan = size(mse.sampen,1);
nscales = size(scind,2);
ntim = size(find(tind),2);
mseavg.dat = [];
% mseavg.dat    = nan(nsub,7,2, nchan,nscales,ntim, 4,4,4); % subj block raw_or_psc    toi channel scale    drug motor diff
mseavg.dat    = nan(nsub,3, nchan,nscales,ntim, 4,4,4); % subj block raw_psc_r    toi channel scale    drug motor diff
mseavg.time = mse.time(tind);
mseavg.timescales = mse.timescales(scind);
if loadtrialinfo
  behavperrun.criterion = nan(nsub, 10, 4); %sub block cond
  behavperrun.dprime = nan(nsub, 10, 4); %sub block cond
end
for itrig = 1:length(triggers)
  for isub = 1:length(SUBJ)
    %         subjdir = fullfile(PREIN, SUBJ(isub).name);
    
    
    %         if strcmp(SUBJ(isub).name, 'NK15')
    %             continue
    %         end
    
    %         fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
    %         cd(subjdir)
    for ises = 1:4
      runlist = dir(sprintf('%s_ses*_%s_run*_%s_data.mat', SUBJ{isub}, sesdirs{ises},  triggers{itrig}));
      
      %             runlist = dir(fullfile(subjdir, ['*' sesdirs{ises} '*'] ));
      if length(runlist) < 2
        fprintf('too few runs\n')
        continue
      end
      
      runsampen = nan(length(runlist), nchan, nscales, ntim);
      run_r = nan(length(runlist), nchan, nscales, ntim);
      
      for irun = 1:length(runlist)
        fprintf('Loading %s: . . .\n', runlist(irun).name)
        load(runlist(irun).name)
        
        if loadtrialinfo % load trialinfo for computing behav per run
          runname = sprintf('*icond%d_istim%d_iresp%d_irun*_mse_%s_%s.mat', 1, istim, iresp, triggers{1}, fun2run);
          run = dir(runname);
          load(run(1).name); % to get the data path
          strtok = tokenize(mse.cfg.inputfile, '/');
          load(fullfile(basepath, strtok{5:end})); % load data
        end
        
        %                 if ises == 1 % add to counter for allocating to dat matrices
        %                     sesadd = irun;
        %                 elseif ises == 2
        %                     sesadd = 3+irun;
        %                 elseif ises == 3
        %                     sesadd = 6+irun;
        %                 end
        
        if ises == 1 % add to counter for allocating to dat matrices
          drugind = 1; motorind = 1;
        elseif ises == 2
          drugind = 1; motorind = 2;
        elseif ises == 3
          drugind = 2; motorind = 1;
        else
          drugind = 2; motorind = 2;
        end
        
        % collect runs
        runsampen(irun,:,:,:) = mse.sampen(:,scind,tind);
        run_r(irun,:,:,:) = mse.r(:,scind,tind);
      end
      % avg over runs
      runsampen = squeeze(mean(runsampen));
      run_r = squeeze(mean(run_r));
      %                 mseavg.dat(isub, irun,1, :,:,:,  drugind, motorind, 3) = mse.sampen(:,scind,tind) ; % chan_scale_time now
      mseavg.dat(isub, 1, :,:,:,  drugind, motorind, 3) = runsampen; % chan_scale_time now
      
      %             basetind = mse.time >= -0.5 & mse.time <= -0.25;
      basetind = mse.time >= -0.25 & mse.time <= 0;
      baseline = squeeze(mean(runsampen(:,:, basetind),3));
      mseavg.baseline = baseline;
      baseline_r = squeeze(mean(mse.r(:,:, basetind),3));
      
      %                 mseavg.dat(isub, irun,2, :,:,:,  drugind, motorind, 3) = ( mse.sampen(:,scind,tind) - baseline) ./ baseline * 100; % chan_scale_time now
      mseavg.dat(isub, 2, :,:,:,  drugind, motorind, 3) = ( runsampen - baseline) ./ baseline * 100; % chan_scale_time now
      mseavg.dat(isub, 3, :,:,:,  drugind, motorind, 3) = run_r;
      mseavg.dat(isub, 4, :,:,:,  drugind, motorind, 3) = ( mse.r(:,scind,tind) - baseline_r) ./ baseline_r * 100; % chan_scale_time now
      
      if loadtrialinfo % load trialinfo for computing behav per run
        trialinfo = data.trialinfo(mse.cfg.trials,:);
        hitrate = sum(trialinfo(:,2) == 1 & trialinfo(:,3) == 1) / ...  % Fig report
          sum(trialinfo(:,2) == 1); % N fig trials
        farate = sum(trialinfo(:,2) == 2 & trialinfo(:,3) == 1) / ...  % Hom report
          sum(trialinfo(:,2) == 2); % N hom trials
        % correct 0 and 1 FA/H rates
        if hitrate == 1;  hitrate = 1 - 1/(2*sum(trialinfo(:,2) == 1));   end
        if farate == 0;   farate  = 1 / (2*sum(trialinfo(:,2) == 2));  end
        if hitrate == 0;  hitrate = 1 / (2*sum(trialinfo(:,2) == 1));   end
        if farate == 1;   farate  = 1 - 1 / (2*sum(trialinfo(:,2) == 2));  end
        
        behavperrun.criterion(isub,sesadd,icond) = -0.5 * (norminv(hitrate) + norminv(farate));
        behavperrun.dprime(isub,sesadd,icond) = norminv(hitrate) - norminv(farate);
        
        %                                 behav.criterion(:,:,4) = behav.criterion(:,:,2) - behav.criterion(:,:,1); % Nomiss -  NOFA
      end
    end
  end
end

mseavg.fsample = mse.fsample;
mseavg.SUBJ = SUBJ;
mseavg.mse_leg = {'raw MSE' 'psc MSE' 'r parameter' 'psc r'}; % only normalized here (2)
% mseavg.dimord = 'subj_block_raw_or_psc_toi_channel_scale_drug motor diff';
mseavg.dimord = 'subj_raw-psc-r_toi_channel_scale_drug_motor_diff';
mseavg.sens = MEG2afc_sensorselection();
mseavg.label = mse.label;
% mseavg.stim_conds = {'Fig' 'Hom' 'allstim' 'Fig - Hom' };
% mseavg.resp_conds = {'Report' 'Noreport' 'Reportcomb' 'Report - Noreport' };
% mseavg.behav_conds = {'Conservative' 'Liberal' 'allbehavconds' 'Lib-Cons'};
% mseavg.sdt_conds = {'Hit' 'Miss' 'Hit+Miss' 'Hit-Miss'; ...
%     'FA'  'CR'  'FA+CR' 'FA-CR' ; ...
%     'Hit+FA' 'Miss+CR' 'allSDT' 'Hit+FA-Miss+CR'; ...
%     'Hit-FA' 'Miss-CR'  'Hit-FA+Miss-CR' 'Hit-FA-Miss-CR'};
mseavg.cfg = mse.cfg;
mseavg.dimordsize = size(mseavg.dat);

mseavg.pharm_conds = {'atomox' 'placebo' 'atxplac_avg' 'atomox - placebo'};
mseavg.motor_conds = {'pressipsi' 'presscontra' 'regimecomb' 'pressipsi - presscontra'};
mseavg.diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};


% % avg over blocks
% mseavg.dat(:,7,:, :,:,:, :,:,:) = nanmean(mseavg.dat,2);

% avg over drugs
mseavg.dat(:,:, :,:,:, 3,:,:) = nanmean(mseavg.dat,6);
% avg over motor
mseavg.dat(:,:, :,:,:, :,3,:) = nanmean(mseavg.dat,7);

% subtract drug-plac
mseavg.dat(:,:, :,:,:, 4,:,:) = mseavg.dat(:,:, :,:,:, 1,:,:) - mseavg.dat(:,:, :,:,:, 2,:,:);

% put behavior in
% behavin = fullfile(basepath, 'projectdata', 'critEEG', 'behavior');
% load(fullfile(behavin, 'behavstruct.mat'));
% mseavg.behavior = behav;
%
% behavperrun.dprime(:,:,4) = behavperrun.dprime(:,:,2) - behavperrun.dprime(:,:,1);
% behavperrun.criterion(:,:,4) = behavperrun.criterion(:,:,2) - behavperrun.criterion(:,:,1);
% mseavg.behavperrun = behavperrun;
%
PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/plots/mse';
mseavg.PREOUT = PREOUT;
