function [behav] = MEG2afc_readbehavior(cfg)
% get various behavior from trl field
% TRL OUTLINE:
% 1 trial start
% 2 trial end
% 3 trial offset
% 4 difficulty condition
% 5 signal location
% 6 button
% 7 correct
% 8 RT
% 9 empty
% 10 trial index
% 11 N button presses in trial (omission = 0)
% 12 drugcond: 1 drug, 2 placebo
% 13 motorcond: 1 ipsi, 2 contra
% 14 run number
% 15 baseline duration

omitdroppedbymeg = 1;

PREIN = cfg.PREIN;
SUBJ = cfg.SUBJ;
sesname = cfg.sesname;
outfile = cfg.outfile;

behav = [];

cd(PREIN)

runlist = dir(sprintf('%s_%s_run*.mat', SUBJ, sesname));
p_repeatunbalanced = nan(9,1);
p_repeatbalanced = nan(9,2);
button_bias = nan(9,2);
dprime = nan(9,2); % easy hard
propcorrect = nan(9,2); % easy hard
criterion = nan(9,2);
bpm = nan(9,1);
basepupil = nan(9,1);
RT= nan(9,2);
RTsd= nan(9,2);
ntrials = nan(9,2);
ddmmat_runs = [];
trl_runs = [];
datakeep={};
trlkeep = [];
disp 'compute behavior per run'
for irun = 1:length(runlist)+1
  if length(runlist) == 0
    disp('no data')
    continue
  end
  if irun < length(runlist)+1
    disp('Loading');    disp(runlist(irun).name)
    load(runlist(irun).name, 'data_eye'); % only eye data!!
    data_eye.fsample = round(data_eye.fsample); % sometimes slightly off
    datakeep{end+1} = data_eye;
  else
    disp 'compute collapsed over runs'
    data_eye = ft_appenddata([], datakeep{:});
    irun = 9; % always put in 9
%     clear datakeep
  end

  if irun < length(runlist)+1
    trl = ft_findcfg(data_eye.cfg, 'trl');
    trlkeep = [trlkeep; trl];
  else
    trl = trlkeep;
  end

  if omitdroppedbymeg
%     trl = trl(ismember(trl(:,10), data_eye.trialinfo(:,7)),:);
    trl = [zeros(length(data_eye.trialinfo),3) data_eye.trialinfo];
  end
  difficulty = trl(:,4);
  stim = trl(:,5);  % the stim presented to the subject
  button = trl(:,6);  % the button pressed by the subject
  choice = button; % actual choice: opposite for contra session
  if trl(:,13) == 2 % presscontra session
    choice = mod( choice,2 ) + 1;
  end


  % use choice to compute p_repeatbalanced
  p_repeatunbalanced(irun,1) = sum(diff(choice) == 0) / (numel(choice(2:end))); % trials with same choice as previous trial / ntrials
  p_repeatbalanced(irun,1) = sum(diff(choice) == 0 & choice(2:end,:) == 1) / (sum(choice(2:end)==1)); 
  p_repeatbalanced(irun,2) = sum(diff(choice) == 0 & choice(2:end,:) == 2) / (sum(choice(2:end)==2));

%   p_repeatunbalanced(irun,1) = sum(diff(button) == 0) / (numel(button(2:end))); % trials with same button as previous trial / ntrials
%   p_repeatbalanced(irun,1) = sum(diff(button) == 0 & button(2:end,:) == 1) / (sum(button(2:end)==1));
%   p_repeatbalanced(irun,2) = sum(diff(button) == 0 & button(2:end,:) == 2) / (sum(button(2:end)==2));

  button_bias(irun,1) = (sum(button == 1)) / sum( button == 1 | button == 2 ); % TODO work this out
  button_bias(irun,2) = (sum(button == 2)) / sum( button == 1 | button == 2 ); % TODO work this out;
  
  disp 'compute SDT d'' and c'
  for idiff = 1:2
    if contains(sesname, 'contra')
      H = sum(stim == 1 & difficulty == idiff & button == 2) / sum(stim == 1 & difficulty == idiff);
      FA = sum(stim == 2 & difficulty == idiff & button == 2) / sum(stim == 2 & difficulty == idiff);
      propcorrect(irun,idiff) = sum(difficulty == idiff & stim ~= button) / sum(difficulty == idiff);
    else
      H = sum(stim == 1 & difficulty == idiff & button == 1) / sum(stim == 1 & difficulty == idiff);
      FA = sum(stim == 2 & difficulty == idiff & button == 1) / sum(stim == 2 & difficulty == idiff);
      propcorrect(irun,idiff) = sum(difficulty == idiff & stim == button) / sum(difficulty == idiff);
    end
    if H == 1; H = 0.99; end
    if FA == 0; FA = 0.01; end
    dprime(irun,idiff) = norminv(H) - norminv(FA);
    criterion(irun,idiff) = -0.5 * (norminv(H) + norminv(FA)); %means??
%     RTok = zscore(trl(:,8)) < 3; % DROP RT's > 3 sd away 
%     RT(irun,idiff) = median(trl(difficulty == idiff & RTok,8)) / 1200;
%     RTsd(irun,idiff) = std(trl(difficulty == idiff & RTok,8)) / 1200;
    RT(irun,idiff) = mean(trl(difficulty == idiff,8)) / 1200; % no drop RTs
    RTsd(irun,idiff) = std(trl(difficulty == idiff,8)) / 1200;
    ntrials(irun,idiff) = sum( difficulty == idiff);  % size(trl,1);
    
  end
  
%   disp 'find out if blink occurred during trial'
%   for itrial = 1:length(data_eye.trial)
%   data_eye.trial
%   end
  
  disp 'make HDDM csv'
  disp 'start with all 160 trials per run'
  trl = ft_findcfg(data_eye.cfg, 'trl');
  difficulty = trl(:,4);
  stim = trl(:,5);  % the stim presented to the subject
  button = trl(:,6);  % the button pressed by the subject
  choice = button; % actual choice: opposite for contra session
  if trl(:,13) == 2 % presscontra session
    choice = mod( choice,2 ) + 1;
  end

  % collect data for history bias and accuracy ddm
  % columns: subj_idx, stimulus, response, prevresp, correct, RT, drug, simon,
  % run_nr, difficulty
  SUBJlist  = {     'NK1'     'NK2' 'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15' 'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'    };
  ddmmat = zeros(size(trl,1),11);
  ddmmat(:,1) = find(strcmp(SUBJlist, SUBJ)) - 1; %subj_idx;
  ddmmat(:,2) = trl(:,5)-1;
  ddmmat(:,3) = button-1;
  ddmmat(:,4) = [NaN; button(1:end-1)-1]; %shifted resp
  ddmmat(:,5) = trl(:,7); %0 incorrect, 1 correct
  ddmmat(:,6) = trl(:,8) / 1200; %RT
  ddmmat(:,7) = contains(sesname, 'drug'); % drug is 1
  ddmmat(:,8) = contains(sesname, 'contra'); % contra is 1
  ddmmat(:,9) = irun;
  ddmmat(:,10) = trl(:,4) - 1; % difficulty, 0 easy, 1 hard
  ddmmat(:,11) = trl(:,10); % trial counter
  ddmmat(:,12) = ismember(trl(:,10), data_eye.trialinfo(:,7)); % 1 = included in MEG
  
  %   minRT = 0.2; % in sec
  %   ddmmat = ddmmat(ddmmat(:,6) > (minRT*1200),:); % drop RT's > 3SD
  %   ddmmat = ddmmat(zscore(ddmmat(:,6))<3,:); % drop RT's > 3SD
  
  %   ddmmat = ddmmat(~isnan(ddmmat(:,4)),:); % remove trials without prevresp
  %   ddmmat = ddmmat(ddmmat(:,3)>-1,:); % remove missed responses
  %   ddmmat = ddmmat(ddmmat(:,4)>-1,:); % remove missed previous resp
  
  if irun < length(runlist)+1
    ddmmat_runs = [ddmmat_runs; ddmmat];
    trl_runs = [trl_runs; trl];
  end
  %   disp 'get heartbeats'
  %   cfg=[];
  %   cfg.continuous = 'no';
  %   cfg.artfctdef.ecg.feedback = 'no';
  %   cfg.artfctdef.ecg.channel = 'EEG059';
  %   cfg.artfctdef.ecg.inspect = 'EEG059'; %Nx1 list of channels which will be shown in a QRS-locked average
  %   [cfg, artifact] = ft_artifact_ecg(cfg, data_eye);
  %   netdur = sum(cellfun(@length, data_eye.time)) / data_eye.fsample / 60;
  %   bpm(irun,1) = length(artifact) / netdur;
  
  %% Get heartrate from epoched or continous data (as done in MEG2afc_preproc)
  getHRfrom ='continous'; %{'epoched' 'continous'};
  switch getHRfrom
    case 'epoched'
      disp 'Detecting heartbeats . . .' % updated to address 1 funny subject
      cfg=[];
      %   cfg.trl = [(data_eye.fsample*3) length(data_eye.trial{1})-(data_eye.fsample*3) 0];
      cfg.continuous = 'no';
      cfg.artfctdef.ecg.channel = {'EEG059'};
      if ismac
        cfg.artfctdef.ecg.feedback = 'no';
      else
        cfg.artfctdef.ecg.feedback = 'no';
      end
      cfg.artfctdef.ecg.inspect = 'EEG059';
      [cfg_heartbeats, heartbeats] = ft_artifact_ecg(cfg, data_eye);
      
      cfg_heartbeats.heartbeats = heartbeats;
      inter_beat_durs = diff(heartbeats(:,1)) / data_eye.fsample; % in sec
      inter_beat_durs = inter_beat_durs(zscore(inter_beat_durs) < 3); % remove outlier durs, indicates ECG is loose
      bpm(irun,1) = length(inter_beat_durs) / (sum(inter_beat_durs)/60);
      cfg_heartbeats.bpm = bpm;
    case 'continous'
      disp 'Loading heartbeats . . .' % updated to address 1 funny subject
      if ismac
        path = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/preproczapline-plus/heartbeats';
      else
        path = '/home/mpib/kloosterman/projectdata/MEG2afc/preproczapline-plus/heartbeats';
      end
      if irun < length(runlist)+1
        load(fullfile(path, [runlist(irun).name(1:end-12) '.mat'])) % cfg_heartbeats comes out
        bpm(irun,1) = cfg_heartbeats.bpm;
      else
        bpm(irun,1) = nanmean(bpm); % avg over runs
      end
      
      
  end
  

  %%
  disp 'get baseline pupil'
  cfg=[];
  cfg.channel = 'EYE_DIAMETER';
  timelock = ft_timelockanalysis(cfg, data_eye);
%     figure; plot(timelock.time, timelock.avg)
  cfg=[];
  cfg.latency = [-0.75 -0.25];
  cfg.avgovertime = 'yes';
  timelock = ft_selectdata(cfg, timelock);
  basepupil(irun,1) = timelock.avg;
end
behav.dprime = dprime;
behav.criterion = criterion;
behav.RT = RT;
behav.RTsd = RTsd;
behav.p_repeatunbalanced = p_repeatunbalanced;
behav.p_repeatbalanced = p_repeatbalanced;
behav.button_bias = button_bias;
behav.ddmmat_runs = ddmmat_runs;
behav.bpm = bpm;
behav.basepupil = basepupil;
behav.trl_runs = trl_runs; % just all trials, nothing removed
behav.ntrials = ntrials;
behav.propcorrect = propcorrect;

disp(outfile);
save(outfile, 'behav');


