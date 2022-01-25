function [trl, event] = sortTrials_MEGhh_2afc(cfg)
% 2AFC: baseline 0.6 to 1.4 s ITI = 0.25 to 0.75 s, baseline

% columns trl matrix 2afc:
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
% 16 response on previous trial
% 17 difficulty on previous trial

% self.p_run_start = 127
% self.p_run_end = 126
% self.p_block_start = 31
% self.p_block_end = 30
% self.p_baseline = 8
% self.p_choice_missed = 2
% self.p_stimulus_left_easy = 16
% self.p_stimulus_left_diff = 17
% self.p_stimulus_right_easy = 18
% self.p_stimulus_right_diff = 19
% self.p_left_error = 32
% self.p_left_correct = 33
% self.p_right_error = 34
% self.p_right_correct = 35

hdr    = cfg.headerfile;
fsr    = cfg.fsample;         % in Hz
% trg    = cfg.trialdef.trg;    % 'stim' or 'resp' or 'baseline'
begtrl = cfg.trialdef.begtim; % in seconds
endtrl = cfg.trialdef.endtim; % in seconds
datatype = cfg.datatype; % MEG or eye
irun   = cfg.irun;
ses    = cfg.ses;

if strcmp(ses, 'B'),     ipharm = 1;   imotor=1; %ipharm: 1 =drug, 2 =plac. motor: 1 = contra, 2 = ipsi
elseif strcmp(ses, 'D'), ipharm = 1;   imotor=2;
elseif strcmp(ses, 'A'), ipharm = 2;   imotor=1;
elseif strcmp(ses, 'C'), ipharm = 2;   imotor=2;
end

% sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
% sesdir_codes = ['B' 'D' 'A' 'C'];

if strcmp(datatype, 'MEG')
  [event]= ft_read_event(hdr); %normal
  
  trgval1 = [event(find(strcmp('UPPT001',{event.type}))).value];
  trgsmp1 = [event(find(strcmp('UPPT001',{event.type}))).sample];
  
  trgval2 = [event(find(strcmp('UPPT002',{event.type}))).value];
  trgsmp2 = [event(find(strcmp('UPPT002',{event.type}))).sample];
  
elseif strcmp(datatype, 'eye')
  hdr = ft_read_header(hdr);  % filename_eye
  msg = hdr.orig.msg; % msg cell struct, find trigger lines there
  
  trigind = cellfun(@(x) not(isempty(x)), (strfind(msg, 'trigger')));
  msg_trig = msg(trigind);  

  trgval1 = [];
  trgsmp1 = [];
  for i = 1:length(msg_trig)
    strtok = tokenize(msg_trig{i});
    if strcmp(strtok{end}, 'None')
      continue
    end
    if i > 1 && str2num( strtok{end} ) == trgval1(end)
      continue
    end
    trgval1 = [trgval1 str2num( strtok{end} )];
    trgsmp1 = [trgsmp1 str2num( strtok{2} )];
  end
%   trgoi = trgval1 ~= 31; % remove block start trigger
%   trgval1 = trgval1(trgoi);
%   trgsmp1 = trgsmp1(trgoi);
  
  trgsmp1 = trgsmp1 - hdr.orig.dat(1,1) + 1; %reference to 1st sample rec data
  
  %   trgsmp2 has the button presses only
  trgsmp2 = trgsmp1(trgval1 > 31 & trgval1 < 36);
  
  event = [];
  for istim = 1:length(trgval1)
    event(istim).type      = 'MSG';
    event(istim).sample    = trgsmp1(istim);  %expressed in samples, the first sample of a recording is 1
    event(istim).value     = trgval1(istim); % number or string
    event(istim).offset    = 0; %expressed in samples
    event(istim).duration  = 1;  %expressed in samples
    %   event(istim).timestamp = expressed in timestamp units, which vary over systems (optional)
  end
%   version = msg( cellfun(@(x) not(isempty(x)), (strfind(msg, 'version'))) );
%   if ~isempty(version)
%     version = version{1}(end);
%     % sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
%     % sesdir_codes = ['B' 'D' 'A' 'C'];
%     if strcmp(version, 'B'),     ipharm = 1;   imotor=1;
%     elseif strcmp(version, 'D'), ipharm = 1;   imotor=2;
%     elseif strcmp(version, 'A'), ipharm = 2;   imotor=1;
%     elseif strcmp(version, 'C'), ipharm = 2;   imotor=2;
%     end
%     trl(:,12) = ipharm;
%     trl(:,13) = imotor;
%   else
%     disp('version not present! TODO Sort this out')
%   end
end

trig_ind_stim = find(trgval1 > 15 & trgval1 < 20); % stim onset triggers
% trig_ind_resp = find(trgval1 > 31 & trgval1 < 36); % stim onset triggers
% trig_ind_conf = find(trgval1 > 63 & trgval1 < 68); % stim onset triggers

trl = zeros(length(trig_ind_stim), 15);
trl(:,1) = trgsmp1(trig_ind_stim)' + begtrl*fsr; % find stim onset samples
trl(:,3) = begtrl*fsr;
trl(:,12) = ipharm;  % drug cond
trl(:,13) = imotor;  % motor cond
trl(:,14) = irun;

% find out if all baseline trig have a stim trig after, sometimes not the case
% at run start etc
trig_ind_base = find(trgval1 == 8); 
valid_ind = trgval1(trig_ind_base+1) > 15 & trgval1(trig_ind_base+1) < 20; % 0 if not a stim right after
basestartsmp = trgsmp1(trig_ind_base(valid_ind));
stimstartsmp = trgsmp1(trig_ind_stim); %trig_ind_stim = find(trgval1 > 15 & trgval1 < 20); % stim onset triggers
trl(:,15) = stimstartsmp - basestartsmp; % baseline duration

ctr=1;
for istim = trig_ind_stim % each stim onset trigger
  
  cur_stim_smp = trgsmp1(istim);
  cur_stim_val = trgval1(istim);
  
  % add stimulus type easy or hard
  if cur_stim_val == 16 || cur_stim_val == 18 % easy
    trl(ctr,4) = 1;
  elseif cur_stim_val == 17 || cur_stim_val == 19 % hard
    trl(ctr,4) = 2;
  end
  
  if cur_stim_val == 16 || cur_stim_val == 17  % signal left or right
    trl(ctr,5) = 1; % left
  elseif cur_stim_val == 18 || cur_stim_val == 19
    trl(ctr,5) = 2; % right
  end
  
  % find response that belongs to this trial
  button_smp2 = trgsmp2(find(trgsmp2 > cur_stim_smp, 1, 'first')); % find first resp after stim onset
  if (button_smp2 - cur_stim_smp) > 3*fsr
    warning('No response found, omission? Skipping trial')
%     trl = trl(1:end-1,:);
    trl = trl([1:ctr-1, (ctr+1):end],:);
    continue
    %     trl(ctr,2) = cur_stim_smp + (3*fsr) + endtrl*fsr;
  else % we have a response
    trl(ctr,2) = button_smp2 + endtrl*fsr;
    trl(ctr,8) = button_smp2 - cur_stim_smp; % RT in col 8
    % find out resp specifics in chan 1
    resp_in_stimchan_ind = find(trgsmp1 >= button_smp2(1), 1, 'first'); % comes right after resp in respchan
    
    resp_val = trgval1(resp_in_stimchan_ind);    % button pressed NOT choice yet: left or right button
    if resp_val == 32 || resp_val == 33
      trl(ctr,6) = 1;
    elseif resp_val == 34 || resp_val == 35
      trl(ctr,6) = 2;
    end
    
    if resp_val == 32 || resp_val == 34
      trl(ctr,7) = 0;
    elseif resp_val == 33 || resp_val == 35
      trl(ctr,7) = 1;
    end
    trl(ctr,11) = length(button_smp2); % N button presses, 1 for 1 and up
  end
    
  ctr = ctr+1;
end

trl(:,10) = 1:size(trl,1); % trial index

button=trl(1:end-1,6);
trl(:,16) =[NaN; button]; % prevresp 
% 1 NaN  e.g.
% 2 1
% 2 2
% 1 2

difficulty=trl(2:end,4); % TODO fix this, still not correct!!!!
trl(:,17) = [NaN; difficulty]; 
