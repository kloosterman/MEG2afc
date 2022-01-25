function MEG2afc_merge_freqs( subjdir, analysistype, trigger, phaselock_type )
%Concatenate runs for each subject and save subjrespavg
% columns trialinfo:
% 1 difficulty condition
% 2 signal location
% 3 button
% 4 correct
% 5 RT
% 6 empty
% 7 trial index
% 8 N button presses in trial (omission = 0)
% added here:
% 9 pupil value
% 10 pupil hi lo
% 11 RT bin
% 12 correct

baselinetype = 'respavg_baselinepersession'; % or respavg_baseline_acrosssessions

if ismac
%   PRE = '/Users/kloosterman/gridmaster2012/kloosterman';
  PRE = '/Users/kloosterman/beegfs/projectdata/MEG2afc';
%   PRE = '/Volumes/LNDG/user/Niels/MEG2afc'; % on the server
  
else
%   PRE = '/home/mpib/kloosterman/';
  PRE = '/home/beegfs/kloosterman/projectdata/MEG2afc/';
end

cd('~/qsub');

[~,SUBJ] = fileparts(subjdir);

sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
sesdir_codes = ['B' 'D' 'A' 'C'];
% analysistype = 'low';
baseline='trial';
% phaselock_type = 'totalpow';
% phaselock_type = 'evoked'; % see MEG2afc_merge_freqs_evokedandplot
% phaselock_type = 'induced';

respavgout = fullfile(fileparts(subjdir), baselinetype, phaselock_type);
mkdir(respavgout)
filesave = fullfile(respavgout, [SUBJ '_' baseline '.mat']);

fprintf('%s; %s; %s %s\n', trigger, analysistype, baseline, phaselock_type)

examplefreqpath = fullfile(subjdir, sesdirs{1});
w = what(examplefreqpath);
load(fullfile(examplefreqpath, w.mat{1}))
% freq.powspctrm = single(freq.powspctrm);

% set up arrays for analysis
%--------------------------------------------------------------------------
CUTLO  = 0;  CUTHI  = 200;       % cutoff for reading in data
switch trigger
  case 'stim'
    basetind = find((freq.time>= -0.25) & (freq.time<=0)); %TODO: make -0.2 to 0!!
    TIMLO = -0.5; TIMHI = 0.5; % stim
  case 'resp'
    TIMLO          = -0.5 ;        TIMHI          = 0.3; %maximum
  case 'baseline'
    TIMLO          = -0.3 ;        TIMHI          = -0.3; %maximum
end

freq.time=round(freq.time*100)/100; %get rid of tiny differences in time axis
frind = find((freq.freq>=CUTLO) & (freq.freq<=CUTHI));
faxis = freq.freq(frind);
tind = freq.time>=TIMLO & freq.time<=TIMHI;
taxis = freq.time(tind);
basetind = taxis>= -0.25 & taxis<=0; %TODO: make -0.2 to 0!! only used for stimlocked
chlabel = freq.label;

% 1     2    3      4        5          6       7       8       9       10           11     12
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) RT bin correct
% subjrespavg = nan( length(chlabel),length(frind),length(taxis),2,2,3,3,3,3,4,3, 'single' );
% subjrespavg = nan( 6,length(frind),length(taxis),2,2,3,3,3,3,4,3, 'single' );
subjrespavg = nan( length(chlabel),length(frind),length(taxis),2,2,3,3,3, 'single' ); % only main conditions
subjpowavg = nan( length(chlabel),length(frind),length(taxis),2,2,3,3,3, 'double' ); % only main conditions
% subjBLavg = nan( length(chlabel),length(frind),2,2,3,3,3, 'double' ); % only main conditions, raw pow

% pow_singletrial = nan(2000, 3, 3, 3, 3, 3, 1);  % pow_singletrial(1:length(trial_inds),ipharm, imotor, idiff, istim, iresp, nmeasures)
pow_singletrial = {};  % pow_singletrial(1:length(trial_inds),ipharm, imotor, idiff, istim, iresp, nmeasures)

for ises = 1:4
  fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
  fprintf('Concatenating runs Subject %s Session %d: %s: %s  %s  . . .\n', SUBJ, ises, sesdir_codes(ises), sesdirs{ises}, phaselock_type)
  
  if ises==1, ipharm = 1;   imotor=1;
  elseif ises==2, ipharm = 1;   imotor=2;
  elseif ises==3, ipharm = 2;   imotor=1;
  elseif ises==4, ipharm = 2;   imotor=2;
  end
  
  sesdir = fullfile(subjdir, sesdirs{ises});
  runlist = dir(fullfile(sesdir, sprintf('*%s_freq.mat', phaselock_type)));
  if isempty(runlist);     continue;            end
  cd(sesdir)
  
  powdat = [];    trialinfo = [];
  fprintf('Loading freqs and concatenating . . .\n')
  
  % import pupil data single trials
  %     path = fullfile(PRE, 'projectdata/MEG2afc/pupil/within_subjects');
  path = fullfile(PRE, '/pupil/within_subjects');
  filename = sprintf('%s_%s_pupil_data.csv',  lower(SUBJ), sesdir_codes(ises));
  pupil = MEG2afc_import_pupil_file(fullfile(path, filename));
  
  % remove blink trials detected by Eyelink that were missed by the EOG
  % load EL blinks from csv
  %     path = fullfile(PRE, 'projectdata/MEG2afc/pupil/EyeLink_blink_data');
  path = fullfile(PRE, 'pupil/EyeLink_blink_data');
  filename = sprintf('%s_%s_pupil_data.csv',  lower(SUBJ), sesdir_codes(ises));
  
  blinkdat = MEG2afc_import_blink_file2(fullfile(path, filename));
  
  for irun = 1:length(runlist)
    fprintf('Loading run %d: %s . . .\n', irun, runlist(irun).name)
    load(runlist(irun).name); % subj pharma diff pupil SDT
    
    % remove blink trials detected by Eyelink that were missed by the EOG
    % remove ghost trials of 0.005 s RT, keep trial ind as is
    ghostfreetrials = blinkdat{irun}.rt > 0.01;
    blinkdat{irun}.blinks_nr =  blinkdat{irun}.blinks_nr(ghostfreetrials);
    blinkdat{irun}.rt =  blinkdat{irun}.rt(ghostfreetrials);
    
    %         % set RT's > 3 in blink rt data to 0 as for the MEG data
    %         blinkdat{irun}.rt(blinkdat{irun}.rt > 2.85) = 0; %TODO instead don't take them into account
    
    trial_ind = blinkdat{irun}.trial_nr; % chronological trial #
    blinkfree_trials = trial_ind(blinkdat{irun}.blinks_nr == 0); %select zero blink trials
    
    cfg = []; % % select freq data that is ELapproved
    cfg.trials = ismember(freq.trialinfo(:,7), blinkfree_trials)';
    fprintf('%d additional EyeLink blinks found . . .', size(freq.trialinfo,1)-length(find(cfg.trials)))
    fprintf('%d trials remain. \n', length(find(cfg.trials)))
    freq = ft_selectdata(cfg, freq);
    
    % write pupil baseline from blinkdat into trialinfo
    pupildat = blinkdat{irun}.baseline_pupil(blinkfree_trials); % pupil of blinkfree trials
    clean_trials = ismember(blinkfree_trials , freq.trialinfo(:,7) )'; % cleaned of both blinks and meg artifacts
    freq.trialinfo(:,13) = pupildat(clean_trials);
    
    %check if any trials remain (NK5)
    if size(freq.trialinfo,1) < 2
      warning('No trials remain!')
      continue
    end
    
    %check if rt's match
    trial_rt = blinkdat{irun}.rt; % chronological trial #
    rt_nonzero_ind = trial_rt(freq.trialinfo(:,7)) > 0; % only take nonzero RT meg trials
    r = corr(trial_rt(freq.trialinfo(rt_nonzero_ind,7)), freq.trialinfo(rt_nonzero_ind,5)/1200*1000);
    %         r = corr(trial_rt(freq.trialinfo(:,7)), freq.trialinfo(:,5)/1200*1000);
    if r > 0.98
      disp('RT''s blink and meg data match');
    else
      figure;plot([trial_rt(freq.trialinfo(:,7)), freq.trialinfo(:,5)/1200])
      legend({'meg', 'EL'})
      warning('RT''s blink and meg data DO NOT match')
      %             return
    end
    
    % remove trials with zero RT
    cfg=[];
    cfg.trials = freq.trialinfo(:,5)/1200 > 0;
    freq = ft_selectdata(cfg, freq);
    
    powdat = [powdat; freq.powspctrm(:,:,:,tind)]; % convert to single after normalization
    
    % Add pupil to trialinfo
    pupil_temp = pupil([pupil.run_nr] == irun); % get exp nr instead of irun
    % throw out weird sma+ll RT (<0.01) trials in pupil data which are not in the MEG data and update
    % trial_nr field
    disp('N dropped trials:'); disp(length(find([pupil_temp.rt] < 0.01)))
    pupil_temp = pupil_temp([pupil_temp.rt] > 0.01);
    %         [pupil_temp.trial_nr] = deal({1:length(pupil_temp)});
    for i=1:length(pupil_temp)
      pupil_temp(i).trial_nr = i;
    end
    % set RT's > 3 in pupil data to 0 as for the MEG data
    [pupil_temp([pupil_temp.rt] > 3).rt] = deal(0);
    
    % find pupil trial indices that are still in the meg data
    [~, pupil_ind] = intersect([pupil_temp.trial_nr], freq.trialinfo(:,7)');
    pupil_intersect = pupil_temp(pupil_ind);
    rt_difference = [pupil_intersect.rt]' - (freq.trialinfo(:,5)/1200);
    disp('Mean:')
    disp(mean(rt_difference))
    disp('Corr:')
    disp(corr([pupil_intersect.rt]' , (freq.trialinfo(:,5)/1200)))
    if any(rt_difference > 0.1)
      %         if  corr([pupil_intersect.rt]' , (freq.trialinfo(:,5)/1200)) < 0.99
      %             figure; scatter([pupil_intersect.rt], freq.trialinfo(:,5)'/1200')
      figure; plot([[pupil_intersect.rt]; freq.trialinfo(:,5)'/1200]')
      %             figure; plot(rt_difference)
      %             continue
      fprintf('Subject %s Session %d: %s  . . .\n', SUBJ, ises, sesdirs{ises})
      %             error('Corr < 0.99. Pupil MEG matching went wrong')
      warning('Corr < 0.99. Pupil MEG matching went wrong') % TODO fix pupil meg matching
    end
    freq.trialinfo(:,9) = [pupil_intersect.decision]'; % add pupil to trialinfo(:,9): baseline decision or feedback
    %save freq.trialinfo for JW behavior analysis
    %         trl_path = fullfile(PRE, 'projectdata/MEG2afc/pupil/trl_no_artf');
    trl_path = fullfile(PRE, 'pupil/trl_no_artf');
    trl_outfile = sprintf('%s_%s_run%d_trl_data.csv', lower(SUBJ), sesdir_codes(ises), irun);
    csvwrite(fullfile(trl_path, trl_outfile), freq.trialinfo)
    
    freq.trialinfo(:,6) = irun; % keep track of runno in col 6 (empty before)
    trialinfo = [trialinfo; freq.trialinfo];
  end
  % add pupil hi lo to trial info
  [~, pup_sort_ind] = sort(trialinfo(:,9), 'descend'); %trial inds: hi then low
  prctl_borders = round(prctile(1:length(pup_sort_ind), [40 60])); % use lowest and highest 40%
  trialinfo(pup_sort_ind(1:prctl_borders(1)), 10) = 1; % hi pupil
  trialinfo(pup_sort_ind(prctl_borders(2):end), 10) = 2; % low, mid range stays 0
  
  % add RT bin: slow medium fast
  [~, pup_sort_ind] = sort(trialinfo(:,5), 'descend'); %trial inds: hi then low
  prctl_borders = round(prctile(1:length(pup_sort_ind), [33 66])); % use lowest and highest 40%
  trialinfo(pup_sort_ind(1:prctl_borders(1)), 11) = 1; % slow rt
  trialinfo(pup_sort_ind(prctl_borders(1):prctl_borders(2)), 11) = 2; % medium
  trialinfo(pup_sort_ind(prctl_borders(2):end), 11) = 3; % high
  
  %     % add correct: correct(1) incorrect (0) ALREADY in col4!
  %     if  imotor == 1 % ipsi
  %         error_trls = trialinfo(:,2) ~= trialinfo(:,3); % stim and resp on same side
  %     else
  %         error_trls = trialinfo(:,2) == trialinfo(:,3); % on opposite side
  %     end
  %     trialinfo(:,12) = error_trls + 1; %1 correct, 2 error
  
  % add choice: left vs right
  if  imotor == 1 % ipsi
    trialinfo(:,12) = trialinfo(:,3) - 1; % 0 is choice left
  else
    trialinfo(:,12) = mod(trialinfo(:,3),2); % 0 is choice left
  end
  
  ntrl = size(powdat,1);
  clear freq
  
  %compute basespec within session
  if ismac
%     PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/stim';
    PREIN = '/Users/kloosterman/beegfs/projectdata/MEG2afc/freq/stim';
  else
    PREIN = '/home/beegfs/kloosterman/projectdata/MEG2afc/freq/stim';
  end
  filesave2 = fullfile(fileparts(PREIN), analysistype, 'stim', baselinetype, phaselock_type, sprintf('basespec_%s_ses%d.mat', SUBJ, ises));
  if strcmp(trigger, 'stim')
    fprintf('Computing baseline spectrum . . .\n')
    basespec_trials = nanmean(powdat(:,:,:,basetind),4);     %average over basetind timebins
    basespec = squeeze(nanmean(basespec_trials, 1));
    %         basespec = squeeze(nanmean(powdat));
    %         basespec = nanmean(basespec(:,:,basetind),3);     %average over basetind timebins
    fprintf('Saving %s . . . \n', filesave2)
    save(filesave2, 'basespec', 'basespec_trials') %
  elseif strcmp(trigger, 'baseline')
    % no baseline correction for finelow baseline data
  else
    fprintf('Loading %s . . . \n', filesave2)
    load(filesave2)
  end
  
  %     % load all 4 basespecs and take avg for all the trials
  %     % WARNING make sure that all 4 ses basespecs are already there!
  %     basespec = nan([4,size(basespec)]);
  %     for ises2 = 1:4
  %         fileload = fullfile(fileparts(PREIN), 'low', 'stim', 'respavg', sprintf('basespec_%s_ses%d.mat', SUBJ, ises2));
  %         if exist(fileload)
  %             temp = load(fileload);
  %         end
  %         basespec(ises2,:,:) = temp.basespec;
  %     end
  %     basespec = squeeze(nanmean(basespec)); % avg over 4 ses basespecs
  %     clear temp
  %
  %normalize
  normalize = 1;
  if ~normalize || strcmp(trigger, 'baseline')
    disp('No normalization applied!!!')
    respdat = powdat;
  else
    respdat = nan(ntrl,length(chlabel),length(frind),length(taxis), 'single');
    %         basespec_temp = nan([2 2 size(basespec)]);
    
    %                 basespec_within = squeeze(nanmean(basespec_trials(diff_ind & pup_ind,:,:), 1)); % TODO use single trial baseline?
    %                 basespec_temp(idiff, ipup,:,:) = basespec_within; % saved below
    
    ntrials = size(powdat,1);
    fprintf('Normalizing %d trials . . .\n', ntrials)
    ft_progress('init', 'etf',     'Please wait...');
    for ich = 1:length(chlabel)
      ft_progress(ich/length(chlabel), 'Processing channel %d from %d', ich, length(chlabel));
      basedat = squeeze(repmat(squeeze(basespec(ich,:)), [1,size(taxis)]));
      %                     basedat_within = squeeze(repmat(squeeze(basespec_within(ich,:)), [1,size(taxis)]));
      for itrial = 1:ntrials
        %                         respdat(itrial,ich,:,:) = single( (squeeze(powdat(itrial,ich,:,:)) - basedat_within) ./ basedat );
        %                 basespec_within = squeeze(basespec_trials(itrial,:,:)); % TODO use single trial baseline
        %                 basedat_within = squeeze(repmat(squeeze(basespec_within(ich,:)), [1,size(taxis)]));
        
        basedat_within = squeeze(repmat(squeeze(basespec_trials(itrial,ich,:)), [1,size(taxis)])); % single trial baseline
        respdat(itrial,ich,:,:) = single( (squeeze(powdat(itrial,ich,:,:)) - basedat_within) ./ basedat );
      end
    end
    ft_progress('close')
    %         filesave = fullfile(fileparts(fileparts(subjdir)), 'stim', 'respavg', sprintf('basespec_within_%s_ses%d.mat', SUBJ, ises));
    %         fprintf('Saving %s . . . \n', filesave)
    %         save(filesave, 'basespec_temp') %
  end
  
  % take TF roi and export values
  %     % beta1: sig cluster drug-plac
  %     freq = [];
  %     freq.time = taxis;    freq.freq = faxis;    freq.dimord = 'rpt_chan_freq_time'; freq.label= chlabel;
  %     freq.powspctrm = respdat;
  %
  %     cfg=[];
  %     cfg.trials = 'all';     cfg.avgovertime = 'no';    cfg.avgoverfreq = 'yes';
  %     cfg.avgoverchan = 'yes';
  %     cfg.latency     = [-0.5 0]; % resp-locked
  %     cfg.frequency      = [13 17];
  %     NKsensorselection;  cfg.channel = sens.ind{1}; % occipital
  %     freq1 = ft_selectdata(cfg, freq);
  %     freq1.powspctrm = squeeze(nanmean(freq1.powspctrm,4)); % nanmean for time dim, var length trials
  %
  %     % TF roi 2: corr with driftrate
  %     cfg.latency     = [-0.9 -0.25]; % resp-locked
  %     cfg.frequency      = [19 25];
  %     freq2 = ft_selectdata(cfg, freq);
  %     freq2.powspctrm = squeeze(nanmean(freq2.powspctrm,4)); % nanmean for time dim, var length trials
  
  freq = [];
  freq.time = taxis;
  freq.freq = faxis;
  freq.dimord = 'rpt_chan_freq_time';
  freq.label= chlabel;
  cfg=[];
  cfg.trials = 'all';
  cfg.avgoverfreq = 'yes';
  cfg.avgoverchan = 'yes';
  sens = MEG2afc_sensorselection();
  % raw power baseline
  switch analysistype
    case 'lowfine'
      freq.powspctrm = powdat;
      cfg.frequency      = [8 13];
      cfg.channel = sens.ind{1}; % occipital pooling based on gamma mod
      freq_alpha = ft_selectdata(cfg, freq);
      
      cfg.frequency      = [13 18];
      freq_beta = ft_selectdata(cfg, freq);  % occipital pooling based on gamma mod
      
      cfg.channel = sens.ind{6}; % occipital sensors based on drug effect
      cfg.frequency      = [17 30];
      freq_beta2 = ft_selectdata(cfg, freq); % occipital high beta
    case 'full'
      freq.powspctrm = respdat;
      cfg.frequency      = [40 100];
      cfg.latency        = [0.15 0.4];
      cfg.avgovertime = 'yes';
      cfg.channel = sens.ind{1}; % occpital
      freq_gamma = ft_selectdata(cfg, freq); % occipital gamma
  end
  
  %%     % export csv's for HDDM and single trial analyses
  %     meg_path = fullfile(PRE, 'projectdata/MEG2afc/export/', sprintf('%s', trigger));
  meg_path = fullfile(PRE, 'projectdata/MEG2afc/export/', sprintf('%s', 'baseline'));
  mkdir(meg_path)
  trl_outfile = sprintf('%s_%s_run%d_trl_freq.csv', lower(SUBJ), sesdir_codes(ises), irun);
  switch analysistype % write the whole thing for lowfine
    case 'lowfine'
      
      % columns: subj_idx, correct, RT, session_nr, run_nr, drug
      SUBJlist  = {     'NK1'     'NK2' 'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15' 'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'    };
      outmat = zeros(size(trialinfo,1),10);
      outmat(:,1) = find(strcmp(SUBJlist, SUBJ)) - 1;
      outmat(:,2) = trialinfo(:,4); %0 incorrect, 1 correct
      outmat(:,3) = trialinfo(:,5) / 1200; %RT
      outmat(:,4) = ises; % not chronological anymore!!
      outmat(:,5) = trialinfo(:,6);
      outmat(:,6) = ~isempty(strfind(sesdirs{ises}, 'drug')); % drug is 1
      outmat(:,7) = ~isempty(strfind(sesdirs{ises}, 'contra'));
      outmat(:,8) = trialinfo(:,12); % choice, 0 left
      outmat(:,9) = trialinfo(:,1) - 1; % difficulty, 0 easy, 1 hard
      outmat(:,10) = trialinfo(:,13); % raw baseline pupil
      %             outmat = [outmat freq_alpha.powspctrm freq_beta.powspctrm freq_beta2.powspctrm];
      outmat(:,11) = freq_alpha.powspctrm; % better to index explicitly
      outmat(:,12) = freq_beta.powspctrm;
      outmat(:,13) = freq_beta2.powspctrm;
      %             % make bins based on baseline pupil, percentiles
      %             bin_ranges = 0:20:100; % 10 bins
      %             bin_edges = prctile(trialinfo(:,13), bin_ranges); % binedges raw baseline pupil
      
      % OLD
      %             % bin based on max and min: not all brain states sampled equally often,
      %             minpupil = min(trialinfo(:,13));
      %             maxpupil = max(trialinfo(:,13));
      %             pupilrange = maxpupil-minpupil;
      %             nbins = 5;
      %             bin_edges = minpupil:pupilrange/nbins:maxpupil; % binedges raw baseline pupil
      
      % 5 pupil bins: 30%, 15%, 10%, 15%, 30%
      range = linspace(min(trialinfo(:,13)), max(trialinfo(:,13)), 100);
      %             bin_edges = prctile(range, [0 30 45 55 70 100]);
      bin_edges = prctile(range, [0 40 60 100]);
      
      outmat(:,14) = discretize(trialinfo(:,13), bin_edges, 0:length(bin_edges)-2 );
      
      % make bins based on occipital alpha
      %             bin_ranges = 0:20:100; % 5 bins
      %             bin_edges = prctile(outmat(:,11), bin_ranges); % outmat(:,11) = occ alpha
      freq_alpha.powspctrm = log(freq_alpha.powspctrm); % take log first to make normal
      %             %OLD:
      %             minalpha = min(freq_alpha.powspctrm);
      %             maxalpha = max(freq_alpha.powspctrm);
      %             alpharange = maxalpha-minalpha;
      %             nbins = 3;
      %             bin_edges = minalpha:alpharange/nbins:maxalpha; % binedges raw baseline pupil
      
      % plot alpha distribution
      %             figure; histogram(freq_alpha.powspctrm, 100); title(filename)
      %             ax=gca;
      %             ax.FontSize =18;
      %             xlabel('log alpha power')
      %             ylabel('Frequency of occurrence')
      %             print out -dpdf
      % 5 alpha bins: 30%, 15%, 10%, 15%, 30%
      range = linspace(min(freq_alpha.powspctrm), max(freq_alpha.powspctrm), 100);
      %             bin_edges = prctile(range, [0 33 47 53 67 100]);
      bin_edges = prctile(range, [0 33 67 100]);
      %             bin_edges = prctile(range, [0 40 60 100]);
      %             bin_edges = prctile(range, [0 30 45 55 70 100]);
      
      outmat(:,15) = discretize(freq_alpha.powspctrm, bin_edges, 0:length(bin_edges)-2 );
      
      %             % remove zero RT trials, NK11, bug???
      %             figure; hist( outmat(:,3), 100)
      %             outmat = outmat(outmat(:,3) > 0.1 ,:);
      
      %             % write header and matrix to file
      %             fid = fopen(fullfile(meg_path, trl_outfile), 'w') ;
      %             fprintf(fid, 'subj_idx,response,rt,session_nr,run_nr,drug,simon,choice,diff,baseline_pupil,front_alpha,front_beta,occ_beta,pupilbin0_4,occ_alphabin\n');
      %             fclose(fid);
      %             dlmwrite(fullfile(meg_path, trl_outfile), outmat, '-append', 'precision', 10)
      
    case 'full' % read in the csv from lowfine and append gamma
      outmat = dlmread(fullfile(meg_path, trl_outfile), ',', 1,0);
      % append gamma
      outmat(:,14) = freq_gamma.powspctrm;
      % write out the whole thing inc gamma
      fid = fopen(fullfile(meg_path, trl_outfile), 'w') ;
      fprintf(fid, 'subj_idx, correct, RT, session_nr, run_nr, drug, simon, choice, diff, baseline_pupil, front_alpha, front_beta, occ_beta, occ_gamma\n');
      fclose(fid);
      dlmwrite(fullfile(meg_path, trl_outfile), outmat, '-append', 'precision', 10)
  end
  
  % keep powdat
  clear temp basespec basespec_trials basespec_temp basedat basedat_within freq
  
  %     %             %% ff plotten
  %     close all
  %     sensind = match_str(chlabel,ft_channelselection(SOI{3},chlabel)); % occ sensors
  %     test = nanmean(respdat(:,sensind,:,:),2);
  %     test = squeeze(nanmean(test));
  %     figure
  %     imagesc(taxis,faxis,test, [-0.05 0.05]);
  %     colorbar
  %     set(gca,'Box','off','XTick',[-1 -0.5 0 0.5 1],...    [-0.5,0,0.5,1,1.5,2,2.5]
  %         'YDir','normal','YTick',[0:20:140],...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
  %         'TickDir','out', 'FontSize', 12);
  
  %     MEG2afc_compute_lateralization(); % done in MEG2afc_load_respavg now!
  
  %% Split up by condition and average over trials
  % freq.trialinfo   % 1 difficulty condition% 2 signal location% 3 button% 4 correct % 5 RT% 6 empty% 7 trial index% 8 N button presses in trial (omission = 0)
  
  % 1     2    3      4        5          6       7       8       9       10           11        12
  % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) RT bin(4) correct (3)
  fprintf('Split up by condition and average over trials\n')
  ctr=0;
  %     diffs = 1:3;    stims = 1:3;    resps = 1:3;    pups = 1:3;    rts = 1:4;    cors = 1:3;
  diffs = 1:3;    stims = 1:3;    resps = 1:3;    pups = 3;    rts = 4;    cors = 3;
  ft_progress('init', 'etf',     'Please wait...');
  nconds=length(diffs) * length(stims) * length(resps) * length(pups) * length(rts) * length(cors);
  for idiff = diffs % easy, hard
    for istim = stims % 1 = left, 2 = right
      for iresp = resps % left or right
        for ipup = pups % hi, lo pupil
          for irt = rts
            for icor = cors
              ctr=ctr+1;
              diff_ind = trialinfo(:,1) == idiff;
              if ~any(diff_ind), diff_ind(:)=true; end
              stim_ind = trialinfo(:,2) == istim;
              if ~any(stim_ind), stim_ind(:)=true; end
              resp_ind = trialinfo(:,3) == iresp;
              if ~any(resp_ind), resp_ind(:)=true; end
              pup_ind = trialinfo(:,10) == ipup;
              if ~any(pup_ind), pup_ind(:)=true; end
              rt_ind = trialinfo(:,11) == irt;
              if ~any(rt_ind), rt_ind(:)=true; end
              cor_ind = trialinfo(:,12) == icor;
              if ~any(cor_ind), cor_ind(:)=true; end
              
              trial_inds = diff_ind & stim_ind & resp_ind & pup_ind & rt_ind & cor_ind;
              
              ft_progress(ctr/nconds, 'Processing diff %d stim %d resp %d pup %d rt %d cor %d: %d trials', idiff, istim, iresp, ipup, irt, icor, ...
                length(find(trial_inds)));
              
              if length(find(trial_inds)) < 1
                warning('No trials remain for this condition')
                continue
              end
              
              %subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp, ipup, irt, icor) = nanmean(respdat(trial_inds,:,:,:),1);
              if ~strcmp(trigger, 'baseline')
                subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp) = nanmean(respdat(trial_inds,:,:,:),1); % only main conditions
              end
              subjpowavg(:,:,:, ipharm, imotor, idiff, istim, iresp) = nanmean(powdat(trial_inds,:,:,:),1); % only main conditions
            end
          end
        end
      end
    end
  end
  ft_progress('close')
  clear respdat powdat
end % ises

% save respavg per subj
fprintf('Saving %s . . . \n', filesave)
save(filesave) % save all vars



