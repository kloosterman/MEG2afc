function MEG2afc_concat_runs( subjdir, analysistype, trigger )
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

if ismac
    PRE = '/Users/kloosterman/gridmaster2012/kloosterman';
else
    PRE = '/home/mpib/kloosterman/';
end

cd('~/qsub');

[~,SUBJ] = fileparts(subjdir);

sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
sesdir_codes = ['B' 'D' 'A' 'C'];
% analysistype = 'low';
baseline='trial';
AVG = 'totalpow';

fprintf('%s; %s; %s %s\n', trigger, analysistype, baseline, AVG)

examplefreqpath = fullfile(subjdir, sesdirs{1});
w = what(examplefreqpath);
load(fullfile(examplefreqpath, w.mat{1}))
freq.powspctrm = single(freq.powspctrm);

% set up arrays for analysis
%--------------------------------------------------------------------------
CUTLO  = 0;  CUTHI  = 200;       % cutoff for reading in data
switch trigger;
    case 'stim';
                basetind = find((freq.time>= -0.25) & (freq.time<=0)); %TODO: make -0.2 to 0!!
        %         TIMLO = -0.5; TIMHI = 2.5; % stim
                TIMLO = -0.25; TIMHI = 1.5; % stim
        %         TIMLO          = -0.2;    TIMHI          =  0.4; % for TFR's
%         TIMLO          = -0.2;    TIMHI          =  1.7; % for ramping analysis
    case 'resp';
        %         TIMLO = -1.5; TIMHI = 0.3; % resp
        %         TIMLO = -0.5; TIMHI = 0.4; % resp
        %         TIMLO          = -0.4 ;        TIMHI          =  0.2; % for TFR's
%         TIMLO          = -1.5 ;        TIMHI          =  0.2; % for ramping analysis
        TIMLO          = -1.5 ;        TIMHI          = 0.3; % 
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
% % subjlatrvar = nan( 2, 2,2,length(frind)+3,length(taxis),2,2, 'single' ); % var/mean isoi iwrt freq time drug regime
% subjlatrvar = nan( 2, 2,2,length(frind)+3,length(taxis),2,2,3 ); % var/mean isoi iwrt freq time drug regime diff

for ises = 1:4
    fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
    fprintf('Concatenating runs Subject %s Session %d: %s: %s  . . .\n', SUBJ, ises, sesdir_codes(ises), sesdirs{ises})
    
    if ises==1, ipharm = 1;   imotor=1;
    elseif ises==2, ipharm = 1;   imotor=2;
    elseif ises==3, ipharm = 2;   imotor=1;
    elseif ises==4, ipharm = 2;   imotor=2;
    end

    sesdir = fullfile(subjdir, sesdirs{ises});
%     runlist = dir([sesdir '/*_freq.mat']);
    runlist = dir(fullfile(sesdir, '*_freq.mat'));
    if isempty(runlist);     continue;            end
    cd(sesdir)
    
    powdat = [];    trialinfo = [];
    fprintf('Loading freqs and concatenating . . .\n')
    
    % import pupil data single trials
    
    path = fullfile(PRE, 'projectdata/MEG2afc/pupil/within_subjects');
    filename = sprintf('%s_%s_pupil_data.csv',  lower(SUBJ), sesdir_codes(ises));
    pupil = MEG2afc_import_pupil_file(fullfile(path, filename));
    
    for irun = 1:length(runlist)
        fprintf('Loading run %d: %s . . .\n', irun, runlist(irun).name)
        load(runlist(irun).name); % subj pharma diff pupil SDT
        powdat = [powdat; freq.powspctrm(:,:,:,tind)]; % convert to single after normalization
        
        %         pupil_temp = pupil([pupil.run_nr] == ft_findcfg(freq.cfg, 'runcfg.batch.exp')); % get exp nr instead of irun
        pupil_temp = pupil([pupil.run_nr] == irun); % get exp nr instead of irun
        % throw out weird small RT (<0.01) trials in pupil data which are not in the MEG data and update
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
        freq.trialinfo = [freq.trialinfo [pupil_intersect.decision]']; % add pupil to trialinfo(:,9): baseline decision or feedback
        %save freq.trialinfo for JW behavior analysis
        trl_path = fullfile(PRE, 'projectdata/MEG2afc/pupil/trl_no_artf');
        trl_outfile = sprintf('%s_%s_run%d_trl_data.csv', lower(SUBJ), sesdir_codes(ises), irun);
        csvwrite(fullfile(trl_path, trl_outfile), freq.trialinfo)
        
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
    
    % add correct: correct(1) incorrect (0)
    if  imotor == 1 % ipsi
        error_trls = trialinfo(:,2) ~= trialinfo(:,3); % stim and resp on same side
    else
        error_trls = trialinfo(:,2) == trialinfo(:,3); % on opposite side
    end
    trialinfo(:,12) = error_trls + 1; %1 correct, 2 error
    
    ntrl = size(powdat,1);
    clear freq
    
    %compute basespec within session
    if ismac
        PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/stim';
    else
        PREIN = '/home/mpib/kloosterman/projectdata/MEG2afc/freq/stim';
    end
    filesave = fullfile(fileparts(PREIN), 'low', 'stim', 'respavg', sprintf('basespec_%s_ses%d.mat', SUBJ, ises));
    if strcmp(trigger, 'stim')
        fprintf('Computing baseline spectrum . . .\n')
        basespec_trials = nanmean(powdat(:,:,:,basetind),4);     %average over basetind timebins
        basespec = squeeze(nanmean(basespec_trials, 1));
        %         basespec = squeeze(nanmean(powdat));
        %         basespec = nanmean(basespec(:,:,basetind),3);     %average over basetind timebins
        fprintf('Saving %s . . . \n', filesave)
        save(filesave, 'basespec', 'basespec_trials') %
    else
        fprintf('Loading %s . . . \n', filesave)
        load(filesave) 
    end
    
    %normalize
    normalize = 1;
    if ~normalize
        disp('No normalization applied!!!')
        respdat = powdat;
    else
        respdat = nan(ntrl,length(chlabel),length(frind),length(taxis), 'single');
%         basespec_temp = nan([2 2 size(basespec)]);
        for idiff = 3 % put 3 to select all trials
            for ipup = 3
                diff_ind = trialinfo(:,1) == idiff;
                if ~any(diff_ind), diff_ind(:)=true; end
                pup_ind = trialinfo(:,10) == ipup;
                if ~any(pup_ind), pup_ind(:)=true; end
                
%                 basespec_within = squeeze(nanmean(basespec_trials(diff_ind & pup_ind,:,:), 1)); % TODO use single trial baseline?
%                 basespec_temp(idiff, ipup,:,:) = basespec_within; % saved below
                
                trialind = find(diff_ind & pup_ind);
                fprintf('Normalizing %d trials . . .\n', length(trialind))
                ft_progress('init', 'etf',     'Please wait...');
                for ich = 1:length(chlabel)
                    ft_progress(ich/length(chlabel), 'Processing channel %d from %d', ich, length(chlabel));
                    basedat = squeeze(repmat(squeeze(basespec(ich,:)), [1,size(taxis)]));
                    %                     basedat_within = squeeze(repmat(squeeze(basespec_within(ich,:)), [1,size(taxis)]));
                    for itrial = trialind'
                        %                         respdat(itrial,ich,:,:) = single( (squeeze(powdat(itrial,ich,:,:)) - basedat_within) ./ basedat );
                        basespec_within = squeeze(basespec_trials(itrial,:,:)); % TODO use single trial baseline
                        basedat_within = squeeze(repmat(squeeze(basespec_within(ich,:)), [1,size(taxis)]));
                        
                        respdat(itrial,ich,:,:) = single( (squeeze(powdat(itrial,ich,:,:)) - basedat_within) ./ basedat );
                    end
                end
                ft_progress('close')
            end
        end
        filesave = fullfile(fileparts(fileparts(subjdir)), 'stim', 'respavg', sprintf('basespec_within_%s_ses%d.mat', SUBJ, ises));
        fprintf('Saving %s . . . \n', filesave)
        save(filesave, 'basespec_temp') %
    end
    % pool sensors into occ and motor
    %     TODO add average over all channels for spectra
    %     TODO keep channels for respavg_topo, dimord: [subj, chan, freq, time,ipharm, imotor, idiff, istim, iresp]
%     NKsensorselection
%     temp=[];
%     temp(:,1,:,:) = mean(respdat(:,occind,  :,:),2);
%     temp(:,2,:,:) = mean(respdat(:,motorind,:,:),2);
%     temp(:,3,:,:) = mean(respdat(:,leftoccind,  :,:),2);
%     temp(:,4,:,:) = mean(respdat(:,rightoccind,:,:),2);
%     temp(:,5,:,:) = mean(respdat(:,leftmotorind,  :,:),2);
%     temp(:,6,:,:) = mean(respdat(:,rightmotorind,:,:),2);
%     respdat = temp;
    
    %     % take gamma effect in LAteralization drug-placebo for taking variance
    %     % across trials
    %     freq = [];
    %     freq.time = taxis;    freq.freq = faxis;    freq.dimord = 'rpt_chan_freq_time'; freq.label= chlabel;
    %     freq.powspctrm = powdat;
    %
    %     cfg=[];
    %     cfg.trials = 'all';cfg.avgovertime = 'yes'; cfg.avgoverfreq = 'yes'; cfg.avgoverchan = 'yes';
    %     cfg.latency     = [-0.2 0.2];
    %     cfg.foilim      = [50 60];
    %     NKsensorselection;  cfg.channel = occ;
    %     freq = ft_selectdata(cfg, freq);
    
    clear temp basespec basespec_trials basespec_temp powdat basedat basedat_within
    
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
    
    % Split up by condition and average over trials
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
                            
%                             subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp, ipup, irt, icor) = nanmean(respdat(trial_inds,:,:,:),1);
                            subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp) = nanmean(respdat(trial_inds,:,:,:),1); % only main conditions
                            
                        end
                    end
                end
            end
        end
    end
    ft_progress('close')
    clear powdat respdat
end % ises
respavgout = fullfile(fileparts(subjdir), 'respavg');
mkdir(respavgout)

% save respavg per subj
filesave = fullfile(respavgout, [SUBJ '_' baseline '.mat']);
fprintf('Saving %s . . . \n', filesave)
% save(filesave, 'subjrespavg')
save(filesave) % save all vars



