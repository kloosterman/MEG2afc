function [ subjrespavg ] = MEG2afc_concat_runs( subjdir, trigger )
%Concatenate runs for each subject and save subjrespavg
%   Detailed explanation goes here

bandoi= [12 30; 40 80];


addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
ft_defaults
% rmpath(fileparts(which('nanmean'))) % want to use matlab's own nanmean b/c 'single' support

[~,SUBJ] = fileparts(subjdir);

sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj

analysistype = 'full';
baseline='trial';
AVG = 'totalpow';

fprintf('%s; %s; %s %s\n', trigger, analysistype, baseline, AVG)

PREIN = fullfile('/mnt/homes/home022/nkloost1/projectdata/2afc/freq/', trigger);

examplefreqpath = fullfile(subjdir, sesdirs{1});
w = what(examplefreqpath);
load(fullfile(examplefreqpath, w.mat{end}))

freq.powspctrm = single(freq.powspctrm);

respavgout = fullfile(PREIN, 'respavg');
mkdir(respavgout)

% set up arrays for analysis
%--------------------------------------------------------------------------
CUTLO  = 0;  CUTHI  = 200;       % cutoff for reading in data
switch trigger;
    case 'stim';
        basetind = find((freq.time>=-0.25) & (freq.time<=0));
        TIMLO = -0.5; TIMHI = 2.5; % resp
    case 'resp';
        basetind = find((freq.time>=-1.5) & (freq.time<=-1));
        TIMLO = -1.5; TIMHI = 0.3; % resp
end

freq.time=round(freq.time*100)/100; %get rid of tiny differences in time axis
frind = find((freq.freq >= CUTLO) & (freq.freq <= CUTHI));
faxis = freq.freq(frind);
tind = find((freq.time>=TIMLO) & (freq.time<=TIMHI));
taxis = freq.time(tind);
chlabel = freq.label;

            % 1     2    3      4        5          6       7       8       9       10   
            % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) 
subjrespavg = nan( length(chlabel),length(frind),length(taxis),2,2,3,3,3,1, 'single' );

for ises = 1:4
    fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
    fprintf('Concatenating runs Subject %s Session %d: %s  . . .\n', SUBJ, ises, sesdirs{ises})

    if ises==1, ipharm = 1;   imotor=1;
    elseif ises==2, ipharm = 1;   imotor=2;
    elseif ises==3, ipharm = 2;   imotor=1;
    elseif ises==4, ipharm = 2;   imotor=2;
    end
    
    sesdir = fullfile(subjdir, sesdirs{ises});
    
    runlist = dir([sesdir '/*.mat']);
    if isempty(runlist);     continue;            end
    
    cd(sesdir)
    
    powdat = [];
    trialinfo = [];
    %concat runs
    fprintf('Loading freqs and concatenating . . .\n')
    for irun = 1:length(runlist)
        fprintf('Loading run %d . . .\n', irun)
        load(runlist(irun).name); % subj pharma diff pupil SDT
%         powdat = [powdat; single(freq.powspctrm)];
        powdat = [powdat; freq.powspctrm]; % convert to single after normalization
        trialinfo = [trialinfo; freq.trialinfo];
    end
    ntrl = size(powdat,1);
    clear freq
    
    filesave = fullfile(fileparts(PREIN), 'stim', 'respavg', sprintf('basespec_%s_ses%d.mat', SUBJ, ises));
    if strcmp(trigger, 'stim')
        %compute basespec within session
        fprintf('Computing baseline spectrum . . .\n')
        if length(size(powdat)) == 4 && size(powdat,1) > 1
            basespec = squeeze(nanmean(powdat));
        else
            basespec = squeeze(powdat); %if only 1 trial
        end
        basespec = nanmean(basespec(:,:,basetind),3);     %average over basetind timebins
        fprintf('Saving %s . . . \n', filesave)
        save(filesave, 'basespec') %
    else
        fprintf('Loading %s . . . \n', filesave)
        load(filesave) % save all vars
    end
    
    %normalize
    respdat = nan(ntrl,length(chlabel),length(frind),length(taxis));
    fprintf('Normalizing trials . . .\n')
    ft_progress('init', 'etf',     'Please wait...');
    for ich = 1:length(chlabel)
        ft_progress(ich/length(chlabel), 'Processing channel %d from %d', ich, length(chlabel));
        basedat = squeeze(repmat(squeeze(basespec(ich,:)),1,size(taxis)));
        for itrial = 1:ntrl
            respdat(itrial,ich,:,:) = (squeeze(powdat(itrial,ich,:,:)) - basedat) ./ basedat;
        end
    end
    ft_progress('close')
    clear powdat basespec
    respdat = single(respdat); % now it is safe to convert to single
    
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
    
    
    %% Split up by condition and average over trials
    fprintf('Split up by condition and average over trials\n')
    % freq.trialinfo   % 1 difficulty condition% 2 signal location% 3 button% 4 correct % 5 RT% 6 empty% 7 trial index% 8 N button presses in trial (omission = 0)
    
    % 1     2    3      4        5          6       7       8       9       10
    % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
    itype=0;
    for idiff=1:2 % easy, hard
        for istim = 1:2 % 1 = left, 2 = right
            for ichoice = 1:2 % left or right
                fprintf('.')
                trialind = find( trialinfo(:,2) == istim &  trialinfo(:,1) == idiff & trialinfo(:,3) == ichoice );
                subjrespavg(:,:,:, ipharm, imotor, idiff, istim, ichoice, 1) = nanvar(respdat(trialind,:,:,:)); %in resp 3
            end
        end
    end
    fprintf('\ncollapse over diff, stim and choice . . .\n')
    subjrespavg(:,:,:, ipharm, imotor, 3, 3, 3, 1) = nanvar(respdat);
    
    fprintf('\ncollapse over diff and stim . . .\n')
    for ichoice = 1:2 % left or right
        fprintf('.')
        trialind = find(  trialinfo(:,3) == ichoice );
        subjrespavg(:,:,:, ipharm, imotor, 3, 3, ichoice, 1) = nanvar(respdat(trialind,:,:,:));
    end
    fprintf('\ncollapse over choice and stim . . .\n')
    for idiff = 1:2 % left or right
        fprintf('.')
        trialind = find(  trialinfo(:,1) == idiff );
        subjrespavg(:,:,:, ipharm, imotor, idiff, 3, 3, 1) = nanvar(respdat(trialind,:,:,:));
    end
    fprintf('\ncollapse over diff and choice . . .\n')
    for istim = 1:2 % left or right
        fprintf('.')
        trialind = find(  trialinfo(:,2) == istim );
        subjrespavg(:,:,:, ipharm, imotor, 3, istim, 3, 1) = nanvar(respdat(trialind,:,:,:));
    end
    
    fprintf('\ncollapse over diff . . .\n')
    for istim = 1:2 % left or right
        for ichoice = 1:2 % left or right
            fprintf('.')
            trialind = find(  trialinfo(:,2) == istim &  trialinfo(:,3) == ichoice);
            subjrespavg(:,:,:, ipharm, imotor, 3, istim, ichoice, 1) = nanvar(respdat(trialind,:,:,:));
        end
    end
    fprintf('\ncollapse over stim . . .\n')
    for idiff = 1:2 % left or right
        for ichoice = 1:2 % left or right
            fprintf('.')
            trialind = find(  trialinfo(:,1) == idiff &  trialinfo(:,3) == ichoice);
            subjrespavg(:,:,:, ipharm, imotor, idiff, 3, ichoice, 1) = nanvar(respdat(trialind,:,:,:));
        end
    end
    fprintf('\ncollapse over choice . . .\n')
    for istim = 1:2 % left or right
        for idiff = 1:2 % left or right
            fprintf('.')
            trialind = find(  trialinfo(:,2) == istim &  trialinfo(:,1) == idiff);
            subjrespavg(:,:,:, ipharm, imotor, idiff, istim, 3, 1) = nanvar(respdat(trialind,:,:,:));
        end
    end
    
end % ises
clear respdat

% save respavg per subj
filesave = fullfile(respavgout, [SUBJ '_' baseline 'variance.mat']);
fprintf('Saving %s . . . \n', filesave)
% save(filesave, 'subjrespavg')
save(filesave) % save all vars



