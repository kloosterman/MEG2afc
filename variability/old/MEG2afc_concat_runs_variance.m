function [ subjrespavg, subjrespvar ] = MEG2afc_concat_runs_variance( subjdir, trigger )
%Concatenate runs for each subject, compute variance across trials and save 
%   Detailed explanation goes here

if nargin==0
    subjdir = '/mnt/homes/home022/nkloost1/projectdata/2afc/preproc/NK3';
    trigger = 'stim';
end
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

PREIN = '/mnt/homes/home022/nkloost1/projectdata/2afc/preproc/variance/';
respavgout = fullfile(PREIN, 'respavg');
mkdir(respavgout)

exampledatapath = fullfile(subjdir, sesdirs{1});
w = what(exampledatapath);
load(fullfile(exampledatapath, w.mat{end-1}))

switch trigger;
    case 'stim';
%         basetind = find((freq.time>=-0.25) & (freq.time<=0));
        TIMLO = -0.5; TIMHI = 2.5; 
        BASELO = -0.25;
        BASEHI = 0;
    case 'resp';
%         basetind = find((freq.time>=-1.5) & (freq.time<=-1));
        TIMLO = -1.5; TIMHI = 0.3; % resp
end

% if resp, ft_redefinetrial

cfg_megplanar = [];
cfg_megplanar.method = 'distance';
cfg_megplanar.channel = 'MEG';
cfg_megplanar.neighbours = ft_prepare_neighbours(cfg_megplanar, data);
cfg_megplanar.planarmethod = 'sincos';
data = ft_megplanar(cfg_megplanar, data);


%downsample to 100 Hz, put on fixed axis?
resamplefs = 100;
TIMLO = -0.5;
TIMHI = 1.5;
ntrials = length(data.trial);
timemat = repmat(TIMLO:1/resamplefs:TIMHI,ntrials,1); % specify a time axis on which you want the data to be resampled.

cfg_downsample=[];
cfg_downsample.fsample = data.fsample;
cfg_downsample.time = mat2cell(timemat, ones(ntrials,1), size(timemat,2))';
cfg_downsample.demean = 'no';
cfg_downsample.detrend = 'no';
data = ft_resampledataNK(cfg_downsample,data); %edit: put nan where no data at cfg_downsample.time
data.cfg.trl = []; %remove trl, is useless after resampling

% convert to timelock
cfg_timelock=[];
% cfg_timelock.vartrllength = 2;
cfg_timelock.keeptrials = 'yes';
timelock= ft_timelockanalysis(cfg_timelock,data);

timelock = ft_combineplanar([], timelock);

% set up arrays for analysis
%--------------------------------------------------------------------------

% freq.time=round(freq.time*100)/100; %get rid of tiny differences in time axis
% frind = find((freq.freq>=CUTLO) & (freq.freq<=CUTHI));
% faxis = freq.freq(frind);
tind = find((timelock.time>=TIMLO) & (timelock.time<=TIMHI));
taxis = timelock.time(tind);
chlabel = timelock.label;

basetind = find(taxis >= BASELO & taxis <= BASEHI);

            % 1     2    3      4        5          6       7       8       9       10   
            % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) 
subjrespavg = nan( length(chlabel),length(taxis),2,2,3,3,3,1);
subjrespvar = nan( length(chlabel),length(taxis),2,2,3,3,3,1);

for ises = 1:4
    fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
    fprintf('Concatenating runs Subject %s Session %d: %s  . . .\n', SUBJ, ises, sesdirs{ises})

    if ises==1, ipharm = 1;   imotor=1;
    elseif ises==2, ipharm = 1;   imotor=2;
    elseif ises==3, ipharm = 2;   imotor=1;
    elseif ises==4, ipharm = 2;   imotor=2;
    end
    
    sesdir = fullfile(subjdir, sesdirs{ises});
    
    runlist = dir([sesdir '/*data.mat']);
    if isempty(runlist);     continue;            end
    
    cd(sesdir)
    
    trialdat = [];
    trialinfo = [];
    %concat runs
    fprintf('Loading preprocessed data and concatenating . . .\n')
    for irun = 1:length(runlist)
        fprintf('Loading run %d . . .\n', irun)
        load(runlist(irun).name); % subj pharma diff pupil SDT

        data = ft_megplanar(cfg_megplanar, data);

        ntrials = length(data.trial);
        timemat = repmat(TIMLO:1/resamplefs:TIMHI,ntrials,1); % specify a time axis on which you want the data to be resampled.
        cfg_downsample.time = mat2cell(timemat, ones(ntrials,1), size(timemat,2))';
        data = ft_resampledataNK(cfg_downsample,data); %edit: put nan where no data at cfg_downsample.time
        
        timelock= ft_timelockanalysis(cfg_timelock,data);

        trialdat = [trialdat; timelock.trial]; 
        trialinfo = [trialinfo; timelock.trialinfo];
    end
    ntrl = size(trialdat,1);
    
    %% Split up by condition and average over trials
    fprintf('Split up by condition and average over trials\n')
    % freq.trialinfo   % 1 difficulty condition% 2 signal location% 3 button% 4 correct % 5 RT% 6 empty% 7 trial index% 8 N button presses in trial (omission = 0)
    
    % 1     2    3      4        5          6       7       8       9       10
    % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)

    %     for idiff=1:2 % easy, hard
%         for istim = 1:2 % 1 = left, 2 = right
%             for ichoice = 1:2 % left or right
%                 fprintf('.')
%                 trialind = find( trialinfo(:,2) == istim &  trialinfo(:,1) == idiff & trialinfo(:,3) == ichoice );
%                 timelock.trial = var(trialdat(trialind,:,:));
%                 timelock_cmb = ft_combineplanar([], timelock);
%                 
%                 subjrespavg(:,:, ipharm, imotor, idiff, istim, ichoice, 1) = timelock_cmb.trial; %in resp 3
% %                 subjrespavg(:,:,:, ipharm, imotor, idiff, istim, ichoice, 1) = nanmean(respdat(trialind,:,:,:)); %in resp 3
%             end
%         end
%     end
    fprintf('\ncollapse over diff, stim and choice . . .\n')
    timelock.trial = nanvar(trialdat);
    timelock_cmb = ft_combineplanar([], timelock);

    bl=mean(timelock_cmb.avg(:,basetind),2);
    baselinedat = repmat(bl, size(timelock_cmb.time));
    
    subjrespvar(:,:, ipharm, imotor, 3, 3, 3, 1) = timelock_cmb.avg - baselinedat; %in resp 3
    % compute ERF:
    timelock.trial = nanmean(trialdat);
    timelock_cmb = ft_combineplanar([], timelock);
    subjrespavg(:,:, ipharm, imotor, 3, 3, 3, 1) = timelock_cmb.avg; %in resp 3

%     
%     fprintf('\ncollapse over diff and stim . . .\n')
%     for ichoice = 1:2 % left or right
%         fprintf('.')
%         trialind = find(  trialinfo(:,3) == ichoice );
%         subjrespavg(:,:,:, ipharm, imotor, 3, 3, ichoice, 1) = nanmean(respdat(trialind,:,:,:));
%     end
%     fprintf('\ncollapse over choice and stim . . .\n')
%     for idiff = 1:2 % left or right
%         fprintf('.')
%         trialind = find(  trialinfo(:,1) == idiff );
%         subjrespavg(:,:,:, ipharm, imotor, idiff, 3, 3, 1) = nanmean(respdat(trialind,:,:,:));
%     end
%     fprintf('\ncollapse over diff and choice . . .\n')
%     for istim = 1:2 % left or right
%         fprintf('.')
%         trialind = find(  trialinfo(:,2) == istim );
%         subjrespavg(:,:,:, ipharm, imotor, 3, istim, 3, 1) = nanmean(respdat(trialind,:,:,:));
%     end
%     
%     fprintf('\ncollapse over diff . . .\n')
%     for istim = 1:2 % left or right
%         for ichoice = 1:2 % left or right
%             fprintf('.')
%             trialind = find(  trialinfo(:,2) == istim &  trialinfo(:,3) == ichoice);
%             subjrespavg(:,:,:, ipharm, imotor, 3, istim, ichoice, 1) = nanmean(respdat(trialind,:,:,:));
%         end
%     end
%     fprintf('\ncollapse over stim . . .\n')
%     for idiff = 1:2 % left or right
%         for ichoice = 1:2 % left or right
%             fprintf('.')
%             trialind = find(  trialinfo(:,1) == idiff &  trialinfo(:,3) == ichoice);
%             subjrespavg(:,:,:, ipharm, imotor, idiff, 3, ichoice, 1) = nanmean(respdat(trialind,:,:,:));
%         end
%     end
%     fprintf('\ncollapse over choice . . .\n')
%     for istim = 1:2 % left or right
%         for idiff = 1:2 % left or right
%             fprintf('.')
%             trialind = find(  trialinfo(:,2) == istim &  trialinfo(:,1) == idiff);
%             subjrespavg(:,:,:, ipharm, imotor, idiff, istim, 3, 1) = nanmean(respdat(trialind,:,:,:));
%         end
%     end
%     
end % ises
% clear respdat
% 
% save respavg per subj
filesave = fullfile(respavgout, [SUBJ '_' baseline '.mat']);
fprintf('Saving %s . . . \n', filesave)
% save(filesave, 'subjrespavg')
save(filesave) % save all vars
% 


% %% 
% close all
% load('sensorselection.mat')
% % SOINsel = { 'occ-l'    'occ-r'  'occ-spec' 'frontal'}; %  'motor-l' 'motor-r' %         SOIN = {'occ-l', 'occ-r', 'occ-spec', 'frontal'};
% % SOINind = find(strcmp(SOINsel{1}, {chans.group}));
% sensind = chans(1).sens; %occ
% 
% resamplefs = 100;
% TIMLO = -0.5;
% TIMHI = 1.5;
% taxis = TIMLO:1/resamplefs:TIMHI;
% 
% dat = squeeze(mean(subjrespvar(sensind,:,:,:,3,3,3)));
% dat(:,:,3) = squeeze(mean(dat,3));
% figure; hold on
% COL = ['g', 'k']
% for ipharm=1:2
%     for imotor=3 %1:2 %2%:2
%         plot(taxis, dat(:,ipharm, imotor), COL(ipharm))
%     end
% end
%         legend({'drug', 'placebo'})
%         
