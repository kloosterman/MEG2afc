function MEG2afc_dva_analysis( subjdir, trigger )
%Run Schurger dva analysis

cd('~/qsub');

[~,SUBJ] = fileparts(subjdir);

sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
% sesdir_codes = ['B' 'D' 'A' 'C'];
% 
% fprintf('%s; %s; %s %s\n', trigger)
% 
% exampledatapath = fullfile(subjdir, sesdirs{1});
% w = what(exampledatapath);
% load(fullfile(exampledatapath, w.mat{1}))
% 
% % move trigger if resplocked analysis
% if strcmp(trigger, 'resp')
%     cfg_resp=[];
%     cfg_resp.offset = -data.trialinfo(:,4);
%     cfg_resp.trials = find(cfg_resp.offset < 1);
%     data = ft_redefinetrial(cfg_resp, data);
% end

% specify a time axis on%   which you want the data to be resampled.
fsample_desired = 250;
if strcmp(trigger, 'stim')
    taxis = -1.0:1/fsample_desired:3;
else
    taxis = -1.5:1/fsample_desired:2;
end

% time = data.time;
% [time{:}] = deal(taxis);

% cfg=[];
% cfg.resample = 'yes';
% cfg.fsample = 500;
% cfg.time = time; % instead of % cfg.resamplefs = 250;
% data = ft_resampledata(cfg,data);
% 
% cfg=[];
% cfg.channel = 'MEG';
% cfg.vartrllength       =  2;
% cfg.keeptrials         = 'yes';
% timelock = ft_timelockanalysis(cfg, data);
% 
% tind = 1:length(timelock.time);
% taxis = timelock.time(tind);
% chlabel = timelock.label;

% 1     2    3      4        5          6       7       8       9       10           11     12
subjrespavg = nan( 4, 3, length(taxis),2,2,3,3,3, 'single' ); % across var, across norm, within var, within norm

for ises = 1:4
    fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
    fprintf('Concatenating runs Subject %s Session %d: . . .\n', SUBJ, ises)
    
    if ises==1, ipharm = 1;   imotor=1;
    elseif ises==2, ipharm = 1;   imotor=2;
    elseif ises==3, ipharm = 2;   imotor=1;
    elseif ises==4, ipharm = 2;   imotor=2;
    end
    
    sesdir = fullfile(subjdir, sesdirs{ises});
    runlist = dir(fullfile(sesdir, '*_data.mat'));
    if isempty(runlist);     continue;            end
    cd(sesdir)
    
    fprintf('Loading data and concatenating . . .\n')
    
    data = cell(size(runlist,1),1);
    for irun = 1:length(runlist)
        temp = load(runlist(irun).name);
        data{irun} = temp.data;
        clear temp
    end
    data = ft_appenddata([], data{:});

    if strcmp(trigger, 'resp')
        cfg_resp=[];
        cfg_resp.offset = -data.trialinfo(:,5);
        cfg_resp.trials = find(cfg_resp.offset < 1);
        data = ft_redefinetrial(cfg_resp, data);
    end
    
    % specify a time axis on%   which you want the data to be resampled.
%     fsample_desired = 250;
%     if strcmp(trigger, 'stim')
%         taxis = -1.0:1/fsample_desired:3;
%     else
%         taxis = -1.5:1/fsample_desired:2;
%     end
    time = data.time;
    [time{:}] = deal(taxis); % defined above
    cfg=[];
    cfg.resample = 'yes';
    cfg.fsample = 500;
    cfg.time = time; % instead of % cfg.resamplefs = 250;
    data = ft_resampledata(cfg,data);
    
    % select trials that go at least until 2 s
    % first, cut off all the trials at 2 s
    cfg=[];
    cfg.latency = [-1 2];
    data = ft_selectdata(cfg, data);
    % then, take only trials containing data at all tp
    cfg=[];
    cfg.channel = 'MEG';
    cfg.vartrllength       =  1; % 1: only trials having data in the avg
    cfg.keeptrials         = 'yes';
    timelock = ft_timelockanalysis(cfg, data);
    taxis = timelock.time;
    
    %%
    % Split up by condition and average over trials
    % freq.trialinfo   % 1 difficulty condition% 2 signal location% 3 button% 4 correct % 5 RT% 6 empty% 7 trial index% 8 N button presses in trial (omission = 0)
    
    trialinfo = timelock.trialinfo;
    NKsensorselection
    dvaleg = {'across dva', 'across norm', 'within dva', 'within norm'};
    % 1     2    3      4        5          6       7       8       9       10           11        12
    % subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) RT bin(4) correct (3)
    fprintf('Split up by condition and average over trials\n')
    ctr=0;
    diffs = 1:3;    stims = 1:3;    resps = 1:3;    isois = 1:3;
    ft_progress('init', 'etf',     'Please wait...');
    nconds=length(diffs) * length(stims) * length(resps) * length(isois);
    for isoi=isois % occipital, motor, all
        for idiff = diffs % easy, hard
            for istim = stims % 1 = left, 2 = right
                for iresp = resps % left or right
                    ctr=ctr+1;
                    diff_ind = trialinfo(:,1) == idiff;
                    if ~any(diff_ind), diff_ind(:)=true; end
                    stim_ind = trialinfo(:,2) == istim;
                    if ~any(stim_ind), stim_ind(:)=true; end
                    resp_ind = trialinfo(:,3) == iresp;
                    if ~any(resp_ind), resp_ind(:)=true; end
                    
                    trial_inds = diff_ind & stim_ind & resp_ind;
                    
                    ft_progress(ctr/nconds, 'Processing diff %d stim %d resp %d soi %d: %d trials', idiff, istim, iresp, isoi, ...
                        length(find(trial_inds)));
                    
                    if length(find(trial_inds)) < 1
                        warning('No trials remain for this condition')
                        continue
                    end
                    
                    %select channels
                    chans = sens.ind{isoi};
                    
                    %dva analysis per condition
                    % across-trial variability
                    dvadat = [];
                    for it = 1:length(taxis)
                        vecSet = squeeze(timelock.trial(trial_inds,chans,it))';
                        [dvadat(it), mn(it)] = dva(vecSet);
                    end
                    %                             figure; hold on
                    %                             plot(taxis, dvadat(:,:))
                    subjrespavg(1,isoi, :, ipharm, imotor, idiff, istim, iresp) = dvadat;
                    subjrespavg(2,isoi, :, ipharm, imotor, idiff, istim, iresp) = mn;
                    
                    % within-trial variability
                    winSize = 25; % 100 ms @ 250 Hz
                    ntrials = size(timelock.trial,1);
                    dvadat = nan(ntrials, length(timelock.time));
                    L2norm = nan(ntrials, length(timelock.time));
                    for itrial= find(trial_inds)'
                        vecSet = squeeze(timelock.trial(itrial,chans,:));
                        [swv, swn] = fast_sw_dva(vecSet,winSize);
                        dvadat(itrial,:) = swv;
                        L2norm(itrial,:) = swn;
                    end
                    %                             figure; hold on
                    %                             plot(timelock.time, nanmean(dvadat(:,:)))
                    subjrespavg(3,isoi, :, ipharm, imotor, idiff, istim, iresp) = nanmean(dvadat);
                    subjrespavg(4,isoi, :, ipharm, imotor, idiff, istim, iresp) = nanmean(L2norm);
                    
                end
            end
        end
    end
    ft_progress('close')
    clear data timelock L2norm dvadat vecset
end % ises
% 
% figure; hold on; plot(taxis, subjrespavg(3,:,1,1,3,3,3))
% plot(taxis, subjrespavg(3,:,2,1,3,3,3))
% figure; hold on; plot(taxis, subjrespavg(1,:,1,1,3,3,3))
% plot(taxis, subjrespavg(1,:,2,1,3,3,3))

respavgout = fullfile(fileparts(subjdir), 'dva', trigger);
mkdir(respavgout)

% save respavg per subj
filesave = fullfile(respavgout, [SUBJ '_dva.mat']);
fprintf('Saving %s . . . \n', filesave)
% save(filesave, 'subjrespavg')
save(filesave) % save all vars



