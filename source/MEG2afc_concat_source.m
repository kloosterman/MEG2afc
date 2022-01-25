function [ subj_sourceavg ] = MEG2afc_concat_source( subjdir, trigger )
%Concatenate source for each subject and save subjrespavg

if ismac
    matlabpath = '/Users/kloosterman/gridmaster2012/MATLAB'; % '/mnt/homes/home022/nkloost1';
else
    matlabpath = '/home/mpib/kloosterman/MATLAB'; % '/mnt/homes/home022/nkloost1';
end
addpath(genpath(fullfile(matlabpath, 'MEG_HH_analysis')))
addpath(fullfile(matlabpath, 'fieldtrip-20141231'))
ft_defaults

[~,SUBJ] = fileparts(subjdir);

sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
sesdir_codes = ['B' 'D' 'A' 'C'];

analysistype = 'full';
baseline='trial';
AVG = 'totalpow';

if ismac
    PREIN = fullfile('/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/source', trigger);
else
    PREIN = fullfile('/home/mpib/kloosterman/projectdata/MEG2afc/source', trigger);
end
respavgout = fullfile(PREIN, 'respavg');
mkdir(respavgout)

fprintf('%s; %s; %s %s\n', trigger, analysistype, baseline, AVG)

% source analysis settings:
cfg=[];
cfg.trial_duration = [-0.4 1];
% cortical stimulus response, gamma
cfg.bsl_interval = [-0.4 -0.2];
cfg.exp_interval = [0.2 0.4];
cfg.tapsmofrq = 20;
cfg.foi = 50;


for ises = 1:4
    fprintf('\n\nSubject directory: %s  . . .\n', subjdir)
    fprintf('Concatenating source Subject %s Session %d: %s: %s  . . .\n', SUBJ, ises, sesdir_codes(ises), sesdirs{ises})

    if ises==1, ipharm = 1;   imotor=1;
    elseif ises==2, ipharm = 1;   imotor=2;
    elseif ises==3, ipharm = 2;   imotor=1;
    elseif ises==4, ipharm = 2;   imotor=2;
    end
    
    sesdir = fullfile(subjdir, sesdirs{ises});
    runlist = dir([sesdir '/*_source_*']);
    if isempty(runlist);     continue;            end    
    cd(sesdir)
    
%     powdat = [];    trialinfo = [];
    fprintf('Loading source and averaging . . .\n')
    
    sourceall = {};
    for irun = 1:length(runlist)
        fprintf('Loading run %d: %s . . .\n', irun, runlist(irun).name)
        load(runlist(irun).name);
        for ichoice = 1:3
            powdat(irun,:,:,:,ichoice) = source_diff_int{ichoice}.avg.pow;
        end
        %         sourceall{irun} = source_diff_int{3} % subj pharma diff pupil SDT
    end
    
    if size(powdat,1) > 1
        subj_sourceavg(:,:,:,ipharm,imotor,:) = squeeze(mean(powdat)); %average over runs
    else % only 1 run
        subj_sourceavg(:,:,:,ipharm,imotor,:) = squeeze(powdat); %average over runs
    end        
end
outputfile = sprintf('%s_source_%dto%d_vs_%dto%d_foi%d_smo%d.mat', SUBJ{1}, cfg.exp_interval*1000, cfg.bsl_interval*1000, cfg.foi, cfg.tapsmofrq);
filesave = fullfile(respavgout, outputfile);
fprintf('Saving %s . . . \n', filesave)
save(filesave) % save all vars

% %     
% %     %     cfg=[];
% %     %     cfg.keepindividual = 'no';
% %     % %     cfg.parameter = 'pow';
% %     %     sourceall_avg = ft_sourcegrandaverage(cfg, sourceall{:});
% %     %     sourceall_avg.pow = reshape(sourceall_avg.pow, sourceall_avg.dim );
% %     %     sourceall_avg.anatomy = source_diff_int{3}.anatomy;
    
%     sourceall_avg = source_diff_int{3};
%     sourceall_avg.avg.pow = squeeze(mean(pow));
%     
%     load('colormap170613.mat'); % colormap(cmap);  %cmap = get(gcf, 'Colormap')
%     
%     close all;
%     cfg               = [];
%     cfg.method        = 'slice'; % 'slice',  'ortho', 'surface',
%     cfg.funparameter  = 'pow';  %avg.pow
%     cfg.maskparameter = cfg.funparameter;
%     cfg.funcolorlim   = [0.0 0.5];
%     cfg.opacitylim    = [0.0 0.25];
%     cfg.opacitymap    = 'rampup';
%     %             cfg.funcolormap = cmap(129:256,:);
%     
%     cfg.camlight ='yes';
%     ft_sourceplot(cfg,sourceall_avg);
%             view([-45 0])

    



