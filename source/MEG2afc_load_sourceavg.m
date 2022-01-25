% function [ output_args ] = MEG2afc_load_respavg( trigger )
%MEG2afc_load_respavg Load subjrespavg for each subject, create condition
%labels, contrast to prepare for plotting etc. occ and motor poolings

% subjrespavg made by MEG2afc_concat_runs
% MEG hh conditions:
% drug placebo
% easy hard,
% pupil hi lo

close all

trigger = 'stim';
disp(trigger)

SUBJ  = { 
    'NK1'     'NK3'   'NK4'   'NK5'   'NK7'   'NK8'   'NK9'  'NK11' ...    
    'NK12'   'NK13'  'NK15'...
         'NK16'  'NK17'           'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK14' no mri 'NK6' out,  'NK10' don't exist, 'NK2' incomplete  'NK12' 'NK15' have bad ddm fits

% SUBJ  = { 
%        'NK20'
%     }; %  'NK6' out,  'NK10' don't exist, 'NK2' incomplete

nsub=length(SUBJ);

baseline = 'trial'; analysistype = 'full';

% set paths
% restoredefaultpath
if ismac
    matlabpath = '/Users/kloosterman/gridmaster2012/MATLAB'; % '/mnt/homes/home022/nkloost1';
    basepath = '/Users/kloosterman/gridmaster2012'; % '/mnt/homes/home022/nkloost1';
else
    matlabpath = '/home/mpib/kloosterman/MATLAB'; % '/mnt/homes/home022/nkloost1';
    basepath = '/home/mpib/kloosterman'; % '/mnt/homes/home022/nkloost1';
end
addpath(genpath(fullfile(matlabpath, 'MEG_HH_analysis')))
addpath(fullfile(matlabpath, 'fieldtrip-20141231'))
ft_defaults

PREIN = fullfile(basepath, 'projectdata/MEG2afc/source', trigger);
PREOUT = fullfile(basepath, 'plots'); mkdir(PREOUT)

% example subjrespavg, put in dat
subjrespavgin = sprintf('%s/respavg/%s_%s.mat', PREIN, SUBJ{1}, baseline);
fprintf('loading subrespavg %s\n', subjrespavgin)
dat = load(subjrespavgin);

% set up arrays for analysis
%--------------------------------------------------------------------------
% FREQLO = 5; FREQHI = 149;
% faxis = dat.faxis;
% frind = find((faxis>=FREQLO) & (faxis<=FREQHI));
% if strcmp(trigger, 'stim') % this timewindow is cut out from the subjrespavg timewin
% %     TIMLO          = -0.2;    TIMHI          =  0.4; % TFR's
%     TIMLO          = -0.2;    TIMHI          =  1; %ramping
% else
% %     TIMLO          = -0.4 ;   TIMHI          =  0.2;
%     TIMLO          = -1 ;   TIMHI          =  0.2;
% end
% tind = find((dat.taxis>=TIMLO) & (dat.taxis<=TIMHI));
% taxis = dat.taxis(tind);
% chlabel = dat.chlabel;
% 1     2    3      4        5          6       7       8       9         10          
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
% respavg = nan(length(SUBJ),length(chlabel),length(frind),length(taxis),2,2,3,3,3,3, 'single');
% respavg = nan(length(SUBJ),6,length(frind),length(taxis),2,2,3,3,3,3,4,3, 'single'); % only poolings
% respavg = nan(length(SUBJ),6,length(frind),length(taxis),2,2,3,3,3,3,4, 'single'); % only poolings,  % collapse over correct
% respavg = nan(length(SUBJ),6,length(faxis),length(taxis),2,2,3,3,3,3, 'single'); % only poolings, collapse over pupil and RT
sourceavg = nan([length(SUBJ), size(dat.subj_sourceavg)]); %  SUBJ x y z idrug ireg ichoice
% respvar = nan(length(SUBJ),2,2,3,3,3,3); % only poolings, collapse over RT and correct
%         latrvar = nan(length(SUBJ), 2, 2,2,length(frind)+3,length(taxis),2,2,3); % var/mean isoi iwrt freq time drug regime
NKsensorselection

% Condition labels
%--------------------------------------------------------------------------
pharm_conds = {'atomox' 'placebo' 'drugscomb' 'atomox - placebo'};
motor_conds = {'pressipsi' 'presscontra' 'regimecomb' 'pressipsi - presscontra'};
% diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};
% pupilconds = {'pupilhigh' 'pupillow' 'pupilcomb' 'pupilhi-lo'};
% stim_conds = {'stimleft' 'stimright' 'allstim' 'stimleft-right'};
choice_conds = {'choiceleft' 'choiceright' 'allchoice' 'choiceleft-right'};
% rt_conds = {'slow' 'medium' 'fast' 'allrts'};
% correct_conds = {'correct' 'error' 'corr+err'};
sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj

for isub = 1:length(SUBJ)
    subjrespavgin = sprintf('%s/respavg/%s_%s.mat', PREIN, SUBJ{isub}, baseline);
    fprintf('loading subrespavg %s\n', subjrespavgin)
    dat = load(subjrespavgin);
    
    % MEG2afc_concat_runs structure:
%     subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp, ipup, irt, icor) 
    
%     respavg(isub,1,:,:, :,:,:,:,:,:) = mean(dat.subjrespavg(occind,  :,tind, :,:,:,:,:,:),1);
%     respavg(isub,2,:,:, :,:,:,:,:,:) = mean(dat.subjrespavg(motorind,:,tind, :,:,:,:,:,:),1);
%     respavg(isub,:,:,:, :,:,:, :,:,: ,:) = dat.subjrespavg(:,:,:, :,:,:, :,:,: ,:,3); % collapse over correct
%     respavg(isub,:,:,:, :,:,:, :,:,: ,:) = dat.subjrespavg(:,frind,tind, :,:,:, :,:,: ,4,3); % collapse over RT and correct
      sourceavg(isub,:,:,:, :,:,:) = dat.subj_sourceavg;
%     latrvar(isub,:,:, :,:,:,:,:,:) = subjlatrvar;
%     respvar(isub, :,:,:, :,:,:) = dat.subjrespvar(:,:,:,:,:,:,4,3);
    %         length(find(isnan(subjrespavg)))
%     clear subjrespavg %subjlatrvar
end

sourceavg(:,:,:,:, 3,:,:) = mean(sourceavg(:,:,:,:, 1:2,:,:),5);
sourceavg(:,:,:,:, 4,:,:) = sourceavg(:,:,:,:, 1,:,:) - sourceavg(:,:,:,:, 2,:,:);
sourceavg(:,:,:,:, :,3,:) = mean(sourceavg(:,:,:,:, :,1:2,:),6);
sourceavg(:,:,:,:, :,:,4) = sourceavg(:,:,:,:, :,:,1) - sourceavg(:,:,:,:, :,:,2);


source_avg_dimord = 'nsub_x_y_z_regime_drug_choice'


%% plotting average
SAV=1;
source = dat.source_diff_int{1};
%     sourceall_avg = source_diff_int{3};
%     sourceall_avg.avg.pow = squeeze(mean(pow));

idrug = 2; ireg = 3;

viewpoints = [-45 0; 45 0];
close all;
iview = 1;
    load('colormap170613.mat'); % colormap(cmap);  %cmap = get(gcf, 'Colormap')

cfg               = [];
cfg.method        = 'ortho'; % 'slice',  'ortho', 'surface', 'glassbrain'
cfg.funparameter  = 'pow';  %avg.pow
cfg.surfinflated = 'surface_inflated_both_caret.mat';
%     cfg.surfinflated = 'surface_inflated_right_caret.mat';
%     cfg.surfinflated = 'surface_pial_both.mat';
%     cfg.surfinflated = 'surface_white_both.mat';
cfg.maskparameter = cfg.funparameter;
%     cfg.maskparameter = stat.mask;
%     cfg.funcolorlim   = [-0.12 0.12];
%     cfg.opacitylim    = [0 0.12];
%     cfg.opacitymap    = stat.mask(:); %'rampup';
cfg.opacitymap    = 'rampup';
cfg.funcolormap = cmap;
cfg.camlight ='yes';

for ichoice = 4
    load('colormap170613.mat'); % colormap(cmap);  %cmap = get(gcf, 'Colormap')
    
    source.avg.pow = squeeze(mean(sourceavg(:,:,:,:,idrug,ireg,ichoice)));
    if strcmp(cfg.method, 'surface')
        for iview=1%:2
            
            
            ft_sourceplot(cfg,source);
            view([0 0])
            camlight('headlight')
            view(viewpoints(iview, :))
            title({sprintf('%s %s %slocked', pharm_conds{idrug}, choice_conds{ichoice}, trigger) ;
                sprintf('%g-%g s, %g Hz, smoothing %d', 0.2, 0.4, ft_findcfg(source.cfg, 'foi'), ft_findcfg(source.cfg, 'tapsmofrq') )}); % 'FontWeight','bold'
            if SAV
                outpath = fullfile(PREOUT, 'source'); warning off; mkdir(outpath); warning on
                outfile = fullfile(outpath, sprintf( 'source_%slocked_%s_%s_view%d_%s', trigger, pharm_conds{idrug}, choice_conds{ichoice}, iview, cfg.method )); % motor_conds{imotor}, pupilconds{ipup}
                disp(outfile)
                export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',                
            end
        end
    else
        ft_sourceplot(cfg,source);
        title({sprintf('%s %s %slocked', pharm_conds{idrug}, choice_conds{ichoice}, trigger) ;
            sprintf('%g-%g s, %g Hz, smoothing %d', ft_findcfg(source.cfg, 'exp_interval'), ft_findcfg(source.cfg, 'foi'), ft_findcfg(source.cfg, 'tapsmofrq') )}); % 'FontWeight','bold'
    end
    if SAV && ~strcmp(cfg.method, 'surface')
        outpath = fullfile(PREOUT, 'source'); warning off; mkdir(outpath); warning on
        outfile = fullfile(outpath, sprintf( 'source_%slocked_%s_%s_view%d_%s', trigger, pharm_conds{idrug}, choice_conds{ichoice}, iview, cfg.method)); % motor_conds{imotor}, pupilconds{ipup}
        disp(outfile)
        export_fig(outfile, '-png', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
    end
    
end

%% run statistics over subjects %
source_drug = {}; source_plac = {};
for isub=1:nsub
    source_drug{isub} = source;
    source_drug{isub}.avg.pow = sourceavg(isub,:,:,:, 3,3,3);
    
    source_plac{isub} = source;
    source_plac{isub}.avg.pow = zeros(size(sourceavg(isub,:,:,:, 3,3,3)));
end
    
 %%   
cfg=[];
cfg.dim         = source_drug{1}.dim;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.parameter   = 'avg.pow';
cfg.correctm    = 'cluster';
cfg.numrandomization = 100;
cfg.alpha       = 0.05; % note that this only implies single-sided testing
cfg.tail        = 0;

nsubj=numel(source_drug);
cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj) ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfg, source_drug{:}, source_plac{:});


