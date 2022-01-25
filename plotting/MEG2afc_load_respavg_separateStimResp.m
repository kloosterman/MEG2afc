% function [ output_args ] = MEG2afc_load_respavg( trigger )
%MEG2afc_load_respavg Load subjrespavg for each subject, create condition
%labels, contrast to prepare for plotting etc. occ and motor poolings

% subjrespavg made by MEG2afc_concat_runs
% MEG hh conditions:
% drug placebo
% easy hard,
% pupil hi lo

close all
disp(trigger)

% ALL Subjects:
SUBJ  = { 
    ... 'NK1' 
    'NK2' 'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
         'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' have bad ddm fits

% % Subjects who have good ddm fits:
% SUBJ  = { 
%     'NK1'  'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK13'   'NK14' ...  
%          'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
%     }; % 'NK2' incomplete 'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' have bad ddm fits

nsub=length(SUBJ);

baseline = 'trial'; 
% analysistype = 'full';
analysistype = 'low';
baselinetype = 'respavg_baselinepersession';

%%%set paths
if ismac
    MATLABPATH = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/';
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/'; % on the cluster
else
    MATLABPATH = '/home/mpib/kloosterman/MATLAB/';
    basepath = '/home/mpib/kloosterman'; % on the cluster
end
addpath(genpath(fullfile(MATLABPATH, 'MEG_HH_analysis')))
addpath(fullfile(MATLABPATH, 'toolbox', 'qsub_tardis'))
addpath(fullfile(MATLABPATH, 'toolbox/fieldtrip-20150803')) %inc JJ edit ft_artifact_zvalue
ft_defaults

PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq', analysistype, trigger);
PREOUT = fullfile(basepath, 'plots'); mkdir(PREOUT)

% example subjrespavg, put in dat
subjrespavgin = sprintf('%s/%s/%s_%s.mat', PREIN, baselinetype, SUBJ{1}, baseline);
fprintf('loading subrespavg %s\n', subjrespavgin)
dat = load(subjrespavgin);

% set up arrays for analysis
%----------------- ---------------------------------------------------------
FREQLO = 2; FREQHI = 149;
faxis = dat.faxis;
frind = find((faxis>=FREQLO) & (faxis<=FREQHI));
TIMLO = dat.TIMLO;
TIMHI = dat.TIMHI;
% if strcmp(trigger, 'stim') % this timewindow is cut out from the subjrespavg timewin
% %     TIMLO          = -0.2;    TIMHI          =  0.4; % TFR's
%     TIMLO          = -0.2;    TIMHI          =  1; %ramping
% else
%     TIMLO          = -0.4 ;   TIMHI          =  0.2;
% %     TIMLO          = -1 ;   TIMHI          =  1;
% end
tind = find((dat.taxis>=TIMLO) & (dat.taxis<=TIMHI));
taxis = dat.taxis(tind);
chlabel = dat.chlabel;
% 1     2    3      4        5          6       7       8       9         10          
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
% respavg = nan(length(SUBJ),6,length(faxis),length(taxis),2,2,3,3,3,3, 'single'); % only poolings, collapse over pupil and RT

% subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp) 
% respavg = nan(length(SUBJ),length(chlabel),length(faxis),length(taxis),2,2,3,3,3, 'single'); % only poolings, collapse over pupil and RT
respavg = nan(length(SUBJ),length(chlabel),length(faxis),length(taxis),4,3,3,3,3, 'single'); % only poolings, collapse over pupil and RT

NKsensorselection

% Condition labels
%--------------------------------------------------------------------------
pharm_conds = {'atomox' 'placebo' 'drugscomb' 'atomox - placebo'};
motor_conds = {'pressipsi' 'presscontra' 'regimecomb' 'pressipsi - presscontra'};
diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};
pupilconds = {'pupilhigh' 'pupillow' 'pupilcomb' 'pupilhi-lo'};
stim_conds = {'stimleft' 'stimright' 'allstim' 'stimleft-right'};
choice_conds = {'choiceleft' 'choiceright' 'allchoice' 'choiceleft-right'};
rt_conds = {'slow' 'medium' 'fast' 'allrts'};
correct_conds = {'correct' 'error' 'corr+err'};
sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj

for isub = 1:length(SUBJ)
    subjrespavgin = sprintf('%s/%s/%s_%s.mat', PREIN, baselinetype, SUBJ{isub}, baseline);

    fprintf('loading subrespavg %s\n', subjrespavgin)
    dat = load(subjrespavgin);
    
    % MEG2afc_concat_runs structure:
%     subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp, ipup, irt, icor) 
    
%     respavg(isub,1,:,:, :,:,:,:,:,:) = mean(dat.subjrespavg(occind,  :,tind, :,:,:,:,:,:),1);
%     respavg(isub,2,:,:, :,:,:,:,:,:) = mean(dat.subjrespavg(motorind,:,tind, :,:,:,:,:,:),1);
%     respavg(isub,:,:,:, :,:,:, :,:,: ,:) = dat.subjrespavg(:,:,:, :,:,:, :,:,: ,:,3); % collapse over correct
%     respavg(isub,:,:,:, :,:,:, :,:,: ,:) = dat.subjrespavg(:,frind,tind, :,:,:, :,:,: ,4,3); % collapse over RT and correct
      respavg(isub,:,:,:, 1:2,1:2,1:3,1:3,1:3) = dat.subjrespavg(:,frind,tind,:, :,:, :,:);
%     latrvar(isub,:,:, :,:,:,:,:,:) = subjlatrvar;
%     respvar(isub, :,:,:, :,:,:) = dat.subjrespvar(:,:,:,:,:,:,4,3);
    %         length(find(isnan(subjrespavg)))
%     clear dat
end

% [nsub, nchan,nfrq,ntim,npharm,nmotor,ndiff,nstim, nchoice,npup]=size(respavg)
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nchoice_npup_nRT_ncorrect'
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nresp_npup'
respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nresp'
size(respavg)

XLIM = [TIMLO TIMHI];
YLIM = [FREQLO FREQHI];
faxis = dat.faxis

%% LATERALIZATION 
NKsensorselection

pooling=[];
% average over sensor groups FIRST.
pooling(:,1,:,:, :,:,:,:,:) =  nanmean(respavg(:,leftoccind,:,:, :,:,:,:,:), 2);
pooling(:,2,:,:, :,:,:,:,:) =  nanmean(respavg(:,rightoccind,:,:, :,:,:,:,:), 2);

pooling(:,3,:,:, :,:,:,:,:) =  nanmean(respavg(:,leftmotorind,:,:, :,:,:,:,:), 2);
pooling(:,4,:,:, :,:,:,:,:) =  nanmean(respavg(:,rightmotorind,:,:, :,:,:,:,:), 2);

% LATERALIZATION OCC WRT CHOICE
% choice 1 - choice 2
% ses1 (ipsi): stim1 and resp1
% ses1 (ipsi): stim2 and resp2

% ses2 (contra): stim2 and resp1
% ses2 (contra): stim1 and resp2
temp=[];
% ipsi ses:
% if choice is left, do pooling right -left
temp(:,1,:,:, :,1,:,:) = pooling(:,2,:,:, :,1,:,:,1) - pooling(:,1,:,:, :,1,:,:,1);
% if choice is right, do pooling left - right
temp(:,2,:,:, :,1,:,:) = pooling(:,1,:,:, :,1,:,:,2) - pooling(:,2,:,:, :,1,:,:,2);

% same for the contra ses
% if choice is left, do pooling right -left
temp(:,1,:,:, :,2,:,:) = pooling(:,2,:,:, :,2,:,:,2) - pooling(:,1,:,:, :,2,:,:,2);
% if choice is right, do pooling left - right
temp(:,2,:,:, :,2,:,:) = pooling(:,1,:,:, :,2,:,:,1) - pooling(:,2,:,:, :,2,:,:,1);

respavg(:,269,:,:, :,1:2,:,:,3) = mean(temp,2); %average across pooling subtractions
% average across motor mapping down below

clear temp
chlabel{end+1} = 'lateralization_occ_WRTchoice';
sens.leg{end+1} = 'lateralization_occ_WRTchoice';
sens.ind{end+1} = 269;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATERALIZATION OCC WRT button press
temp=[];
% respright: left - right pooling 
temp(:,1,:,:, :,:,:,:) = pooling(:,1,:,:, :,:,:,:,2) - pooling(:,2,:,:, :,:,:,:,2); % respright: left - right pooling 
% respleft: right - left pooling 
temp(:,2,:,:, :,:,:,:) = pooling(:,2,:,:, :,:,:,:,1) - pooling(:,1,:,:, :,:,:,:,1); % respleft: right - left pooling 
respavg(:,270,:,:, :,:,:,:,3) = nanmean(temp,2); %average across pooling subtractions
chlabel{end+1} = 'lateralization_occ_WRTbutton';
sens.leg{end+1} = 'lateralization_occ_WRTbutton';
sens.ind{end+1} = 270;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATERALIZATION OCC WRT stim location
temp=[];
% stimright: left - right pooling 
temp(:,1,:,:, :,:,:,:) = pooling(:,1,:,:, :,:,:,2,:) - pooling(:,2,:,:, :,:,:,2,:); % stimright: left - right pooling 
% stimleft: right - left pooling 
temp(:,2,:,:, :,:,:,:) = pooling(:,2,:,:, :,:,:,1,:) - pooling(:,1,:,:, :,:,:,1,:); % stimleft: right - left pooling 

respavg(:,271,:,:, :,:,:,:,3) =  nanmean(temp,2); %average across pooling subtractions
clear temp
chlabel{end+1} = 'lateralization_occ_WRTtarget';
sens.leg{end+1} = 'lateralization_occ_WRTtarget';
sens.ind{end+1} = 271;

%% DO THE SAME FOR MOTOR CORTEX

% LATERALIZATION MOTOR WRT CHOICE
% choice 1 - choice 2
% ses1 (ipsi): stim1 and resp1
% ses1 (ipsi): stim2 and resp2

% ses2 (contra): stim2 and resp1
% ses2 (contra): stim1 and resp2
temp=[];
% ipsi ses:
% if choice is left, do pooling right -left
temp(:,1,:,:, :,1,:,:) = pooling(:,4,:,:, :,1,:,:,1) - pooling(:,3,:,:, :,1,:,:,1);
% if choice is right, do pooling left - right
temp(:,2,:,:, :,1,:,:) = pooling(:,3,:,:, :,1,:,:,2) - pooling(:,4,:,:, :,1,:,:,2);

% same for the contra ses
% if choice is left, do pooling right -left
temp(:,1,:,:, :,2,:,:) = pooling(:,4,:,:, :,2,:,:,2) - pooling(:,3,:,:, :,2,:,:,2);
% if choice is right, do pooling left - right
temp(:,2,:,:, :,2,:,:) = pooling(:,3,:,:, :,2,:,:,1) - pooling(:,4,:,:, :,2,:,:,1);

respavg(:,272,:,:, :,1:2,:,:,3) = nanmean(temp,2); %average across pooling subtractions
% average across motor mapping down below

clear temp
chlabel{end+1} = 'lateralization_motor_WRTchoice';
sens.leg{end+1} = 'lateralization_motor_WRTchoice';
sens.ind{end+1} = 272;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATERALIZATION MOTOR WRT button press
temp=[];
% respright: left - right pooling 
temp(:,1,:,:, :,:,:,:) = pooling(:,3,:,:, :,:,:,:,2) - pooling(:,4,:,:, :,:,:,:,2); % respright: left - right pooling 
% respleft: right - left pooling 
temp(:,2,:,:, :,:,:,:) = pooling(:,4,:,:, :,:,:,:,1) - pooling(:,3,:,:, :,:,:,:,1); % respleft: right - left pooling 
respavg(:,273,:,:, :,:,:,:,3) = nanmean(temp,2); %average across pooling subtractions
chlabel{end+1} = 'lateralization_motor_WRTbutton';
sens.leg{end+1} = 'lateralization_motor_WRTbutton';
sens.ind{end+1} = 273;

%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATERALIZATION MOTOR WRT stim location
temp=[];
% stimright: left - right pooling 
temp(:,1,:,:, :,:,:,:) = pooling(:,3,:,:, :,:,:,2,:) - pooling(:,4,:,:, :,:,:,2,:); % stimright: left - right pooling 
% stimleft: right - left pooling 
temp(:,2,:,:, :,:,:,:) = pooling(:,4,:,:, :,:,:,1,:) - pooling(:,3,:,:, :,:,:,1,:); % stimleft: right - left pooling 

respavg(:,274,:,:, :,:,:,:,3) =  nanmean(temp,2); %average across pooling subtractions
clear temp
chlabel{end+1} = 'lateralization_motor_WRTtarget';
sens.leg{end+1} = 'lateralization_motor_WRTtarget';
sens.ind{end+1} = 274;





% 1     2    3      4        5          6       7       8       9         10
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
%
respavg(:,:,:,:, :,3,:,:,:) = nanmean(respavg(:,:,:,:, :,1:2,:,:,:),6);
respavg(:,:,:,:, 3,:,:,:,:) = nanmean(respavg(:,:,:,:, 1:2,:,:,:,:),5);
respavg(:,:,:,:, 4,:,:,:,:) = respavg(:,:,:,:, 1,:,:,:,:) - respavg(:,:,:,:, 2,:,:,:,:); % drug contrast
% respavg(:,:,:,:, :,4,:,:,:) = respavg(:,:,:,:, :,1,:,:,:) - respavg(:,:,:,:, :,2,:,:,:); % motor ipsi - contra
respavg(:,:,:,:, :,:,4,:,:) = respavg(:,:,:,:, :,:,1,:,:) - respavg(:,:,:,:, :,:,2,:,:); % easy-hard
% respavg(:,:,:,:, :,:,:,4,:) = respavg(:,:,:,:, :,:,:,1,:) - respavg(:,:,:,:, :,:,:,2,:); % stim left-right
% respavg(:,:,:,:, :,:,:,:,4) = respavg(:,:,:,:, :,:,:,:,1) - respavg(:,:,:,:, :,:,:,:,2); % resp left-right




















% 
% % LATERALIZATION OCC WRT CHOICE
% % choice 1 - choice 2
% % ses1: stim1 and resp1
% % ses1: stim2 and resp2
% 
% % ses2: stim2 and resp1
% % ses2: stim1 and resp2
% respavg_size = size(respavg);
% temp = nan( respavg_size(1:end-2) ); % leave out stim and resp dims
% temp(:,:,:,:, :,1,:) = squeeze(respavg(:,:,:,:, :,1,:,3,1) - respavg(:,:,:,:, :,1,:,3,2)); % ses1 choice left-right, independent of stim
% temp(:,:,:,:, :,2,:) = squeeze(respavg(:,:,:,:, :,2,:,3,2) - respavg(:,:,:,:, :,2,:,3,1)); % ses2 choice left-right, independent of stim
% temp(:,:,:,:, :,3,:) = mean(temp(:,:,:,:, :,1:2,:), 6); % average across sessions
% % add lateralization as extra sensor
% respavg(:,269,:,:, :,:,:,3,3) = mean(temp(:,rightoccind,:,:, :,:,:), 2) - mean(temp(:,leftoccind,:,:, :,:,:), 2); 
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTchoice';
% sens.leg{3} = 'lateralization_occ_WRTchoice';
% sens.ind{3} = 269;
% 
% % LATERALIZATION OCC WRT button press
% % resp 1 -  2
% % respavg_size = size(respavg);
% % temp = nan( respavg_size(1:end-2) ); % leave out stim and resp dims
% temp = squeeze(respavg(:,:,:,:, :,:,:,3,1) - squeeze(respavg(:,:,:,:, :,:,:,3,2))); % button left-right, independent of stim
% % add lateralization as extra sensor
% respavg(:,270,:,:, :,:,:,3,3) = mean(temp(:,rightoccind,:,:, :,:,:), 2) - mean(temp(:,leftoccind,:,:, :,:,:), 2); 
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTbutton';
% sens.leg{4} = 'lateralization_occ_WRTbutton';
% sens.ind{4} = 270;
% 
% % LATERALIZATION OCC WRT stim location
% % resp 1 -  2
% temp = squeeze(respavg(:,:,:,:, :,:,:,1,3)) - squeeze(respavg(:,:,:,:, :,:,:,2,3)); % stim left-right, independent of resp
% % add lateralization as extra sensor
% respavg(:,271,:,:, :,:,:,3,3) = mean(temp(:,rightoccind,:,:, :,:,:), 2) - mean(temp(:,leftoccind,:,:, :,:,:), 2); 
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTtarget';
% sens.leg{5} = 'lateralization_occ_WRTtarget';
% sens.ind{5} = 271;



% %% add lateralization
% % wrt choice, depends on SR mapping
% % pressipsi session
% contra = []; ipsi = [];
% contra(:,1,:,:,: ,:,:,:) = mean(respavg(:,leftoccind,:,:,: ,1,:,:,2),2); %contra1: left sensors, choice right, pressipsi session
% contra(:,2,:,:,: ,:,:,:) = mean(respavg(:,rightoccind,:,:,: ,1,:,:,1),2); % contra2: right sensors, choice left
% ipsi(:,1,:,:,: ,:,:,:) =   mean(respavg(:,leftoccind,:,:,: ,1,:,:,1),2); %ipsi1
% ipsi(:,2,:,:,: ,:,:,:) =   mean(respavg(:,rightoccind,:,:,: ,1,:,:,2),2); % ipsi2
% temp(1,:,:,:,:, :,:,:) =  mean(contra,2) - mean(ipsi,2);
% % presscontra session
% contra(:,1,:,:,: ,:,:,:) = mean(respavg(:,leftoccind,:,:,: ,2,:,:,1),2); %contra1: left sensors, choice right, presscontra session
% contra(:,2,:,:,: ,:,:,:) = mean(respavg(:,rightoccind,:,:,: ,2,:,:,2),2); % contra2: right sensors, choice left
% ipsi(:,1,:,:,: ,:,:,:) =   mean(respavg(:,leftoccind,:,:,: ,2,:,:,2),2); %ipsi1
% ipsi(:,2,:,:,: ,:,:,:) =   mean(respavg(:,rightoccind,:,:,: ,2,:,:,1),2); % ipsi2
% temp(2,:,:,:,:, :,:,:) =  mean(contra,2) - mean(ipsi,2);
%     
% latr = squeeze(mean(temp));
% 
% %% compute Lateralization
% fprintf('\nComputing lateralizations .')
% clear respavglat contra ipsi
% leftind = {leftoccind leftmotorind};
% rightind = {rightoccind rightmotorind};
% 
% % latrleg = {'occ wrt target' 'occ wrt buttonpress' 'motor wrt target' 'motor wrt buttonpress'};
% latrleg = {'occ wrt target' 'occ wrt buttonpress' 'occ wrt choice' 'motor wrt target' 'motor wrt buttonpress' 'motor wrt choice'};
% 
% % Lateralization 10 dims
% sois = [3 4 ; 5 6]; % LR occ then motor
% respavglat= nan(length(SUBJ),6,length(faxis),length(taxis),4,4,4,4,4, 'single'); % only poolings, collapse over RT and correct 
% latctr = 0; 
% for isoi=1:2
%     latctr = latctr+1; contra=[]; ipsi=[];
%     % wrt target loc
%     contra(:,1,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,1),:,:,: ,:,:,2,:,:); %contra1: left sensors, right bp
%     contra(:,2,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,2),:,:,: ,:,:,1,:,:); % contra2: right sensors, left bp
%     ipsi(:,1,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,1),:,:,: ,:,:,1,:,:); %ipsi1
%     ipsi(:,2,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,2),:,:,: ,:,:,2,:,:); % ipsi2
%     respavglat(:,latctr,:,:,:, :,:,:,:) =  mean(contra,2) - mean(ipsi,2); % dims: subj_soi_bandoi_toi_npharm_nmotor_ndiff_stimorresp_npup_nRT
%     
%     latctr = latctr+1;
%     % wrt button press
%     contra(:,1,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,1),:,:,: ,:,:,:,2,:); %contra1: left sensors, target right
%     contra(:,2,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,2),:,:,: ,:,:,:,1,:); % contra2: right sensors, target left
%     ipsi(:,1,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,1),:,:,: ,:,:,:,1,:); %ipsi1
%     ipsi(:,2,:,:,: ,:,:,:,:) = respavg(:,sois(isoi,2),:,:,: ,:,:,:,2,:); % ipsi2
%     respavglat(:,latctr,:,:,:, :,:,:,:) =  mean(contra,2) - mean(ipsi,2); % dims: subj_soi_bandoi_toi_npharm_nmotor_ndiff_stimorresp_npup_nRT
% 
%     latctr = latctr+1; contra=[]; ipsi=[];
%     % wrt choice, depends on SR mapping
%     % pressipsi session
%     contra(:,1,:,:,: ,:,:,:) = respavg(:,sois(isoi,1),:,:,: ,1,:,:,2,:); %contra1: left sensors, choice right, pressipsi session
%     contra(:,2,:,:,: ,:,:,:) = respavg(:,sois(isoi,2),:,:,: ,1,:,:,1,:); % contra2: right sensors, choice left
%     ipsi(:,1,:,:,: ,:,:,:) =   respavg(:,sois(isoi,1),:,:,: ,1,:,:,1,:); %ipsi1
%     ipsi(:,2,:,:,: ,:,:,:) =   respavg(:,sois(isoi,2),:,:,: ,1,:,:,2,:); % ipsi2
%     temp(1,:,:,:,:, :,:,:) =  mean(contra,2) - mean(ipsi,2); 
%     % presscontra session
%     contra(:,1,:,:,: ,:,:,:) = respavg(:,sois(isoi,1),:,:,: ,2,:,:,1,:); %contra1: left sensors, choice right, presscontra session
%     contra(:,2,:,:,: ,:,:,:) = respavg(:,sois(isoi,2),:,:,: ,2,:,:,2,:); % contra2: right sensors, choice left
%     ipsi(:,1,:,:,: ,:,:,:) =   respavg(:,sois(isoi,1),:,:,: ,2,:,:,2,:); %ipsi1
%     ipsi(:,2,:,:,: ,:,:,:) =   respavg(:,sois(isoi,2),:,:,: ,2,:,:,1,:); % ipsi2
%     temp(2,:,:,:,:, :,:,:) =  mean(contra,2) - mean(ipsi,2); 
%     
%     respavglat(:,latctr,:,:,:, 3,:,:,:) =  mean(temp);
%     fprintf(' . ')
% end
% outfile = ['respavglat_' trigger];
% save(fullfile(PREIN, 'respavg', outfile), 'respavglat');
% 
% fprintf(' Done. \n')
% 


% % LATERALIZATION OCC WRT CHOICE
% % choice 1 - choice 2
% % ses1: stim1 and resp1
% % ses1: stim2 and resp2
% 
% % ses2: stim2 and resp1
% % ses2: stim1 and resp2
% respavg_size = size(respavg);
% temp = nan( respavg_size(1:end-2) ); % leave out stim and resp dims
% temp(:,:,:,:, :,1,:) = squeeze(respavg(:,:,:,:, :,1,:,3,1) - respavg(:,:,:,:, :,1,:,3,2)); % ses1 choice left-right, independent of stim
% temp(:,:,:,:, :,2,:) = squeeze(respavg(:,:,:,:, :,2,:,3,2) - respavg(:,:,:,:, :,2,:,3,1)); % ses2 choice left-right, independent of stim
% temp(:,:,:,:, :,3,:) = mean(temp(:,:,:,:, :,1:2,:), 6); % average across sessions
% % add lateralization as extra sensor
% respavg(:,269,:,:, :,:,:,3,3) = mean(temp(:,rightoccind,:,:, :,:,:), 2) - mean(temp(:,leftoccind,:,:, :,:,:), 2); 
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTchoice';
% sens.leg{3} = 'lateralization_occ_WRTchoice';
% sens.ind{3} = 269;
% 
% % LATERALIZATION OCC WRT button press
% % resp 1 -  2
% % respavg_size = size(respavg);
% % temp = nan( respavg_size(1:end-2) ); % leave out stim and resp dims
% temp = squeeze(respavg(:,:,:,:, :,:,:,3,1) - squeeze(respavg(:,:,:,:, :,:,:,3,2))); % button left-right, independent of stim
% % add lateralization as extra sensor
% respavg(:,270,:,:, :,:,:,3,3) = mean(temp(:,rightoccind,:,:, :,:,:), 2) - mean(temp(:,leftoccind,:,:, :,:,:), 2); 
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTbutton';
% sens.leg{4} = 'lateralization_occ_WRTbutton';
% sens.ind{4} = 270;
