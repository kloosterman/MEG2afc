% function [ output_args ] = MEG2afc_load_respavg( trigger )
%MEG2afc_load_respavg Load subjrespavg for each subject, create condition
%labels, contrast to prepare for plotting etc. occ and motor poolings

% subjrespavg made by MEG2afc_concat_runs
% MEG hh conditions:
% drug placebo
% easy hard,
% pupil hi lo

close all
% clear all
% disp(trigger)

% ALL Subjects:
SUBJ  = { 
%     'NK1' ... has funny driftrate
    'NK2' 'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
         'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

nsub=length(SUBJ);

baseline = 'trial'; 
% analysistype = 'full';
% analysistype = 'low';
analysistype = {'low' 'full'};
baselinetype = 'respavg_baselinepersession'; % or respavg_baseline_acrosssessions

%%%set paths
if ismac
    MATLABPATH = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/';
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman/'; % on the cluster
else
    MATLABPATH = '/home/mpib/kloosterman/MATLAB/';
    basepath = '/home/mpib/kloosterman'; % on the cluster
end
if ~isdeployed
    addpath(genpath(fullfile(MATLABPATH, 'MEG_HH_analysis')))
    addpath(fullfile(MATLABPATH, 'tools', 'qsub_tardis'))
    addpath(fullfile(MATLABPATH, 'tools/fieldtrip-20150803')) %inc JJ edit ft_artifact_zvalue
    ft_defaults
end

FREQLO(1,1) = 2; %lowfreq
FREQHI(1,1) = 35;
FREQLO(2,1) = 37; % highfreq
FREQHI(2,1) = 149;

% keep N timepoints the same for stim and resp!
TIMLO(1,1) = -0.25; %stim
TIMHI(1,1) = 1; 
TIMLO(2,1) = -1;%resp
TIMHI(2,1) = 0.25; %dat{itrig}.TIMHI;

trigger_leg = {'stim', 'resp'};
freq_leg = {'low', 'full'};
for itrig = 1:2
    for ifreq = 1:2
        
        PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq', analysistype{ifreq}, trigger_leg{itrig});
        PREOUT = fullfile(basepath, 'plots'); mkdir(PREOUT)
        
        % example subjrespavg, put in dat
        subjrespavgin = sprintf('%s/%s/%s_%s.mat', PREIN, baselinetype, SUBJ{1}, baseline);
        fprintf('loading subrespavg %s\n', subjrespavgin)
        dat = load(subjrespavgin);
        
%         faxis{ifreq} = dat{ifreq}.faxis; % take faxis from data at it is now
        
        frind{ifreq} = find((dat.faxis >= FREQLO(ifreq,1)) & (dat.faxis <= FREQHI(ifreq,1))); % select frind of interest
        faxis{ifreq} = dat.faxis(frind{ifreq});

    end
    
    tind{itrig} = find((dat.taxis >= TIMLO(itrig,1)) & (dat.taxis <= TIMHI(itrig,1)));
    taxis{itrig} = dat.taxis(tind{itrig});
    

end

taxis_lengths = cellfun(@length, taxis);
faxis_lengths = cellfun(@length, faxis);
faxis_all = [faxis{1} faxis{2} ];
taxis_all = [taxis{1} taxis{2} ];
chlabel = dat.chlabel;

% 1     2    3      4        5          6       7       8       9         10          
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
% respavg = nan(length(SUBJ),6,length(faxis),length(taxis),2,2,3,3,3,3, 'single'); % only poolings, collapse over pupil and RT

% subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp) 
% respavg = nan(length(SUBJ),length(chlabel),length(faxis),length(taxis),2,2,3,3,3, 'single'); % only poolings, collapse over pupil and RT

% % Combine stim and resp-locked: take longest taxis for respavg
% respavg = nan(length(SUBJ),length(chlabel),length(faxis),max(taxis_lengths),4,3,4,3,3, 2, 'single'); % subj chan freq time ipharm, imotor, idiff, istim, iresp, itrig

% Combine freqlo and hi stim and resp-locked: take longest taxis for respavg+ get rid of
% stim and resp dims
respavg = nan(length(SUBJ),length(chlabel),max(faxis_lengths), max(taxis_lengths),4, 4, 2, 2, 'single'); % subj chan freq time ipharm, idiff, itrig, ifreq

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
    for itrig = 1:2
        for ifreq = 1:2
            
            PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq', analysistype{ifreq}, trigger_leg{itrig}, baselinetype, filesep );
            subjrespavgin = sprintf('%s_%s.mat', SUBJ{isub}, baseline);
            
            fprintf('loading subrespavg %s %s %s\n', subjrespavgin, trigger_leg{itrig}, freq_leg{ifreq})
            dat = load([PREIN subjrespavgin]);
            
            % MEG2afc_concat_runs structure:
            %     subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp, ipup, itrig)
            
            %         respavg(isub,:,:,tind{itrig}, 1:2,1:2,1:3,1:3,1:3, itrig) = dat.subjrespavg(:,frind,tind{itrig},:, :,:, :,:);
            % average over regime on the fly
            respavg(isub,:,1:length(frind{ifreq}),1:length(tind{itrig}), 1:2, 1:3, itrig, ifreq) = squeeze(nanmean(dat.subjrespavg(:,frind{ifreq}, tind{itrig},:, :, :, 3,3), 5));
            
            % TODO compute lateralizations here!
        end
    end
end

% contrasts
% 1     2    3      4        5          6       7       8       9         10
% subj nchan nfreq ntimebins pharma(2) diff(2) itrig  ifreq
%
respavg(:,:,:,:, 3,:,:,:) = nanmean(respavg(:,:,:,:, 1:2,:,:,:),5); % drug avg
respavg(:,:,:,:, 4,:,:,:) = respavg(:,:,:,:, 1,:,:,:) - respavg(:,:,:,:, 2,:,:,:); % drug contrast
respavg(:,:,:,:, :,3,:,:) = nanmean(respavg(:,:,:,:, :,1:2,:,:),6); % diff avg, overwrites diff-collapsed avg from concat_runs
respavg(:,:,:,:, :,4,:,:) = respavg(:,:,:,:, :,1,:,:) - respavg(:,:,:,:, :,2,:,:); % easy-hard

% respavg(:,:,:,:, :,:,:,4,:) = respavg(:,:,:,:, :,:,:,1,:) - respavg(:,:,:,:, :,:,:,2,:); % stim left-right
% respavg(:,:,:,:, :,:,:,:,4) = respavg(:,:,:,:, :,:,:,:,1) - respavg(:,:,:,:, :,:,:,:,2); % resp left-right

% [nsub, nchan,nfrq,ntim,npharm,nmotor,ndiff,nstim, nchoice,npup]=size(respavg)
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nchoice_npup_nRT_ncorrect'
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nresp_npup'
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nresp_ntrig'
respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_ntrig'
size(respavg)

XLIM = [TIMLO TIMHI];
YLIM = [FREQLO FREQHI];


%% LATERALIZATION   TODO add trigger dimension here!
% NKsensorselection

% % for itrig = 1:2
% pooling=[];
% % average over sensor groups FIRST.
% pooling(:,1,:,:, :,:,:,:,:) =  nanmean(respavg(:,leftoccind,:,:, :,:,:,:,:), 2);
% pooling(:,2,:,:, :,:,:,:,:) =  nanmean(respavg(:,rightoccind,:,:, :,:,:,:,:), 2);
% 
% pooling(:,3,:,:, :,:,:,:,:) =  nanmean(respavg(:,leftmotorind,:,:, :,:,:,:,:), 2);
% pooling(:,4,:,:, :,:,:,:,:) =  nanmean(respavg(:,rightmotorind,:,:, :,:,:,:,:), 2);
% 
% % LATERALIZATION OCC WRT CHOICE
% % choice 1 - choice 2
% % ses1 (ipsi): stim1 and resp1
% % ses1 (ipsi): stim2 and resp2
% 
% % ses2 (contra): stim2 and resp1
% % ses2 (contra): stim1 and resp2
% temp=[];
% % ipsi ses:
% % if choice is left, do pooling right -left
% temp(:,1,:,:, :,1,:,:) = pooling(:,2,:,:, :,1,:,:,1) - pooling(:,1,:,:, :,1,:,:,1);
% % if choice is right, do pooling left - right
% temp(:,2,:,:, :,1,:,:) = pooling(:,1,:,:, :,1,:,:,2) - pooling(:,2,:,:, :,1,:,:,2);
% 
% % same for the contra ses
% % if choice is left, do pooling right -left
% temp(:,1,:,:, :,2,:,:) = pooling(:,2,:,:, :,2,:,:,2) - pooling(:,1,:,:, :,2,:,:,2);
% % if choice is right, do pooling left - right
% temp(:,2,:,:, :,2,:,:) = pooling(:,1,:,:, :,2,:,:,1) - pooling(:,2,:,:, :,2,:,:,1);
% 
% respavg(:,269,:,:, :,1:2,:,:,3) = mean(temp,2); %average across pooling subtractions
% % average across motor mapping down below
% 
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTchoice';
% sens.leg{end+1} = 'lateralization_occ_WRTchoice';
% sens.ind{end+1} = 269;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % LATERALIZATION OCC WRT button press
% temp=[];
% % respright: left - right pooling 
% temp(:,1,:,:, :,:,:,:) = pooling(:,1,:,:, :,:,:,:,2) - pooling(:,2,:,:, :,:,:,:,2); % respright: left - right pooling 
% % respleft: right - left pooling 
% temp(:,2,:,:, :,:,:,:) = pooling(:,2,:,:, :,:,:,:,1) - pooling(:,1,:,:, :,:,:,:,1); % respleft: right - left pooling 
% respavg(:,270,:,:, :,:,:,:,3) = nanmean(temp,2); %average across pooling subtractions
% chlabel{end+1} = 'lateralization_occ_WRTbutton';
% sens.leg{end+1} = 'lateralization_occ_WRTbutton';
% sens.ind{end+1} = 270;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % LATERALIZATION OCC WRT stim location
% temp=[];
% % stimright: left - right pooling 
% temp(:,1,:,:, :,:,:,:) = pooling(:,1,:,:, :,:,:,2,:) - pooling(:,2,:,:, :,:,:,2,:); % stimright: left - right pooling 
% % stimleft: right - left pooling 
% temp(:,2,:,:, :,:,:,:) = pooling(:,2,:,:, :,:,:,1,:) - pooling(:,1,:,:, :,:,:,1,:); % stimleft: right - left pooling 
% 
% respavg(:,271,:,:, :,:,:,:,3) =  nanmean(temp,2); %average across pooling subtractions
% clear temp
% chlabel{end+1} = 'lateralization_occ_WRTtarget';
% sens.leg{end+1} = 'lateralization_occ_WRTtarget';
% sens.ind{end+1} = 271;
% 
% %% DO THE SAME FOR MOTOR CORTEX
% 
% % LATERALIZATION MOTOR WRT CHOICE
% % choice 1 - choice 2
% % ses1 (ipsi): stim1 and resp1
% % ses1 (ipsi): stim2 and resp2
% 
% % ses2 (contra): stim2 and resp1
% % ses2 (contra): stim1 and resp2
% temp=[];
% % ipsi ses:
% % if choice is left, do pooling right -left
% temp(:,1,:,:, :,1,:,:) = pooling(:,4,:,:, :,1,:,:,1) - pooling(:,3,:,:, :,1,:,:,1);
% % if choice is right, do pooling left - right
% temp(:,2,:,:, :,1,:,:) = pooling(:,3,:,:, :,1,:,:,2) - pooling(:,4,:,:, :,1,:,:,2);
% 
% % same for the contra ses
% % if choice is left, do pooling right -left
% temp(:,1,:,:, :,2,:,:) = pooling(:,4,:,:, :,2,:,:,2) - pooling(:,3,:,:, :,2,:,:,2);
% % if choice is right, do pooling left - right
% temp(:,2,:,:, :,2,:,:) = pooling(:,3,:,:, :,2,:,:,1) - pooling(:,4,:,:, :,2,:,:,1);
% 
% respavg(:,272,:,:, :,1:2,:,:,3) = nanmean(temp,2); %average across pooling subtractions
% % average across motor mapping down below
% 
% clear temp
% chlabel{end+1} = 'lateralization_motor_WRTchoice';
% sens.leg{end+1} = 'lateralization_motor_WRTchoice';
% sens.ind{end+1} = 272;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % LATERALIZATION MOTOR WRT button press
% temp=[];
% % respright: left - right pooling 
% temp(:,1,:,:, :,:,:,:) = pooling(:,3,:,:, :,:,:,:,2) - pooling(:,4,:,:, :,:,:,:,2); % respright: left - right pooling 
% % respleft: right - left pooling 
% temp(:,2,:,:, :,:,:,:) = pooling(:,4,:,:, :,:,:,:,1) - pooling(:,3,:,:, :,:,:,:,1); % respleft: right - left pooling 
% respavg(:,273,:,:, :,:,:,:,3) = nanmean(temp,2); %average across pooling subtractions
% chlabel{end+1} = 'lateralization_motor_WRTbutton';
% sens.leg{end+1} = 'lateralization_motor_WRTbutton';
% sens.ind{end+1} = 273;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % LATERALIZATION MOTOR WRT stim location
% temp=[];
% % stimright: left - right pooling 
% temp(:,1,:,:, :,:,:,:) = pooling(:,3,:,:, :,:,:,2,:) - pooling(:,4,:,:, :,:,:,2,:); % stimright: left - right pooling 
% % stimleft: right - left pooling 
% temp(:,2,:,:, :,:,:,:) = pooling(:,4,:,:, :,:,:,1,:) - pooling(:,3,:,:, :,:,:,1,:); % stimleft: right - left pooling 
% 
% respavg(:,274,:,:, :,:,:,:,3) =  nanmean(temp,2); %average across pooling subtractions
% clear temp
% chlabel{end+1} = 'lateralization_motor_WRTtarget';
% sens.leg{end+1} = 'lateralization_motor_WRTtarget';
% sens.ind{end+1} = 274;





















% OLD?


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
