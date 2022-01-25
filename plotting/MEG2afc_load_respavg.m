function [respavg] = MEG2afc_load_respavg( )
%MEG2afc_load_respavg Load subjrespavg for each subject, create condition
%labels, contrast to prepare for plotting etc. occ and motor poolings

% subjrespavg made by MEG2afc_concat_runs
% MEG hh conditions:
% drug placebo
% easy hard,
% pupil hi lo

% close all
% clear all
% disp(trigger)

% ALL Subjects:
SUBJ  = { 
    'NK1' ... has funny driftrate
    'NK2' ...
    'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
         'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

%   % ALL Subjects minus outlier beta band NK8:
% SUBJ  = { 
%     'NK1' ... has funny driftrate
%     'NK2' ...
%     'NK3'   'NK4'   'NK5'   'NK7'   'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
%          'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
%     }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

%   % only subjects who get worse from the drug:
% SUBJ  = { 
%     'NK1' ... 'NK2'    'NK3' ... has funny driftrate
%      'NK4'   'NK5'  'NK8'  'NK9'   'NK11'  'NK12'    'NK14'  'NK15'...    'NK7'   'NK13' 
%          'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
%     };  
% SUBJ([2 3 6 11]) ARE PROBABLY the ones to be dropped: 'NK2'    'NK3'    'NK7'    'NK13'
  
nsub=length(SUBJ);

respavg = [];
respavg.SUBJ = SUBJ;
respavg.PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/plots';

baseline = 'trial'; 
analysistype = {'low' 'full' 'lowfine' 'highfine'};
baselinetype = 'respavg_baselinepersession'; % or respavg_baseline_acrosssessions
% phaselock_type = 'totalpow';
phaselock_type = 'induced';

if ismac
%     basepath = '/Users/kloosterman/gridmaster2012/kloosterman'; % on the cluster
    basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc'; % on the cluster
%     basepath = '/Volumes/LNDG/user/Niels/MEG2afc'; % on the server
else
    basepath = '/home/mpib/kloosterman'; % on the cluster
end

FREQLO(1,1) = 2; %lowfreq
FREQHI(1,1) = 35;
% FREQLO(2,1) = 37; % highfreq
FREQLO(2,1) = 5; % highfreq
FREQHI(2,1) = 149;
FREQLO(3,1) = 1; % lowfreq fine res
FREQHI(3,1) = 36;
FREQLO(4,1) = 37; % highfreq
FREQHI(4,1) = 149;

% TIMLO(1,1) = -0.25; %stim
% TIMHI(1,1) = 1.5; 
% TIMLO(2,1) = -1.5;%resp
% TIMHI(2,1) = 0.25; %dat{itrig}.TIMHI;

TIMLO(1,1) = -0.5; %stim
TIMHI(1,1) = 0.4; 
TIMLO(2,1) = -0.4;%resp
TIMHI(2,1) = 0.3; %dat{itrig}.TIMHI;
TIMLO(3,1) = -0.4;%only -0.3 present
TIMHI(3,1) = 0.3; 
TIMLO(4,1) = -0.4;%only -0.3 present
TIMHI(4,1) = 0.3; 

% % keep N timepoints the same for stim and resp!
% TIMLO(1,1) = -0.2; %stim
% TIMHI(1,1) = 0.4; 
% TIMLO(2,1) = -0.4;%resp
% TIMHI(2,1) = 0.2; %dat{itrig}.TIMHI;
% 
trigger_leg = {'stim', 'resp', 'baseline'};
freq_leg = {'low', 'full', 'lowfine', 'highfine'};

% only load baseline freq data
trig_oi = 3;
band_oi = 3%:4;
% % load stim, resp, low high
% trig_oi = 1:3;
% band_oi = 1;
% trig_oi = 3;
% band_oi = 3;

for itrig = trig_oi
    for iband = band_oi
        
        PREIN = fullfile(basepath, 'freq', analysistype{iband}, trigger_leg{itrig});
%         PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq_robdetr', analysistype{iband}, trigger_leg{itrig});
        PREOUT = fullfile(basepath, 'plots'); mkdir(PREOUT)
        
        % example subjrespavg, put in dat
        subjrespavgin = sprintf('%s/%s/%s/%s_%s.mat', PREIN, baselinetype, phaselock_type, SUBJ{1}, baseline);
        fprintf('loading subrespavg %s\n', subjrespavgin)
        dat = load(subjrespavgin);
        
%         faxis{iband} = dat{iband}.faxis; % take faxis from data at it is now
        
        frind{iband} = find((dat.faxis >= FREQLO(iband,1)) & (dat.faxis <= FREQHI(iband,1))); % select frind of interest
        faxis{iband} = dat.faxis(frind{iband});

    end
    
    tind{itrig} = find((dat.taxis >= TIMLO(itrig,1)) & (dat.taxis <= TIMHI(itrig,1)));
    taxis{itrig} = dat.taxis(tind{itrig});
    
end

taxis_lengths = cellfun(@length, taxis);
faxis_lengths = cellfun(@length, faxis);
faxis_all = cat(2,faxis{:});
taxis_all = cat(2,taxis{:});

respavg.label = dat.chlabel;
respavg.time = taxis;
respavg.freq = faxis;
respavg.trigger_leg = trigger_leg;
respavg.freq_leg = freq_leg;

% Condition labels
%--------------------------------------------------------------------------
respavg.pharm_conds = {'atomox' 'placebo' 'drugscomb' 'atomox - placebo'};
respavg.motor_conds = {'pressipsi' 'presscontra' 'regimecomb' 'pressipsi - presscontra'};
respavg.diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};
respavg.pupilconds = {'pupilhigh' 'pupillow' 'pupilcomb' 'pupilhi-lo'};
respavg.stim_conds = {'stimleft' 'stimright' 'allstim' 'stimleft-right'};
respavg.choice_conds = {'choiceleft' 'choiceright' 'allchoice' 'choiceleft-right'};
respavg.rt_conds = {'slow' 'medium' 'fast' 'allrts'};
respavg.correct_conds = {'correct' 'error' 'corr+err'};
respavg.sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
respavg.latr_leg = {'modulation', 'latr wrt resp', 'latr wrt choice', 'latr wrt stim'};

% % % return
if isdeployed || ~ismac
    fprintf('We are on the cluster, not loading respavg\n\n\n')
    return
end
%% load in single subjects



% 1     2    3      4        5          6       7       8       9         10          
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
% respavg = nan(length(SUBJ),6,length(faxis),length(taxis),2,2,3,3,3,3, 'single'); % only poolings, collapse over pupil and RT

% subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp) 
% respavg = nan(length(SUBJ),length(chlabel),length(faxis),length(taxis),2,2,3,3,3, 'single'); % only poolings, collapse over pupil and RT

% % Combine stim and resp-locked: take longest taxis for respavg
% respavg = nan(length(SUBJ),length(chlabel),length(faxis),max(taxis_lengths),4,3,4,3,3, 2, 'single'); % subj chan freq time ipharm, imotor, idiff, istim, iresp, itrig

% Combine freqlo and hi stim and resp-locked: take longest taxis for respavg+ get rid of
% stim and resp dims
respavg.dat = nan(length(SUBJ),length(respavg.label),max(faxis_lengths), max(taxis_lengths),4, 4, 2, 3, 3, 'single'); % subj chan freq time ipharm, idiff, itrig, iband, ilatr
% respavg.pow = nan(length(SUBJ),length(respavg.label),max(faxis_lengths), max(taxis_lengths),4, 4, 2, 3, 3, 'double'); % subj chan freq time ipharm, idiff, itrig, iband, ilatr
% ilatr: modulation; lateralization wrt resp
sens = MEG2afc_sensorselection();
LR_subtract_mat = sens.LR_subtract_mat;

for isub = 1:length(respavg.SUBJ)
%     for itrig = 3 %3%:2
%         for iband = 3:4 %1:3
    for itrig = trig_oi
        for iband = band_oi

%             PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq_robdetr', analysistype{iband}, trigger_leg{itrig}, baselinetype, phaselock_type, filesep );
            PREIN = fullfile(basepath, 'freq', analysistype{iband}, trigger_leg{itrig}, baselinetype, phaselock_type, filesep );
            subjrespavgin = sprintf('%s_%s.mat', SUBJ{isub}, baseline);
            fprintf('loading subrespavg %s %s %s\n', subjrespavgin, trigger_leg{itrig}, freq_leg{iband})
            
            dat = load([PREIN subjrespavgin]);
            
            % raw pow:
            respavg.pow(isub,:,1:length(frind{iband}),1:length(tind{itrig}), 1:2, 1:3, itrig, iband, 1) = ...
                squeeze(nanmean(dat.subjpowavg(:,frind{iband}, tind{itrig},:, :, :, 3,3), 5)); % 1 = modulation, avg over ses

            if itrig == 3
                continue
            end
            % average over regime on the fly
            respavg.dat(isub,:,1:length(frind{iband}),1:length(tind{itrig}), 1:2, 1:3, itrig, iband, 1) = ...
                squeeze(nanmean(dat.subjrespavg(:,frind{iband}, tind{itrig},:, :, :, 3,3), 5)); % 1 = modulation, avg over ses
            
            % lateralization wrt resp:
            sensordat = squeeze(nanmean(dat.subjrespavg(:,frind{iband}, tind{itrig},:, :, :, 3, 1:2), 5)); % keep left and right bp, 6 dims remain
            temp = nan( size(sensordat) );% last dim will be used for subtractions, resp not needed
            % respright: left - right sensordat
            temp(LR_subtract_mat(:,2),:,:, :,:,1) = sensordat(LR_subtract_mat(:,1),:,:, :,:,2) - sensordat(LR_subtract_mat(:,2),:,:, :,:,2); % respright: left - right sensordat
            % respleft: right - left sensordat
            temp(LR_subtract_mat(:,2),:,:, :,:,2) = sensordat(LR_subtract_mat(:,2),:,:, :,:,1) - sensordat(LR_subtract_mat(:,1),:,:, :,:,1); % respleft: right - left sensordat
            temp = squeeze(nanmean(temp,6)); % average across 2 lateralizations - project on RH
            
            respavg.dat(isub,:,1:length(frind{iband}),1:length(tind{itrig}), 1:2, 1:3, itrig, iband, 2) = temp; % L channels are nan
            
            %%%
            % lateralization wrt choice:
            % dims subjrespavg(:,:,:, ipharm, imotor, idiff, istim, iresp)
            sensordat = squeeze( dat.subjrespavg(:,frind{iband}, tind{itrig},:, :, :, 3, 1:2) ); % dims: chan freq time drug ses diff resp
            temp = nan( size(sensordat) );% last dim will be used for latr subtractions wrt choice, resp not needed

            % ipsi ses:
            % if choice is left, do pooling right -left, project on RH
            temp(LR_subtract_mat(:,2),:,:, :,1,:,1) = sensordat(LR_subtract_mat(:,2),:,:, :,1,:,1) - sensordat(LR_subtract_mat(:,1),:,:, :,1,:,1); 
            % if choice is right, do pooling left - right
            temp(LR_subtract_mat(:,2),:,:, :,1,:,2) = sensordat(LR_subtract_mat(:,1),:,:, :,1,:,2) - sensordat(LR_subtract_mat(:,2),:,:, :,1,:,2); 
            
            % same for the contra ses
            % if choice is left, do pooling right -left
            temp(LR_subtract_mat(:,2),:,:, :,2,:,1) = sensordat(LR_subtract_mat(:,2),:,:, :,2,:,2) - sensordat(LR_subtract_mat(:,1),:,:, :,2,:,2); 
            % if choice is right, do pooling left - right
            temp(LR_subtract_mat(:,2),:,:, :,2,:,2) = sensordat(LR_subtract_mat(:,1),:,:, :,2,:,1) - sensordat(LR_subtract_mat(:,2),:,:, :,2,:,1); 
            temp = nanmean(temp,7); % average across 2 lateralizations
            temp = squeeze(nanmean(temp,5)); % average across 2 sessions
            
            respavg.dat(isub,:,1:length(frind{iband}),1:length(tind{itrig}), 1:2, 1:3, itrig, iband, 3) = temp; % L channels are nan
           
            %%%
            % lateralization wrt stimulus location:
            sensordat = squeeze(nanmean(dat.subjrespavg(:,frind{iband}, tind{itrig},:, :, :, 1:2, 3), 5)); % keep left and right stim, 6 dims remain
            temp = nan( size(sensordat) );% last dim will be used for subtractions, stim dim not needed
            % respright: left - right sensordat
            temp(LR_subtract_mat(:,2),:,:, :,:,1) = sensordat(LR_subtract_mat(:,1),:,:, :,:,2) - sensordat(LR_subtract_mat(:,2),:,:, :,:,2); % respright: left - right sensordat
            % respleft: right - left sensordat
            temp(LR_subtract_mat(:,2),:,:, :,:,2) = sensordat(LR_subtract_mat(:,2),:,:, :,:,1) - sensordat(LR_subtract_mat(:,1),:,:, :,:,1); % respleft: right - left sensordat
            temp = squeeze(nanmean(temp,6)); % average across 2 lateralizations - project on RH
            
            respavg.dat(isub, :, 1:length(frind{iband}), 1:length(tind{itrig}), 1:2, 1:3, itrig, iband, 4) = temp; % L channels are nan
            
        end
    end
end

% contrasts
% 1     2    3      4        5          6       7       8       9         10
% subj nchan nfreq ntimebins pharma(2) diff(2) itrig  iband
fprintf('Calculating averages and contrasts . . .\n')
respavg.dat(:,:,:,:, 3,:,:,:,:) = nanmean(respavg.dat(:,:,:,:, 1:2,:,:,:,:),5); % drug avg
respavg.dat(:,:,:,:, 4,:,:,:,:) = respavg.dat(:,:,:,:, 1,:,:,:,:) - respavg.dat(:,:,:,:, 2,:,:,:,:); % drug contrast
respavg.dat(:,:,:,:, :,3,:,:,:) = nanmean(respavg.dat(:,:,:,:, :,1:2,:,:,:),6); % diff avg, overwrites diff-collapsed avg from concat_runs
respavg.dat(:,:,:,:, :,4,:,:,:) = respavg.dat(:,:,:,:, :,1,:,:,:) - respavg.dat(:,:,:,:, :,2,:,:,:); % easy-hard

% respavg.pow(:,:,:,:, 4,:,:,:,:) = respavg.pow(:,:,:,:, 1,:,:,:,:) - respavg.pow(:,:,:,:, 2,:,:,:,:); % drug contrast

respavg.pow(:,:,:,:, 3,:,:,:,:) = mean(respavg.pow(:,:,:,:, 1:2,:,:,:,:), 5);
% respavg.pow(:,:,:,:, 4,:,:,:,:) = (respavg.pow(:,:,:,:, 1,:,:,:,:) - respavg.pow(:,:,:,:, 2,:,:,:,:)) ./ ...
%     respavg.pow(:,:,:,:, 3,:,:,:,:); % drug contrast
% % respavg.pow(:,:,:,:, 4,:,:,:,:) = respavg.pow(:,:,:,:, 4,:,:,:,:) .* 100; % convert to psc

respavg.pow(:,:,:,:, 4,:,:,:,:) = (respavg.pow(:,:,:,:, 1,:,:,:,:) - respavg.pow(:,:,:,:, 2,:,:,:,:)); % drug contrast

% respavg.pow(:,:,:,:, 4,:,:,:,:) = log(respavg.pow(:,:,:,:, 1,:,:,:,:)) - log(respavg.pow(:,:,:,:, 2,:,:,:,:)) 

respavg.dimord = 'subj_chan_freq_time_drug_diff_trig_freqrange_latr'; % subj chan freq time ipharm, idiff, itrig, iband, ilatr
respavg.datsize = size(respavg.dat);
respavg.sens = MEG2afc_sensorselection();

% respavg.driftrates = MEG2afc_load_drifts_regression(); % NK1 is out
respavg.ddmdat = MEG2afc_load_ddm_paras('/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_drug_diff2.csv');
% respavg.ddmdat.v_hard_atx - respavg.ddmdat.v_hard_plac

XLIM = [TIMLO TIMHI];
YLIM = [FREQLO FREQHI];

