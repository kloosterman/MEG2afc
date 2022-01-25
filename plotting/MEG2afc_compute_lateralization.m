function [ output_args ] = MEG2afc_compute_lateralization( input_args )
% MEG2afc_compute_lateralization
% compute lateralization

% Todo input args
NKsensorselection
leftind = {leftoccind leftmotorind};
rightind = {rightoccind rightmotorind};

latsoi = {'occipital', 'motor'};
latrWRT = {'wrtbp' 'wrttarget'};
latrleg = {'contra1', 'contra2', 'ipsi1', 'ipsi2', 'contra-ipsi' };
bandoi = [7 13; 15 30; 40 80];

fprintf('\nCompute variance acrross trial, collapse over diff, stim and choice . . .\n')
% dims subjlatrvar: var/mean isoi iwrt freq time drug regime
for isoi=1:2 % occipital and motor cortex
    for icol = [2,3]         % wrt: stim location (column 2) and bp (col 3)
        fprintf('.')
        for idiff = 1:2 % easy and hard
            trialind = find(  trialinfo(:,icol) == 2 & trialinfo(:,1) == idiff); % target/bp right
            contra = mean(powdat(trialind,leftind{isoi},:,:),2); %contra1: left sensors, right target
            ipsi = mean(powdat(trialind,rightind{isoi},:,:),2); %ipsi1: right sensors, right target
            latr = squeeze(contra - ipsi);
            
            trialind = find(  trialinfo(:,icol) == 1 & trialinfo(:,1) == idiff); % target/bp right
            contra = mean(powdat(trialind,rightind{isoi},:,:),2); %contra2: right sensors, left target
            ipsi = mean(powdat(trialind,leftind{isoi},:,:),2); %ipsi2: left sensors, left target
            latr = [latr; squeeze(contra - ipsi)]; %concat latr trials
            for iband = 1:length(bandoi)
                iband_frind = find(faxis >= bandoi(iband,1) & faxis <= bandoi(iband,2));
                latr(:,length(faxis)+iband,:) = mean(latr(:,iband_frind,:),2); %put bandoi in freq dim
            end
            subjlatrvar(1, isoi, icol-1, :,:, ipharm, imotor, idiff) = squeeze(nanvar(latr)); %compute variance over trials
            subjlatrvar(2, isoi, icol-1, :,:, ipharm, imotor, idiff) = squeeze(nanmean(latr)); %also add mean over trials
        end
        
        trialind = find(  trialinfo(:,icol) == 2 ); % target/bp right
        contra = mean(powdat(trialind,leftind{isoi},:,:),2); %contra1: left sensors, right target
        ipsi = mean(powdat(trialind,rightind{isoi},:,:),2); %ipsi1: right sensors, right target
        latr = squeeze(contra - ipsi);
        
        trialind = find(  trialinfo(:,icol) == 1 ); % target/bp left
        contra = mean(powdat(trialind,rightind{isoi},:,:),2); %contra2: right sensors, left target
        ipsi = mean(powdat(trialind,leftind{isoi},:,:),2); %ipsi2: left sensors, left target
        latr = [latr; squeeze(contra - ipsi)]; %concat latr trials
        
        for iband = 1:length(bandoi)
            iband_frind = find(faxis >= bandoi(iband,1) & faxis <= bandoi(iband,2));
            latr(:,length(faxis)+iband,:) = mean(latr(:,iband_frind,:),2); %put bandoi in freq dim
        end
        subjlatrvar(1, isoi, icol-1, :,:, ipharm, imotor, 3) = squeeze(nanvar(latr)); %compute variance over trials
        subjlatrvar(2, isoi, icol-1, :,:, ipharm, imotor, 3) = squeeze(nanmean(latr)); %also add mean over trials
        %             latr_mean = squeeze(nanmean(latr)); %compute variance over trials
    end
end


end

