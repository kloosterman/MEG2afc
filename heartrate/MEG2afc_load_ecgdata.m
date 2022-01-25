% function [ output_args ] = MEG2afc_load_dvadata( trigger )
%MEG2afc_load_dvadata Load dva data for each subject, create condition
%labels, contrast to prepare for plotting etc. 

% subjrespavg made by MEG2afc_concat_runs
% MEG hh conditions:
% drug placebo
% easy hard,

close all

SUBJ  = { 
     'NK1' 'NK2' ...
    'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
         'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist, 'NK2' incomplete  'NK12' 'NK15' have bad ddm fits

drug_by_order_mat= [2	1	1	2;  %1 == drug, 2 = =placebo
1	2	1	2; % Fabian, only 3 sessions
1	1	2	2;
1	2	2	1;
1	2	2	1;
2	1	2	1;
2	2	1	1;
2	1	1	2;
2	2	1	1;
2	1	2	1;
1	1	2	2;
2	1	2	1;

1   2   2   1; %NK 15

2	1	1	2;
2	1	2	1;
1	2	2	1;
2	1	2	1;
1	2	2	1;
1	2	1	2]; 
% drug_by_order_mat(:,5) = 1
% drug_by_order_mat = logical(drug_by_order_mat);

motor_by_order_mat = [1	1	2	2  % 1==ipsi, 2==contra
1   2   2   1    
1	2	2	1
1	1	2	2
1	2	1	2
2	1	1	2
2	1	1	2
2	2	1	1
1	2	2	1
1	2	2	1
2	1	1	2
1	1	2	2

1	1	2	2 %NK 15

1	2	1	2
2	2	1	1
1	2	1	2
1	2	2	1
2	2	1	1
2	2	1	1];
% motor_by_order_mat = logical(motor_by_order_mat);

nsub=length(SUBJ);

baseline = 'trial'; analysistype = 'full';

%%%set paths
if ismac
    MATLABPATH = '/Users/kloosterman/gridmaster2012/MATLAB/';
    basepath = '/Users/kloosterman/gridmaster2012/'; % on the cluster
else
    MATLABPATH = '/home/mpib/kloosterman/MATLAB/';
    basepath = '/home/mpib/kloosterman'; % on the cluster
end
addpath(genpath(fullfile(MATLABPATH, 'MEG_HH_analysis')))
addpath(fullfile(MATLABPATH, 'toolbox', 'qsub_tardis'))
if exist('ft_defaults') ~= 2
    addpath(fullfile(MATLABPATH, 'toolbox/fieldtrip-20150803')) %inc JJ edit ft_artifact_zvalue
    ft_defaults
end

PREIN = fullfile(basepath, 'projectdata/MEG2afc/preproc');
PREOUT = fullfile(basepath, 'plots'); mkdir(PREOUT)

% example subjrespavg, put in dat
% % subjrespavgin = sprintf('%s/%s_dva.mat', PREIN, SUBJ{1});
% subjrespavgin = fullfile(PREIN, SUBJ{1}, 'drug_contra', '_dva.mat');
% fprintf('loading subrespavg %s\n', subjrespavgin)
% dat = load(subjrespavgin);

% set up arrays for analysis
% if strcmp(trigger, 'stim') % this timewindow is cut out from the subjrespavg timewin
% %     TIMLO          = -1;    TIMHI          =  1.5; % TFR's
%     TIMLO          = -1.5;    TIMHI          =  4; 
% else
% %     TIMLO          = -0.4 ;   TIMHI          =  0.2;
%     TIMLO          = -1.5 ;   TIMHI          = 2;
% end
% tind = find((dat.taxis>=TIMLO) & (dat.taxis<=TIMHI));
% taxis = dat.taxis(tind);
% chlabel = dat.chlabel;

% % respavg = nan(length(SUBJ),length(chlabel),length(faxis),length(taxis),4,3,3,3,3, 'single'); % only poolings, collapse over pupil and RT
% % 1     2    3      4        5          6       7       8       9       10           11     12
respavg = nan(nsub, 8, 2,2 ); % across var, across norm, within var, within norm, baseline-corrected within var
% respavg_accum = nan(nsub, 2, 3, 2,2,3,3,3, 'single' ); % across var, across norm, within var, within norm
% basedatavg = nan([nsub size(dat.basedat)]);
% 
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
    for ises = 1:4
        
        if ises==1, ipharm = 1;   imotor=1;
        elseif ises==2, ipharm = 1;   imotor=2;
        elseif ises==3, ipharm = 2;   imotor=1;
        elseif ises==4, ipharm = 2;   imotor=2;
        end
        
        try
            cd(fullfile(PREIN, SUBJ{isub}, sesdirs{ises}))
        catch
            warning('Name is nonexistent or not a directory')
        end        
        runlist = dir('*heartbeats.mat');
        
        for irun = 1:length(runlist)
            fprintf('loading heartbeats subj %s, run %d\n', SUBJ{isub}, irun)
            
            load(runlist(irun).name)
            
            run_dur = (artifact(end) - artifact(1)) / 1200 / 60;
            bpm = size(artifact,1) / run_dur;
            
            if bpm < 40
                bpm = nan; % channel broken er something
            end
            
            respavg(isub,irun, ipharm, imotor ) = bpm;
            
        end
        
    end
end

% [nsub, nchan,nfrq,ntim,npharm,nmotor,ndiff,nstim, nchoice,npup]=size(respavg)
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nchoice_npup_nRT_ncorrect'
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_nstim_nresp_npup'

% 1     2    3      4        5          6       7       8       9         10
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)
%
respavg(:,:,:,3) = nanmean(respavg,4); % motor avg
respavg(:,:,3,:) = nanmean(respavg,3); % drug avg
respavg(:,9,:,:) = nanmean(respavg,2); % run avg

% respavg(:,:,:,:, :,3,:,:,:) = nanmean(respavg(:,:,:,:, :,1:2,:,:,:),6);
% respavg(:,:,:,:, 3,:,:,:,:) = nanmean(respavg(:,:,:,:, 1:2,:,:,:,:),5);
% respavg(:,:,:,:, 4,:,:,:,:) = respavg(:,:,:,:, 1,:,:,:,:) - respavg(:,:,:,:, 2,:,:,:,:); % drug contrast
% % respavg(:,:,:,:, :,4,:,:,:) = respavg(:,:,:,:, :,1,:,:,:) - respavg(:,:,:,:, :,2,:,:,:); % motor ipsi - contra
% % respavg(:,:,:,:, :,:,4,:,:) = respavg(:,:,:,:, :,:,1,:,:) - respavg(:,:,:,:, :,:,2,:,:); % easy-hard
% % respavg(:,:,:,:, :,:,:,4,:) = respavg(:,:,:,:, :,:,:,1,:) - respavg(:,:,:,:, :,:,:,2,:); % stim left-right
% % respavg(:,:,:,:, :,:,:,:,4) = respavg(:,:,:,:, :,:,:,:,1) - respavg(:,:,:,:, :,:,:,:,2); % resp left-right
% 
% respavg_accum(:,:,:, :,3,:,:,:) =  mean(respavg_accum(:,:,:, :,1:2,:,:,:),5);

respavgdimord = 'nsub_nruns_npharm_nmotor'
respavg_size = size(respavg)
% respavg_accum_size = size(respavg_accum);
% basedatavg_size = size(basedatavg);


%% Make session # respavg, ie look up what condition was when
respavg_ses = nan(nsub, 9, 4); 
% respavg_ses_accum = nan(respavg_accum_size([1:4, 6:8 ])); 
% basedatavg_ses = nan(basedatavg_size([1:3, 5:7 ])); 
for isub=1:length(SUBJ)
    for ises = 1:4
        drug_ind = drug_by_order_mat(isub,ises);
        motor_ind = motor_by_order_mat(isub,ises);
        respavg_ses(isub,:,ises) = respavg(isub, :, drug_ind, motor_ind);
        
%         respavg_ses_accum(isub,:,:, ises,:,:,:) = respavg_accum(isub,:,:,drug_ind, motor_ind,:,:,:);
%         
%         basedatavg_ses(isub,:, ises,:,:,:) = basedatavg(isub,:, drug_ind, motor_ind,:,:,:);
        
%         if ~isempty(find(isnan(respavg(isub,1,3,1,drug_ind, motor_ind,3,3,3))))
%             disp(ises)
%             disp(SUBJ{isub})
%         end
    end
    
% %     % silence session 3
%     drug_ind = drug_by_order_mat(isub,2);
%     motor_ind = motor_by_order_mat(isub,2);
%     respavg(isub,:,:,:, drug_ind, motor_ind,:,:,:) = nan;
% 
%     %     % silence session 3
%     drug_ind = drug_by_order_mat(isub,3);
%     motor_ind = motor_by_order_mat(isub,3);
%     respavg(isub,:,:,:, drug_ind, motor_ind,:,:,:) = nan;
%     
%     % silence session 4
%     drug_ind = drug_by_order_mat(isub,4);
%     motor_ind = motor_by_order_mat(isub,4);
%     respavg(isub,:,:,:, drug_ind, motor_ind,:,:,:) = nan;
end
% % repeat contrasts
% respavg(:,:,:,:, :,3,:,:,:) = nanmean(respavg(:,:,:,:, :,1:2,:,:,:),6);
% respavg(:,:,:,:, 3,:,:,:,:) = nanmean(respavg(:,:,:,:, 1:2,:,:,:,:),5);
% respavg(:,:,:,:, 4,:,:,:,:) = respavg(:,:,:,:, 1,:,:,:,:) - respavg(:,:,:,:, 2,:,:,:,:); % drug contrast
% 
% % basedatavg_ses(basedatavg_ses == 0) = nan; %make NK2 ses 4 nan 
% %    
% % respavg_ses_accum(:,3,:, :,:,:,:) = basedatavg_ses; % put in accum to plot easier

