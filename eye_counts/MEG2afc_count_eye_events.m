function eye_count = MEG2afc_count_eye_events(SAV)
% eye_count.dat = MEG2afc_count_eye_events counts blinks etc and puts them into
% a matrix
if nargin==0
    SAV=1;
end

addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'))
addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/cbrewer')

basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/preproc';

% ALL Subjects:
SUBJ  = { 
%     'NK1' ... has funny driftrate
    'NK2' ...
    'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'   'NK14'  'NK15'...  
         'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

nsub=length(SUBJ);

sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra' 'drug' 'placebo'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj

eye_count = [];
eye_count.SUBJ = SUBJ;
eye_count.ses = sesdirs;
eye_count.eventnames = {'threshold'      'eog_ver'      'eog_hor'         'jump'       'muscle'};
eye_count.dat = nan(nsub,4,8,5);

for isub = 1:nsub
    for ises = 1:4
        fprintf('Subj %s Ses %s  ', SUBJ{isub}, sesdirs{ises})
        try
            cd(fullfile(basepath, SUBJ{isub}, sesdirs{ises}))
%             fprintf(' found!\n')
        catch
            fprintf('not found\n')
            continue
        end
        
        runlist = dir('*preprocinfo.mat');
        nruns = length(runlist);
        for irun = 1:nruns
            load(runlist(irun).name);
            fprintf('%d  ' ,irun)

            eye_count.dat(isub, ises, irun,1) = size(preprocinfo.artfctdef.threshold.artifact,1);
            eye_count.dat(isub, ises, irun,2) = size(preprocinfo.artfctdef.eog_ver.artifact,1);
            eye_count.dat(isub, ises, irun,3) = size(preprocinfo.artfctdef.eog_hor.artifact,1);
            eye_count.dat(isub, ises, irun,4) = size(preprocinfo.artfctdef.jump.artifact,1);
            eye_count.dat(isub, ises, irun,5) = size(preprocinfo.artfctdef.muscle.artifact,1);

        end
        
%         %
%         nblinks.drug(isub,:) = sum(eye_count.dat(isub,5,2));
%         nblinks.plac(isub,:) = sum(eye_count.dat(:,6,2));

        % convert to n events per run        
        eye_count.dat(isub, ises, :,:) = eye_count.dat(isub, ises, :,:) / nruns;
         
        fprintf('\n')

    end
end

% sum over runs
eye_count.dat = squeeze(nansum(eye_count.dat, 3));


eye_count.dat(:,5,:) = sum(eye_count.dat(:,1:2,:), 2); %drug
eye_count.dat(:,6,:) = sum(eye_count.dat(:,3:4,:), 2); % plac

nblinks.drug = sum(eye_count.dat(:,5,2));
nblinks.plac = sum(eye_count.dat(:,6,2));

%% Plot eye events for drug vs plac
close all
figure;
set(gcf, 'Position', [100 150 900 900])

cmap = cbrewer('qual', 'Set1', 2);

for ievent=1:length(eye_count.eventnames)
    subplot(2,3,ievent)
    dat = eye_count.dat(:,5:6,ievent);
    barweb(mean(dat), std(dat)/sqrt(nsub))
    [~,p]= ttest(dat(:,1), dat(:,2));
    title(sprintf('%s, p = %g',  eye_count.eventnames{ievent}, p))
    if ievent ==1; legend(sesdirs{5:6}); end % , 'Location', 'southeast'
    ylabel('# events / run')
end
colormap(cmap(1:2,:))    

if SAV
    outfile = 'eye_events';    
    %                         print('-dpdf', fullfile(outpath, outfile))
    %                         export_fig( fullfile(outpath, outfile), '-pdf')
    export_fig( fullfile(basepath, outfile), '-png')
    cd(basepath)
end


