function [bpm] =  MEG2afc_plot_heartbeats()

PREIN = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/preproczapline-plus/heartbeats';
cd(PREIN)
PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/ABC_project/MEG2afc_paper/behavdata';

drugconds = {'drug' 'plac'};
simonconds = {'contra' 'ipsi'};

SUBJ= [1:5, 7:9, 11:21]; % all
bpm = NaN(8, length(SUBJ), 2, 2);

for idrug=1:2
  for isimon = 1:2
    for isub = 1:length(SUBJ)
      fprintf('NK%d\n', SUBJ(isub))
      runlist = dir(sprintf('NK%d_*%s_%s*.mat', SUBJ(isub),  drugconds{idrug}, simonconds{isimon} ));
      for irun = 1:length(runlist)
        load(runlist(irun).name)
        
        if ~isfield(cfg_heartbeats, 'heartbeats' )
          disp('heartbeats field not found, skipping')
          continue
        end
          
        heartbeats = cfg_heartbeats.heartbeats;
        fprintf('%d beats      ', length(heartbeats))
        if length(heartbeats)< 300
          disp('< 300 heartbeats, skipping')
          continue
        end
        inter_beat_durs = diff(heartbeats(:,1)) / 350; % in sec, assuming downsample to 350 Hz
        inter_beat_durs = inter_beat_durs(zscore(inter_beat_durs) < 3); % remove outlier durs, indicates ECG is loose
        bpm(irun, isub, idrug, isimon) = length(inter_beat_durs) / (sum(inter_beat_durs)/60);

      end
      fprintf('\n')
    end
  end
end

bpm = nanmean(bpm);
bpm = squeeze(nanmean(bpm,4));

[h,p, ~, stats]=ttest(bpm(:,1), bpm(:,2))

bpm(:,1) - bpm(:,2)
mean(bpm)
figure; plot(bpm')

%% make table, export
subjid={};
for isub = 1:length(SUBJ)
  subjid{isub,1} = sprintf('NK%d', SUBJ(isub));
end
t = table(subjid, bpm(:,1), bpm(:,2), 'VariableNames', {'subjID', 'ATX' 'placebo'});
writetable(t, fullfile(PREOUT, 'heartbeats_MEGatx.csv'))

%% also write csv for pupil
load(fullfile(PREOUT, 'behavstruct.mat' ));
pupil = squeeze(behavior.basepupil(:,end,1:2,3))
figure; plot(pupil')

[h,p]=ttest(pupil(:,1), pupil(:,2))
mean(pupil)

t = table(subjid, pupil(:,1), pupil(:,2), 'VariableNames', {'subjID', 'ATX' 'placebo'});
writetable(t, fullfile(PREOUT, 'pupil_MEGatx.csv'))

