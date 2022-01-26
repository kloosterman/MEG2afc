function [bpm] =  MEG2afc_plot_heartbeats()

PREIN = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/preproczapline-plus/heartbeats';
cd(PREIN)

drugconds = {'drug' 'plac'};
simonconds = {'contra' 'ipsi'};

SUBJ= [1:5, 7:9, 11:21]; % all
bpm = NaN(16, length(SUBJ), 2, 2);

for idrug=1:2
  for isimon = 1:2
    for isub = 1:length(SUBJ)
      runlist = dir(sprintf('NK%d_*%s_%s*.mat', SUBJ(isub),  drugconds{idrug}, simonconds{isimon} ));
      for irun = 1:length(runlist)
        load(runlist(irun).name)
        bpm(irun, isub, idrug, isimon) = cfg_heartbeats.bpm;
      end
    end
  end
end

bpm = squeeze(nanmean(bpm));

[h,p]=ttest(bpm(:,1), bpm(:,2))

bpm(:,1) - bpm(:,2)
figure; plot(bpm')

