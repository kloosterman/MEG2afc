function MEG2afc_preproc_inspect()

close all
if ismac
  basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
end

cd(fullfile(basepath, 'preproc'))

runlist = dir('NK*.mat');
nruns = length(runlist);
runlist = runlist(randsample(nruns,nruns));

for irun = 1:nruns
  disp(runlist(irun).name)
  load(runlist(irun).name)
  cfg=[];
  cfg.layout = 'CTF275.lay';
  cfg.metric = 'var';
  cfg.channel = 'MEG';
  ft_rejectvisual(cfg, data)
  
%   ft_databrowser([],data)
end

% 2 chans weird in trial 38: load NK14_drug_contra_run4.mat
% 2e-12
% max 2.5e-12

