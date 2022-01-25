function [megdat] = MEG2afc_mergedva_stats(megdat)
% Run stats, add them to the megdat struct

close all
disp 'stats'

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
% cfg.alpha            = 0.1;
cfg.numrandomization = 1000;
cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = 'ctf275_neighb.mat';
cfg.neighbours       = []; %ft_prepare_neighbours(cfg0_neighb);

nsub = size(megdat.dvalock(3,3).avg,1);
design = zeros(2,2*nsub);
for i = 1:nsub
  design(1,i) = i;
end
for i = 1:nsub
  design(1,nsub+i) = i;
end
design(2,1:nsub)        = 1;
design(2,nsub+1:2*nsub) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;


for idva = 1:2
  for ichan = 1:2 % dva, norm
    cfg.channel = ichan;
    
    dva = megdat.dvalock(idva,3); % dvalock{idva, idiff}
    dvazero = dva; %create zero dva to test against
    dvazero.avg = zeros(size(dva.avg));
    stat= ft_timelockstatistics(cfg, dva, dvazero);
    
    %     stat{ifreq,itrig}.powspctrm = squeeze(mean(dva.powspctrm));
    %     stat{ifreq,itrig}.powspctrm = dva.powspctrm;
    stat.avg = dva.avg;
    
    megdat.stat{idva, ichan} = stat;
  end
end