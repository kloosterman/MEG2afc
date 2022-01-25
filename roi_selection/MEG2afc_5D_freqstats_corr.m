function MEG2afc_5D_freqstats_corr(cfg, respavg)

% Run 5D correlation stats on meg vs driftrate or modulation vs 0
% testdv = {'modulation'; 'correlation'}; % test power modulation or power-drift correlation

% read data
% MEG2afc_load_respavg
driftrates = MEG2afc_load_driftrates161116;

freq = [];
freq.label = cfg.chlabel;
freq.freq = cfg.faxis_all;
freq.time = cfg.taxis{cfg.itrig};
freq.dimord = 'subj_chan_freq_time';
% freq.powspctrm = squeeze(temp(:,:,:,:, cfg.idrug, cfg.idiff, cfg.itrig)); %DIMS: sub chan freq tind
% stack low and high freq vertically
% freq.powspctrm = squeeze(cat(3,  respavg(:,:, 1:length(frind{1}), :, cfg.idrug, cfg.idiff, cfg.itrig, 1), ...
%     respavg(:,:, 1:length(frind{2}), :, cfg.idrug, cfg.idiff, cfg.itrig, 2)  ));
freq.powspctrm = squeeze(cat(3,  respavg(:,:, 1:length(cfg.frind{1}), :, 1), ...
    respavg(:,:, 1:length(cfg.frind{2}), :, 2)  ));

% temp = cat(3,  respavg(:,:, 1:length(frind{1}), :, :, :, :, 1), respavg(:,:, 1:length(frind{2}), :, :, :, :, 2)  );
clear respavg

% added to cfg:
cfg.channel          = {'MEG'};
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
switch cfg.testdv
    case 'correlation'
        cfg.statistic        = 'ft_statfun_correlationT';
        cfg.type             = 'Spearman';
%         cfg.type             = 'Pearson';
    case 'modulation'
        cfg.statistic        = 'ft_statfun_depsamplesT';
end

cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg.avgoverfreq      = 'no';
%     cfg.subj             = 1:28;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'template'; %TODO do based on data (chans missing)
if ismac
    cfg_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/ctf275_neighb.mat';
else
    cfg_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/ctf275_neighb.mat';
end
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb);
% ft_neighbourplot(cfg_neighb, freq)

switch cfg.testdv
    case 'correlation'
        cfg.design   = squeeze(driftrates(:,cfg.idrug,cfg.idiff))';
        cfg.uvar     = [];
        cfg.ivar     = 1;
        stat = ft_freqstatistics(cfg, freq);
        
    case 'modulation'
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
        
        freq2 = freq;
        freq2.powspctrm = zeros(size(freq.powspctrm));
        stat = ft_freqstatistics(cfg, freq, freq2);
end


outfile = sprintf('stat_idrug%d_idiff%d_itrig%d_%s', cfg.idrug, cfg.idiff, cfg.itrig, cfg.testdv);
fprintf('SAVING TO stat_idrug%d_idiff%d_itrig%d_%s\n', cfg.idrug, cfg.idiff, cfg.itrig, cfg.testdv);
save(fullfile(cfg.PREIN, outfile), 'stat');

