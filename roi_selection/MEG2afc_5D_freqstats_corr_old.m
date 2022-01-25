function MEG2afc_5D_freqstats_corr(cfg)

% Run 5D correlation stats on meg vs driftrate or modulation vs 0
% testdv = {'modulation'; 'correlation'}; % test power modulation or power-drift correlation

% read data
MEG2afc_load_respavg
driftrates = MEG2afc_load_driftrates161116(SUBJ, PREOUT);

% stack low and high freq vertically
temp=[];
temp = cat(3,  respavg(:,:, 1:length(frind{1}), :, :, :, :, 1), respavg(:,:, 1:length(frind{2}), :, :, :, :, 2)  );

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
cfg.numrandomization = 100;
cfg.avgoverfreq      = 'no';
%     cfg.subj             = 1:28;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'template'; %TODO do based on data (chans missing)
if ismac
    cfg_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20150803/template/neighbours/ctf275_neighb.mat';
else
    cfg_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20150803/template/neighbours/ctf275_neighb.mat';
end
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb);
% ft_neighbourplot(cfg_neighb, freq)


freq = [];
freq.label = chlabel;
freq.freq = faxis_all;
freq.time = taxis{cfg.itrig};   
freq.dimord = 'subj_chan_freq_time';
freq.powspctrm = squeeze(temp(:,:,:,:, cfg.idrug, cfg.idiff, cfg.itrig)); %DIMS: sub chan freq tind

freq2 = freq;
switch cfg.testdv
    case 'correlation'
        for isub=1:nsub
            freq2.powspctrm(isub,:,:,:) = squeeze(driftrates(isub,cfg.idrug,cfg.idiff));
        end
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

        freq2.powspctrm = zeros(size(freq.powspctrm));
end  

stat = ft_freqstatistics(cfg, freq, freq2);

outfile = sprintf('stat_idrug%d_idiff%d_itrig%d_%s_%s', cfg.idrug, cfg.idiff, cfg.itrig, cfg.testdv, cfg.clusterthreshold);
fprintf('SAVING TO stat_idrug%d_idiff%d_itrig%d_%s_%s\n', cfg.idrug, cfg.idiff, cfg.itrig, cfg.testdv, cfg.clusterthreshold);
save(fullfile(PREIN, outfile), 'stat');

