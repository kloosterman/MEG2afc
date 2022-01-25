load('/Users/kloosterman/beegfs/projectdata/MEG2afc/freqzap/NK3_plac_ipsi_freq.mat')
freqzap=freq;
load('/Users/kloosterman/beegfs/projectdata/MEG2afc/freq/NK3_plac_ipsi_freq.mat')
freqdft=freq;

%%
close all
cfg=[];
cfg.layout = 'CTF275_helmet.mat';
% cfg.baseline = [-0.25 0];
% cfg.baselinetype = 'relchange';
cfg.zlim = 'absmax';
figure; ft_singleplotTFR(cfg, freqzap(2,2,2,1));figure; ft_singleplotTFR(cfg, freqdft(2,2,2,1));
% figure; ft_multiplotTFR(cfg, freqzap(1,1,2,1));figure; ft_multiplotTFR(cfg, freqdft(1,1,2,1));
