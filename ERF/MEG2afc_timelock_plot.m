% function MEG2afc_timelock_plot(respavg)

load colormap_jetlightgray.mat

drugleg = {'drug' 'plac'};

close all
cfg=[];
cfg.layout = 'CTF275.lay';
cfg.colorbar = 'yes';
cfg.xlim = [-0.5 1.5];
cfg.zlim = 'maxabs';  % zeromax
cfg.nanmean = 'yes';
cfg.colormap = cmap;

% layout = ft_prepare_layout(cfg)
% layout.pos(1:132,:) = layout.pos(1:132,:)

figure

iltr = 2;
itrig = 1;
idrug = 3;
idiff = 3;
if iltr > 1
  cfg.outline = 'convex';
  cfg.mask        = 'convex';
end


% cfg2 = [];
% cfg2.latency = [-0.5 1.5];
% respavg{iltr, itrig, 1, idiff} = ft_selectdata(cfg2, respavg{iltr, itrig, 1, idiff});
% respavg{iltr, itrig, 2, idiff} = ft_selectdata(cfg2, respavg{iltr, itrig, 2, idiff});

ft_multiplotER(cfg, respavg{iltr, itrig, idrug, idiff})
legend(drugleg);

% ft_multiplotER(cfg, respavg{1, 4, 3})
% legend({'drug - plac'});

% ft_multiplotER(cfg, respavg{1, 2, 1}, respavg{1, 2, 2})
% legend({'easy' 'hard'})

%% 
figure
cfg=[];
cfg.layout = 'CTF275.lay';
% cfg.outline = 'convex';
% cfg.mask        = 'convex'; 
cfg.parameter = 'individual';
ft_topoplotER(cfg, respavg{iltr, itrig, 1, idiff})


%%
cfg=[];
ft_databrowser(cfg, respavg{2, 3, 3})
%% timelockstatistics 
nsub = size(respavg{1}.individual,1);

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
% cfg.alpha            = 0.1;
cfg.numrandomization = 100;
cfg.neighbours       = []; %in case no channel data present
cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = 'ctf275_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);

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

iltr = 2;
itrig = 1;
idrug = 4;
idiff = 3;
timelock = respavg{iltr, itrig, idrug, idiff};

timelockzero = timelock; %create zero freq to test against
timelockzero.individual = zeros(size(timelock.individual));

% if itrig == 1
%   cfg2 = [];
%   cfg2.latency = [-0.5 1];
%   timelock = ft_selectdata(cfg2, timelock)
% end
poolstat = ft_timelockstatistics(cfg, timelock, timelockzero);

poolstat.posclusters(1:5).prob
poolstat.negclusters(1:5).prob
%% multiplot

close all
cfg=[];
cfg.layout = 'CTF275.lay';
cfg.colorbar = 'yes';
cfg.zlim = 'zeromax';
% cfg.xlim = [-0.2 1];
drug = respavg{iltr, itrig, 4, idiff};
plac = respavg{iltr, itrig, 4, idiff};
% drug.mask = poolstat.mask;
% plac.mask = poolstat.mask;
drug.mask = poolstat.negclusterslabelmat == 1;
plac.mask = poolstat.negclusterslabelmat == 1;
cfg.maskparameter = 'mask';
figure
ft_multiplotER(cfg, drug, plac)

%% ca. 10 hz difference drug and plac!
cfg=[];
cfg.layout = 'CTF275.lay';
cfg.colorbar = 'yes';
cfg.zlim = 'zeromax';
cfg.xlim = [-0.2 1];
ft_multiplotER(cfg, respavg{1, 1, 1:2, 3});

% legend(drugleg);


%% clusterplot
cfg = [];
% cfg.highlightsymbolseries = ['*','*','.','.','.'];
cfg.layout = 'CTF275.lay';
cfg.contournum = 0;
cfg.visible = 'off';
cfg.saveaspng = 'test';
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.parameter='stat';
cfg.zlim = [-5 5];
ft_clusterplot(cfg, poolstat);

%%
load colormap_jetlightgray.mat
% drug = ft_timelockanalysis([], respavg{1, 4, 3}) ;

iltr = 2;
itrig = 1;
idrug = 4;
idiff = 3;
cfg = [];
if iltr > 1
  cfg.outline = 'convex';
  cfg.mask        = 'convex';
end

drug = ft_timelockanalysis([], respavg{iltr, itrig, idrug, idiff}) ;
cfg.colorbar = 'yes';
cfg.parameter = 'avg';
% cfg.xlim = [-0 1];
cfg.xlim = [-0.5 0.25];
cfg.zlim = 'maxabs';
% cfg.layout = 'CTF275_helmet.mat';
cfg.layout = 'CTF275.lay';
f = figure; 
f.Position = [680 284 560 814];
colormap(cmap);
ft_movieplotER(cfg, drug);