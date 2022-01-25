% respavg = MEG2afc_load_respavg()

%% freqstatistics on TFR:
nsub=length(respavg.SUBJ);

loadstat = 0;
outfile = 'poolstat';
if loadstat;
    disp('Loading poolstat . . .')
    load(fullfile(PREIN, 'respavg', outfile));
    return
end
freq=[];
freq.freq = respavg.freq{3};
freq.time = respavg.time{3};
freq.dimord = 'subj_chan_freq_time';
% freq.label = {'custompooling'};
freq.label = respavg.label;

cfg = [];
% cfg.frequency        = YLIM(3,:);
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
%     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% prepare_neighbours determines what sensors may form clusters
cfg0_neighb.method    = 'template'; %TODO do based on data (chans missing)
if ismac
    cfg0_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
else
    cfg0_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
end
cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);
% ft_neighbourplot(cfg0_neighb, freq)

% cfg.neighbours       = []; %in case no channel data present

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

poolstat=[];
for ilatr = 1%:4 % latr_leg = {'modulation', 'latr wrt resp', 'latr wrt choice', 'latr wrt stim'};
    for isoi = 1 %1%:2 %:2 % occ and motor
        for idiff = 3% 3:4 %1:4
            for idrug = 4 %3:4 %1:4 %1:4% 1:4
                for itrig =  3%:2 %1:4 %1:4% 1:4
                    %                     cfg.latency          = XLIM(itrig,:);
                    %                     freq.time = taxis{itrig};
                    
                    %average over SENSORS
                    %                     freq.powspctrm = squeeze(nanmean(respavg.dat(:,sens.ind{isoi}, frind{3}, 1:length(tind{itrig}) ,idrug, idiff, itrig, 3, ilatr), 2)); %tind{itrig}
                    %                     freq.powspctrm = shiftdim(freq.powspctrm, -1);
                    %5D: lowfreq
                    freq.powspctrm = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{3}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, 3, ilatr))); %tind{itrig}
                    
                    %                 freq.powspctrm = squeeze(respavg(:,isoi,:,tind{itrig}, idrug,3,idiff, itrig)); %DIMS:
                    
                    freqzero = freq; %create zero freq to test against
                    freqzero.powspctrm = zeros(size(freq.powspctrm));
                    %                     poolstat{ilatr, isoi, idrug, idiff, itrig} = ft_freqstatistics(cfg, freq, freqzero);
                    stat = ft_freqstatistics(cfg, freq, freqzero);
                end
            end
        end
    end
end
% save(fullfile(PREIN, outfile), 'poolstat');

%%
freq.powspctrm = squeeze(mean(freq.powspctrm));
freq.dimord = 'chan_freq_time';

%% correlate sig frontal cluster with drift
pow = reshape(freq.powspctrm, 19, []);
pow = pow(2:19,:); % drop NK1, no drifts
pow = mean(pow(:, stat.mask(:)),2);
[r,p] = corr(pow, respavg.driftrates(:,4,3), 'type', 'Pearson')

close all
figure; scatter(pow, respavg.driftrates(:,4,3))
lsline; box on; axis square
t = title(sprintf('(Drug-plac)./plac vs drift corr;\nr = %g, p = %g', r,p));
t.FontSize = 20;
xlabel('Pow (Drug-plac)./plac (psc)')
ylabel('drift rate Drug-plac regressed')
%% plot multiplot for significant clusters

close all

cfg = [];

% istim=1; iresp=1;
% for icond = 4% 1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
%     for iband = 1%:2%:2
%         for itrig = 1 %1:2

freq.powspctrm = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{3}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, 3, ilatr))); %tind{itrig}

% freq.mask = stat.mask;
freq.mask = stat.posclusterslabelmat == 1;
% freq.mask = stat.prob < 0.05;
% freq.powspctrm(:, ~freq.mask) = 0;

% cfg.maskparameter = 'mask';
% cfg.nanmean = 'yes';
%             freq.mask = stat.posclusterslabelmat == 1;
%             freq.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.negclusterslabelmat == 2;
%             if ~any(freq.mask)
%                 continue;
%             end
%             freq.time = poolstat{iband, icond, istim, iresp}.time;

%             freq=[];
%             freq.label = respavg.label;
%             freq.dimord = 'chan_freq_time';
%             freq.time = respavg.time{itrig};
%             freq.freq = respavg.freq{iband};
% %             freq.powspctrm  = squeeze(respavg.pow(:,:, 1:length(respavg.freq{iband}), 1:23, ...
% %                 iband, itrig, 4,icond, istim, iresp));
%             freq.powspctrm  = squeeze(mean(respavg.pow(:,:, 1:length(respavg.freq{iband}), :, ...
%                 iband, itrig, 4,icond, istim, iresp),1));


% cfg.xlim = [-0.8 0.3];
cfg.shading = 'flat';
% cfg.layout = 'CTF275.lay';       %CTF275_RH.lay neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat
cfg.layout = 'CTF275_helmet.lay';       %CTF275_RH.lay neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat

%             cfg.zlim = [-0.1 0.1];
cfg.zlim = 'maxabs';
% cfg.zlim = [-1 1];
load('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/critEEG_analysis/plotting/colormap_jetlightgray.mat')
cfg.colormap = cmap;

cfg.hotkeys = 'yes';
cfg.fontsize = 18;
cfg.colorbar = 'yes';

f = figure;
f.Position = [ 680   678   1200   1200];

ft_multiplotTFR(cfg, freq);

%% correlate clusters with ddm paras etc

powspctrm = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{3}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, 3, ilatr))); %tind{itrig}
powspctrm = powspctrm(:,:);

mask = stat.posclusterslabelmat == 1;
mask = mask(:);
megdat = mean(powspctrm(:,mask),2);

% vdat = respavg.ddmdat.v_hard_atx - respavg.ddmdat.v_hard_plac;
vdat = respavg.ddmdat.v_easy_atx - respavg.ddmdat.v_easy_plac;

r = corr(megdat, vdat, 'type', 'Spearman');
% close all
figure; scatter(megdat, vdat);
axis square; box on
text(megdat, vdat, respavg.SUBJ);

xlabel('alphabeta Poscluster atx?plac')
ylabel('ddm atx?plac')
title(r)

% tbl = table(megdat,vdat,'VariableNames',{'MEG','vdat'});
% lm = fitlm(tbl,'MEG~vdat', 'quadratic')

fitlm(megdat*1e27, vdat, 'Quadratic')

%% plot
load('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/critEEG_analysis/plotting/colormap_jetlightgray.mat')
close all
f = figure;
f.Position = [ 2004         415        1281         533];

% dimord: 'subj_chan_freq_time_drug_diff_trig_freqrange_latr'

iplot=0;
for idrug = [1,2,4]
    pow = squeeze(mean(respavg.pow(:,respavg.sens.ind{1},: ,:, idrug,3,1,1,1),2));
    pow = squeeze(mean(pow));
    
    
    
    scale = [-max(pow(:)) max(pow(:))];
    iplot = iplot+1;
    subplot(2,3,iplot);
    im = imagesc(respavg.time{1}, respavg.freq{1}, pow, scale );
    ax=gca;
    ax.YDir = 'normal';
    colorbar
    colormap(cmap)
    title(respavg.pharm_conds{idrug})
    xlabel('Time from stim onset')
    %     plot([0 0], ax.YLim)
end

tind = respavg.time{1} > -0.5 & respavg.time{1} < 0;
pow = squeeze(mean(respavg.pow(:,respavg.sens.ind{1},: ,:, idrug,3,1,1,1),2));
pow = squeeze(mean(pow(:,:,tind),3));
iplot = iplot+1;
subplot(2,3,iplot)
% plot(respavg.freq{1}, pow)
shadedErrorBar(respavg.freq{1}, squeeze(mean(pow)), std(double(pow)) / sqrt(length(respavg.SUBJ))  )
r = refline(0,0);
r.Color = 'k';
xlabel('Frequency (Hz)')
ylabel('Power')
title(respavg.pharm_conds{idrug})


