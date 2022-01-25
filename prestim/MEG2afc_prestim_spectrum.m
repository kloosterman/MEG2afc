% respavg = MEG2afc_load_respavg()

% freqstatistics on TFR:
nsub=length(respavg.SUBJ);
% nsub=18;

loadstat = 0;
outfile = 'poolstat';
if loadstat;
    disp('Loading poolstat . . .')
    load(fullfile(PREIN, 'respavg', outfile));
    return
end

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

% cfg.latency = [-0.5 -0.2];
% cfg.avgovertime = 'yes';
% % cfg.channel = respavg.sens.ind{1};
% % pooling increased beta modukation in psc to resp
% cfg.channel = {'MLO11', 'MLO12', 'MLP31', 'MLP41', 'MLP42', 'MLP51', 'MLP52', 'MLP53', 'MLP54', 'MRO12', 'MRO13', 'MRO14', 'MRO24', 'MRP41', 'MRP53', 'MRP54', 'MZO01', 'MZP01'};
% cfg.avgoverchan = 'yes';

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

itrig = 3;

stat=[];
for ifreq = 3%:4 %3:4
    for ilatr = 1%:4 % latr_leg = {'modulation', 'latr wrt resp', 'latr wrt choice', 'latr wrt stim'};
        for idiff = 3 % 3:4 %1:4
            for idrug = 4 %3:4 %1:4 %1:4% 1:4
                
                %5D: lowfreq and highfreq apart
                freq=[];
                freq.freq = respavg.freq{ifreq};
                freq.time = respavg.time{itrig};
                freq.dimord = 'subj_chan_freq_time';
                freq.label = respavg.label;
                freq.powspctrm = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
                
                %                 % hi and low combined
%                 freq=[];
%                 freq.freq = [respavg.freq{3} respavg.freq{4}];
%                 freq.time = respavg.time{itrig};
%                 freq.dimord = 'subj_chan_freq_time';
%                 freq.label = respavg.label;
%                 lowfreqpow = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{3}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, 3, ilatr))); %tind{itrig}
%                 highfreqpow = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{4}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, 4, ilatr))); %tind{itrig}
%                 freq.powspctrm = cat(3, lowfreqpow, highfreqpow);

                freqzero = freq; %create zero freq to test against
                freqzero.powspctrm = zeros(size(freq.powspctrm));
                stat{ifreq} = ft_freqstatistics(cfg, freq, freqzero);
            end
        end
    end
end
% save(fullfile(PREIN, outfile), 'poolstat');

%% multiplot for above cond
cfg = [];
% cfg.xlim = [-0.8 0.3];
cfg.shading = 'flat';
cfg.layout = 'CTF275.lay';       %CTF275_RH.lay neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat

%             cfg.zlim = [-0.1 0.1];
cfg.zlim = 'maxabs';
% cfg.zlim = [-1 1];
load('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/critEEG_analysis/plotting/colormap_jetlightgray.mat')
cfg.colormap = cmap;

% freq.mask = stat{ ifreq}.mask;
freq.mask = stat{ ifreq}.negclusterslabelmat == 3;
% freq.mask = stat{ ifreq}.posclusterslabelmat == 1;

cfg.maskparameter = 'mask';

cfg.hotkeys = 'yes';
cfg.fontsize = 18;
cfg.colorbar = 'yes';

f = figure;
f.Position = [ 680   678   1200   1200];

ft_multiplotTFR(cfg, freq);


%% plot spectrum occ pooling

ifreq = 3;
itrig = 3;
idrug = 4;
idiff=3;
ilatr=1;

% poolingname = 'lateral alpha increase';
% chans = {'MLO41', 'MLO42', 'MLO43', 'MLO51', 'MLO52', 'MLO53', 'MLT11', 'MLT22', 'MLT23', 'MLT24', 'MLT32', 'MLT33', 'MLT34', 'MLT35', 'MLT36', 'MLT42', 'MLT43', 'MLT44', 'MLT45', 'MLT54', 'MRO41', 'MRO43', 'MRO51', 'MRO52', 'MRO53', 'MRT33', 'MRT34', 'MRT42', 'MRT43', 'MRT44', 'MRT45', 'MRT53', 'MZO03'};

poolingname = 'midocc alpha decrease';
chans = {'MLO21', 'MLO31', 'MLP41', 'MLP53', 'MLP54', 'MRO12', 'MRO13', 'MRP21', 'MRP31', 'MRP32', 'MRP33', 'MRP41', 'MRP42', 'MRP51', 'MRP52', 'MRP53', 'MZO02'};

chaninds = match_str(respavg.label, chans);

cols = {'r', 'b', 'g', 'k'}
frind = respavg.freq{3} > 8 & respavg.freq{3} < 12;

powdat = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
% powdat = squeeze(mean(powdat(:,respavg.sens.ind{1},:),2))
powdat = squeeze(mean(powdat(:,chaninds,:),2))

figure;
h(idrug)=shadedErrorBar(respavg.freq{ifreq}, mean(powdat), std(powdat)/sqrt(length(respavg.SUBJ)), cols{idrug}, 1 );
hold on
ssvep = 8.55;
% ssvep = 10;
% plot([ssvep ssvep], get(gca, 'YLim'), '--k' )
% plot([ssvep/2 ssvep/2], get(gca, 'YLim'), '--k' )
% plot([ssvep*1.5 ssvep*1.5], get(gca, 'YLim'), '--k' )
% plot([ssvep*2 ssvep*2], get(gca, 'YLim'), '--k' )
grid on
r = refline(0,0);
r.Color = 'k';
title(poolingname)


%% plot single subj spectra: outliers?
stat=[];
itrig = 3;
ifreq = 3; %:4 %3:4
ilatr = 1; %:4 % latr_leg = {'modulation', 'latr wrt resp', 'latr wrt choice', 'latr wrt stim'};
idiff = 3; % 3:4 %1:4
idrug = 4%:2; %3:4 %1:4 %1:4% 1:4

freq=[];
freq.freq = respavg.freq{ifreq};
%                 freq.time = respavg.time{itrig};
freq.dimord = 'subj_freq'; % subj_chan_freq_time
freq.label = {'custompooling'};
freq.powspctrm = double(mean(squeeze(respavg.pow(:,respavg.sens.ind{1}, 1:length(respavg.freq{ifreq}), ...
    1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr)),2)); % avg over chan
tind = respavg.time{itrig} >= -0.5 & respavg.time{itrig} <= -0.2;
freq.powspctrm = squeeze(mean(freq.powspctrm(:,:,:,tind),4));
freq.powspctrm = freq.powspctrm * 1e21;

% cmap = cbrewer( 'qual', 'Dark2', length(respavg.SUBJ));
% cmap = cbrewer( 'div', 'Spectral', length(respavg.SUBJ));
cmap = jet(nsub);

close all
f = figure; hold on
f.Position =[        1138          53         767        1052];

set(0, 'defaultaxesfontsize', 14)
ax=gca;
% ax.ColorOrder = summer(length(respavg.SUBJ));
ax.ColorOrder = cmap;
plot(freq.freq, freq.powspctrm, 'Linewidth',4);
legend(respavg.SUBJ)

refline(0,0)
legend(respavg.SUBJ)
figure; hold on
pl = plot(freq.freq, mean(freq.powspctrm))
refline(0,0)
% save(fullfile(PREIN, outfile), 'poolstat');

%% spectra drug and plac
figure; hold on
for idrug = 1:2
    freq.powspctrm = double(mean(squeeze(respavg.pow(:,respavg.sens.ind{1}, 1:length(respavg.freq{ifreq}), ...
        1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr)),2)); % avg over chan
    tind = respavg.time{itrig} >= -0.5 & respavg.time{itrig} <= -0.2;
    dum = squeeze(mean(freq.powspctrm(:,:,:,tind),4));
    plot(freq.freq, mean(dum))
end
legend(respavg.pharm_conds(1:2))


%% integrate chans, plot spectrum
% mask = stat{ifreq}.mask;
            load colormap_jetlightgray.mat

close all
sign = {'pos', 'neg'};
SAV=0;
isoi = 1;

for ifreq = 3%:4 %3:4
    
    for isign = 1:2
        fieldname = [sign{isign} 'clusters'];
        % mask = stat{ifreq}.posclusterslabelmat == 2;
%         nclus = find([stat{ifreq}.(fieldname).prob] < 0.025);
        nclus = find([stat{ifreq}.(fieldname).prob] <= 1);
%         nclus = find([stat{ifreq}.(fieldname).prob] < 0.9);
        if isempty(nclus); continue; end
        
%         nclus = length(unique(stat{ifreq}.(fieldname)))-1;
        
        fieldname = [sign{isign} 'clusterslabelmat'];

        for iclus = nclus
            idrug = 2;

            mask = stat{ifreq}.(fieldname) == iclus;
            
            f = figure; %imagesc(powdat)
%             f.Position =  [ 680   310   569   788];
%             f.Position =  [  2327         304         777         700];
            
                        % plot drug and plac raw spectra
            subplot(3,2,1); hold on
            cols = {'r', 'b'};
            clear h
            for idrug = 1:2
                powdat = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
                %integrate chans that are part of the cluster
                powdat = squeeze(trapz(powdat(:,any(mask,2), :),2));
                
%                 powdat = squeeze(trapz(powdat(:,respavg.sens.ind{isoi}, :),2)); % to get the occipital alpha
                
                h(idrug)=shadedErrorBar(respavg.freq{ifreq}, mean(powdat), std(powdat)/sqrt(length(respavg.SUBJ)), cols{idrug}, 1 );
                box on
                ax=gca;
                ax.FontSize = 16;
                h(idrug).mainLine.LineWidth=2;
            end
            legend([h(:).mainLine], {'ATX', 'Placebo'}); legend boxoff
%             legend([h(:).mainLine], respavg.pharm_conds(1:2)); legend boxoff
%             title('Raw power spectra of highlighted channels')
            xlabel('Frequency (Hz)')
            ylabel('MEG power')
%             ax.YLim = ax.YLim/2;
            ax.XLim = [0 35];
            ax.Position = [0.1300 0.68 0.4 0.25]
            
            idrug = 4;
            % %plot heatmap chan by freq
            % powdat = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
            % subplot(2,2,1)
            % ft_plot_matrix(1:length(freq.label), freq.freq, powdat', 'highlight', mask', 'clim', [-0.3 0.3])
            % colorbar
            
            %integrate chans that are part of the cluster
            subplot(2,1,2); hold on
            powdat = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
            
                        powdat = squeeze(trapz(powdat(:,any(mask,2), :),2));
            
%             powdat = squeeze(mean(powdat(:,respavg.sens.ind{isoi}, :),2)); % take occ pooling
            
            ax=gca;
            h=shadedErrorBar(respavg.freq{ifreq}, mean(powdat), std(powdat)/sqrt(length(respavg.SUBJ)), [], 0 );
            rl = refline(0,0); rl.Color = 'k'; rl.LineStyle = '--';
            box on
            % plos sig bar
            barind = find(diff([0 any(mask,1)]));
            if mod(length(barind),2), barind = [barind length(any(mask,1))]; end
            ypos = diff(ax.YLim)*-0.1;
            pl = plot( respavg.freq{ifreq}(barind), [ypos ypos] );
            pl.LineWidth = 10;
            pl.Color = 'k';
            ax.FontSize = 16;
            xlabel('Frequency (Hz)')
%             ylabel(sprintf('ATX-induced power wrt placebo\n(psc, integr.)'))
%             ylabel(sprintf('ATX-induced power change \n(psc, integr.)'))
            ylabel(sprintf('Power change from placebo\n(integrated %%)'))
            title(sprintf('Cluster p = %1.3f', stat{ifreq}.([sign{isign} 'clusters'])(iclus).prob ))
%             legend(h.mainLine, {'(ATX-placebo) ./ placebo'}); legend boxoff
            
            % plot topo
            subplot(3,2,2)
            powdat = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
            topodat = squeeze(mean(powdat,1));
            topodat = trapz(topodat(:, any(mask,1)),2);
            frtop=[];
            frtop.dimord = 'chan';
            frtop.label = respavg.label;
            frtop.freq = 1;
            % frtop.time = respavg.time{itrig};
            frtop.powspctrm = topodat;
            cfg = [];
            cfg.layout = 'CTF275.lay';
            cfg.marker = 'on';
            cfg.markersize = 3;
            cfg.markersymbol = '.'; %'o'
            cfg.comment = 'no';
            cfg.shading = 'flat';
            cfg.style = 'both'; %both 'imsat' straight
            cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
            cfg.parameter = 'powspctrm';
            cfg.interactive = 'no';
            cfg.zlim = 'maxabs';
            cfg.highlight = 'on';
            cfg.highlightchannel = frtop.label(any(mask,2));
%             cfg.highlightchannel = respavg.sens.ind{isoi};
            cfg.highlightsymbol = 'x';
            cfg.highlightsize = 8;
            cfg.colormap = cmap;
            
            ft_topoplotTFR(cfg, frtop)
%             colorbar
%             load colormap_jetlightgray.mat
%             colormap(cmap)
            ax=gca;
            ax.Position = [0.5 0.65 0.5 0.3]
            
        end
    end
    if SAV
        outpath = fullfile(respavg.PREOUT, 'prestim_spectra');
        warning off; mkdir(outpath); warning on
        outfile = sprintf( '%s%sstats_spectra_%s', outpath, filesep, respavg.freq_leg{ifreq}); % motor_conds{imotor}, pupilconds{ipup}
        disp(outfile)
        export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
        export_fig(outfile, '-tiff', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
        cd(outpath)
    end;
end
%% plot interactive topo to look at occipital cluster
close all
cfg.interactive = 'yes';
cfg.highlight = 'off';
% cfg.highlightchannel = respavg.sens.ind{6};
cfg.colorbar = 'yes';
cfg.colormap = cmap;
cfg.comment = 'auto';
cfg.zlim = 'maxabs';
            cfg.layout = 'CTF275.lay';
cfg.xlim = [-0.5 -0.2];
            
            
frtop.dimord = 'chan_freq_time';
idrug=4;
% powdat = double(squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %tind{itrig}
powdat = squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr)); %tind{itrig}

frtop.powspctrm = squeeze(mean(powdat));
frtop.freq = respavg.freq{ifreq};
frtop.time = respavg.time{itrig};
frtop.label = respavg.label;
figure;ft_topoplotTFR(cfg, frtop)

%%
freq.powspctrm = squeeze(mean(freq.powspctrm));
freq.dimord = 'chan_freq_time';

%% correlate sig frontal cluster with drift
pow = reshape(freq.powspctrm, 19, []);
pow = pow(2:19,:); % drop NK1, no drifts
pow = mean(pow(:, stat{ifreq}.mask(:)),2);
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

% % istim=1; iresp=1;
% % for icond = 4% 1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
% %     for iband = 1%:2%:2
% %         for itrig = 1 %1:2
% freq.mask = stat{ifreq}.mask;
% %             freq.mask = stat{ifreq}.posclusterslabelmat == 1;
% %             freq.mask = poolstat{isoi, iband, itrig, icond, istim, iresp}.negclusterslabelmat == 2;
% %             if ~any(freq.mask)
% %                 continue;
% %             end
% %             freq.time = poolstat{iband, icond, istim, iresp}.time;
% % % 
% %             freq=[];
% %             freq.label = respavg.label;
% %             freq.dimord = 'chan_freq_time';
% %             freq.time = respavg.time{itrig};
% %             freq.freq = respavg.freq{iband};
% % %             freq.powspctrm  = squeeze(respavg.pow(:,:, 1:length(respavg.freq{iband}), 1:23, ...
% % %                 iband, itrig, 4,icond, istim, iresp));
% %             freq.powspctrm  = squeeze(mean(respavg.pow(:,:, 1:length(respavg.freq{iband}), :, ...
% %                 iband, itrig, 4,icond, istim, iresp),1));
% % 

cfg = [];
% cfg.xlim = [-0.8 0.3];
cfg.shading = 'flat';
cfg.layout = 'CTF275.lay';       %CTF275_RH.lay neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat

%             cfg.zlim = [-0.1 0.1];
cfg.zlim = 'maxabs';
% cfg.zlim = [-1 1];
load('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/critEEG_analysis/plotting/colormap_jetlightgray.mat')
cfg.colormap = cmap;

% cfg.maskparameter = 'mask';

cfg.hotkeys = 'yes';
cfg.fontsize = 18;
cfg.colorbar = 'yes';

f = figure;
f.Position = [ 680   678   1200   1200];

ft_multiplotTFR(cfg, freq);

%             title( sprintf('%s, %s, %s, %s',  respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}, respavg.trigger_leg{itrig} ))
%         end
%     end
% end



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


