% builds on plotTFR_2AFC_poolings. ie plots in same figure

% MEG2afc_load_respavg

% driftrates = MEG2afc_load_driftrates161116(SUBJ, PREOUT);

driftrates = MEG2afc_load_drifts_regression;


corrtype = 'Pearson'; % 'Spearman'
%% do stats on correlation TFR's and driftrate

loadstat = 0;
outfile = 'corrstat';
if loadstat
    disp('Loading corrstat . . .')
    load(fullfile(PREIN, outfile));
    return
end

freq=[];
freq.freq = faxis{iband};
freq.dimord = 'chan_subj_freq_time';
freq.label = {'custompooling'};

cfg = [];
cfg.frequency        = YLIM(iband,:);
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_correlationT';
cfg.type             = corrtype;
cfg.design   = squeeze(driftrates(:, idrug, idiff))';
cfg.uvar     = [];
cfg.ivar     = 1;

cfg.correctm         = 'cluster';
% % % %     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg.neighbours       = []; %in case no channel data present

corrstat=[];
for isoi = 1:2%:2 %:2 % occ and motor
    for idiff = 3%:4 % 1:4 %[1,2,4] %1:4 %1:4
        for idrug = 4 %1:4 %1:4 %1:4 %[1,2,4] %1:4
            for itrig = 1:2 %1:4 %1:4 %[1,2,4] %1:4
                freq.time = taxis{itrig};
                cfg.latency          = XLIM(itrig,:);
%                 freq.powspctrm = squeeze(mean(respavg(:,sens.ind{isoi},:,tind{itrig}, idrug,3,idiff,itrig),2)); %DIMS:
                
                freq.powspctrm = squeeze(mean(respavg(:,sens.ind{isoi}, 1:length(frind{iband}), :, idrug, idiff, itrig, iband, 1), 2));
                
                %             freq.powspctrm = squeeze(mean(respavg(:,sens.ind{isoi},:,:, idrug,3,idiff,3,3),2)); %DIMS:
                freq.powspctrm = shiftdim(freq.powspctrm, -1);
%                 freqbehav = freq; %create zero freq to test against
%                 for isub=1:nsub
%                     freqbehav.powspctrm(1,isub,:,:) = squeeze(driftrates(isub,idrug,idiff));
%                 end
%                 corrstat{isoi, idrug, idiff, itrig} = ft_freqstatistics(cfg, freq, freqbehav);
                corrstat{isoi, idrug, idiff, itrig} = ft_freqstatistics(cfg, freq);
            end
        end
    end
end
save(fullfile(PREIN, outfile), 'corrstat');


%% Correlate TFR's with driftrate and plot
addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/custom_tools/stats')
addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'))

close all
SAV = 1;
set(0,'DefaultAxesFontsize',12)

showstats = 1;

iband = 3; % faxis{3}

TYP = '2afc';
for isoi = 1:2   %:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
    %     figh = gcf;  
%     figure%(figh)
    figh = figure;
    figh.Position = [ 2310         499         800         600];
    
%     iplot=6; hold on;
    iplot=0; hold on;
%     set(figh, 'Position', [0 -200 375*4 210*3])
%     load('colormap170613.mat');
%     load colormap_jetlightgray.mat
    cmap = cbrewer('div', 'RdBu',256);
    cmap = flipud(cmap);

    
    colormap(cmap);  %cmap = get(gcf, 'Colormap')
    for idiff= 3 %1:3 %[1:2,4]
        for idrug = 4 %[1:2,4] % 1:2
            for itrig = 1:2
                
                ZLIM = [-1 1];
%                 iplot = iplot+1; subplot(4,6,iplot); cla; hold on
%                 iplot = iplot+1; subplot(4,6,iplot); cla; hold on
                iplot = iplot+1; subplot(2,2,iplot); cla; hold on

                dum = squeeze( nanmean(respavg(:,sens.ind{isoi}, 1:length(frind{iband}), :, idrug, idiff, itrig, iband, 1), 2));
                
%                 dum = squeeze( nanmean(respavg(:,sens.ind{isoi},frind,tind{itrig},idrug, imotor, idiff, itrig), 2 ) ); %average over subj
                %             dum = squeeze(respavglat(:,isoi,frind,:,idrug, imotor, idiff, 3, ipup)); %average over subj
%                 rcoeff=[];
%                 for ifreq = 1:size(dum,2)
%                     for itind = 1:size(dum,3)
%                         rcoeff(ifreq,itind) = corr(squeeze(driftrates(:,idrug,idiff)), squeeze(dum(:,ifreq,itind)), 'type', corrtype);
%                     end
%                 end
                
                rho = squeeze(corrstat{isoi, idrug, idiff, itrig}.rho);
                if showstats
                    mask = double(squeeze( corrstat{isoi, idrug, idiff, itrig}.mask));
                    mask(mask==0) = 0.25;
                    %                             mask = double(squeeze(poolstat{isoi, iband, itrig, icond, istim, iresp}.prob < 0.15));
                    %                             mask(mask==0) = 0.25;
                    
                    ft_plot_matrix(taxis{itrig}, faxis{iband}, rho, 'clim', ZLIM, 'box', 'no', ...
                        'highlight', mask, 'highlightstyle', 'opacity'); % opacity
                else
                    ft_plot_matrix(respavg.time{itrig}, respavg.freq{iband}, rho, 'clim', ZLIM, 'box', 'no' ); % opacity
                end
               
                title({sprintf('%s corr Mod-Driftr. %s %s ', corrtype, sens.leg{isoi}, pharm_conds{idrug}) ; 
                    sprintf('%s [%g]',  diff_conds{idiff}, ZLIM(2) )}); % 'FontWeight','bold'
                hold on
                yaxis = [faxis{iband}(1),faxis{iband}(end)];
                plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                
                xlim(XLIM(itrig,:))
                ylim([faxis{iband}(1),faxis{iband}(end)])
                YTICKS = [0:10:200];

                set(gca,'Box','off','XTick',-2:0.5:2,...    [-0.5,0,0.5,1,1.5,2,2.5]
                    'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                    'TickDir','out', 'FontSize', 12);
%                 if itrig == 1
                    ylabel('Frequency (Hz)');
%                 else
%                     set(gca, 'YTickLabel', []);
%                 end
                LAB=1;
                xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
                box on
                                h=colorbar;
%                 if isoi == 1
%                     fois = [19 25];
%                     tois = [0.75 1.1; -0.9 -0.25]; %stim then resp-locked
%                 elseif isoi == 2
%                     fois = [15 22];
%                     tois = [0 0.25; -0.9 -0.4]; %stim then resp-locked
%                 end
%                 rectangle('Position', [tois(itrig,1), fois(1), tois(itrig,2)-tois(itrig,1), fois(2)-fois(1)], 'Linewidth', 1.5)

            end
        end
    end
    
    if SAV
        outpath = fullfile(PREOUT, 'poolings');
        warning off; mkdir(outpath); warning on
        outfile = sprintf( '%s%sTFR_corr_driftrate_%slocked_%s_%sfreq_%s', outpath, filesep, trigger_leg{itrig}, TYP, analysistype{iband}, sens.leg{isoi} ); % motor_conds{imotor}, pupilconds{ipup}
        disp(outfile)
        export_fig(outfile, '-pdf', '-painters', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
%         export_fig(outfile, '-eps') %'-png',  '-pdf',
%         print(outfile, '-dpng', '-painters') %'-png',  '-pdf',
    end;
end
cd(outpath)






%% Correlate drift rate with significant TF clusters
     
addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/toolbox/stats');

% close all
SAV=1;
freq=[];
freq.time = taxis;
freq.freq = faxis{iband};
freq.dimord = 'subj_freq_time';
freq.label = {'custompooling'};

cfg=[];
cfg.trials      = 'all';                cfg.avgovertime = 'yes';          cfg.avgoverfreq = 'yes';

figh = figure;

% figure(figh)
% set(0,'DefaultAxesFontsize',8)

iplot=0;
% set(gcf, 'Position', [0 -100 375*3 210*4])
for isoi = 1 %:2 %:2 % occ and motor
    for imeas = 1 %:2 % modulation or lateralization
        for idiff = 3 % 1:4 %1:3 %[1,2,4] %1:4 %1:4
            for idrug = 4 %[1,2,4] % 1:4 %1:4
                for itrig = 2 %1:2
%                     freq.powspctrm = squeeze(mean(respavg(:,sens.ind{isoi},frind,tind{itrig}, idrug,3,idiff,itrig), 2)); %DIMS:
%                     freq.time = taxis{itrig};
%                     cfg.latency     = tois(itrig,:);
%                     cfg.frequency   = fois;
%                     freqavg = ft_selectdata(cfg, freq);
%                     [r, p] = corr(freqavg.powspctrm, squeeze(driftrates(:,idrug, idiff)), 'type', corrtype);

                    dum = squeeze( mean(respavg(:,sens.ind{isoi}, 1:length(frind{iband}), :, idrug, idiff, itrig, iband, 1), 2));
                    mask = squeeze( corrstat{isoi, idrug, idiff, itrig}.mask);
                    dum = mean(dum(:,mask),2);
                    
                    [r, p] = corr(dum, squeeze(driftrates(:,idrug, idiff)), 'type', corrtype);
                    %                 [p, r] = permuteCorr(freqavg.powspctrm, squeeze(driftrates(:,idrug, idiff)));
                    
                    iplot=iplot+1; subplot(1,1,iplot); cla; hold on; axis square
                    for isub=1:nsub
%                         plot(freqavg.powspctrm(isub), driftrates(isub,idrug, idiff), 'o', 'MarkerSize', 8, ...
%                             'Color', 'white', 'markerfacecolor', [0 0 0], 'Linewidth', 1);
                        plot(dum(isub), driftrates(isub,idrug, idiff), 'o', 'MarkerSize', 16, ...
                            'Color', 'white', 'markerfacecolor', [0 0 0], 'Linewidth', 1);
                        
                        text(double(dum(isub)),  driftrates(isub,idrug, idiff), SUBJ{isub}, 'Fontsize', 12);
                        
                    end
                    %                 scatter(freqavg.powspctrm, squeeze(driftrates(:,idrug, idiff)))
                    if p < 0.05 && strcmp(corrtype, 'Pearson')
                        %                     lsline
                        beta=polyfit(dum, squeeze(driftrates(:,idrug, idiff)),1);
                        h_ls = refline(beta);
                        set(h_ls, 'Linewidth', 1.5, 'Color', [0 0 0])
                    end
%                     title({sprintf('%s, r = %0.2f, p = %0.2f', pharm_conds{idrug}, r,p), ...
%                         sprintf('%s %s %d-%d Hz, %g-%g s', diff_conds{idiff}, sens.leg{isoi}, cfg.frequency, cfg.latency)})
                    title({sprintf('%s, r = %0.2f, p = %0.5f', pharm_conds{idrug}, r,p), ...
                        sprintf('%s %s ', diff_conds{idiff}, sens.leg{isoi})})
                    set(gca, 'FontSize', 12);
                    
                    hline = refline(0,0);
                    hline.Color = 'k';
                    hline.LineStyle = '--';
                    
                    
                    if imeas==1; xlabel('Modulation (%)'), else xlabel('Latr (%)'), end
                    ylabel('Drift rate')
                    box on
                    YLIM = get(gca, 'YLIM');
                    vline = plot([0 0], YLIM);
                    ylim(YLIM)
                    vline.Color = 'k';
                    vline.LineStyle = '--';
                    
                end
            end
        end
    end
end
if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%sscatter_corr_driftrate__%s_%sfreq_%s', outpath, filesep, TYP, analysistype{iband}, sens.leg{isoi} ); % motor_conds{imotor}, pupilconds{ipup}
    disp(outfile)
    export_fig(outfile, '-pdf', '-painters', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
%             export_fig(outfile, '-eps') %'-png',  '-pdf',
%     print(outfile, '-dpsc', '-painters') %'-png',  '-pdf',
end;



