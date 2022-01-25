%% plot TFR's for selected poolings
% function plotTFR_2AFC(~)
% plots SOI-TFRs
% average of selected sensors for multiple subjects
% sorted according to stimulus condition
% % % % % addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
% % % % % addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
% % % % % ft_defaults
cd('/Users/kloosterman/gridmaster2012/kloosterman')
% MEG2afc_setup_paths

% MEG2afc_setup_paths
% trigger = 'stim' % both are now loaded!!

MEG2afc_load_respavg

%% freqstatistics on TFR:
% test motor and occ vs 0
% motor and occ drug vs placebo
nsub=length(SUBJ);

loadstat = 0;
outfile = ['poolstat_' trigger];
if loadstat;
    disp('Loading poolstat . . .')
    load(fullfile(PREIN, 'respavg', outfile));
    return
end
freq=[];
freq.time = taxis;
freq.freq = faxis;
freq.dimord = 'chan_subj_freq_time';
freq.label = {'custompooling'};

cfg = [];
cfg.latency          = XLIM;
cfg.frequency        = YLIM;
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% % % %     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg.neighbours       = []; %in case no channel data present

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
for isoi = 1:2 %:2 % occ and motor
    for idiff = 3 %1:4 %1:4
        for idrug = 1:4 %1:4 %1:4% 1:4
            freq.powspctrm = squeeze(respavg(:,isoi,:,:, idrug,3,idiff,3,3)); %DIMS:
            freq.powspctrm = shiftdim(freq.powspctrm, -1);
            freqzero = freq; %create zero freq to test against
            freqzero.powspctrm = zeros(size(freq.powspctrm));
            poolstat{isoi, idrug, idiff} = ft_freqstatistics(cfg, freq, freqzero);
        end
    end
end
save(fullfile(PREIN, 'respavg', outfile), 'poolstat');

%% plotting: 2AFC TFR drug vs placebo: motor and visual cortex + stats
close all
SAV=1;
set(0,'DefaultAxesFontsize',10)

showstats = 0;
imotor=1;
istim=3;
iresp=3;

TYP = '2afc';
for isoi = 1:2 %:4 %1:3%  1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
        figure;    iplot=0; hold on
        set(gcf, 'Position', [0 -200 375*4 210*3])
        load('colormap170613.mat');
        colormap(cmap);  %cmap = get(gcf, 'Colormap')
        % colormap(jet(256));  %cmap = get(gcf, 'Colormap')
        %         for ipup= [3 1 2 4] %1:4
        for idiff = 1:4
            for idrug = [1:2,4] % 1:2 1:4 %1:4 % 
                for itrig = 1:2
                    %         ZLIM = [-0.05 0.05];
                    ZLIM = [-0.4 0.4];
                    iplot = iplot+1; subplot(4,6,iplot); hold on
                    if iplot > 12;  ZLIM = [-0.1 0.1];  end
                    if idiff == 4;  ZLIM = [-0.05 0.05];  end
                    if idrug == 4;  ZLIM = [-0.05 0.05];  end
                    
%                     dum = squeeze(mean(respavg(:,sens.ind{isoi},frind, tind{itrig} ,idrug, imotor, idiff, istim, iresp, itrig))); %average over subj
                    dum = squeeze(nanmean(respavg(:,sens.ind{isoi},frind, tind{itrig} ,idrug, imotor, idiff, itrig))); %average over subj
                    dum = squeeze(mean(dum)); % average over sensors
                    if showstats
                        dum = showstatsTFR(dum, poolstat{isoi, idrug, idiff}.prob, frind, ZLIM, showstats);
                    end
                    imagesc(taxis{itrig}, faxis, dum,ZLIM);
                    
                    title({sprintf('%s %s %slocked [%g]', motor_conds{imotor}, diff_conds{idiff}, trigger_leg{itrig}, ZLIM(2)) ;
                        sprintf('%s  %s ', sens.leg{isoi}, pharm_conds{idrug} )}); % 'FontWeight','bold'
                    hold on
                    yaxis = [faxis(1), faxis(end)];
                    plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                    
                    xlim(XLIM(itrig,:))
                    ylim([FREQLO 35])
                    YTICKS = [0:5:200];
                    set(gca,'Box','off','XTick',-1:0.5:2,...    [-0.5,0,0.5,1,1.5,2,2.5]
                        'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                        'TickDir','out', 'FontSize', 8);
                    if itrig == 1
                        ylabel('Frequency (Hz)');
                    else
                        set(gca, 'YTickLabel', []);
                    end                    
                    xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
                    %                                         h=colorbar;
                end
            end 
        end
    
        if SAV
            outpath = fullfile(PREOUT, 'poolings');            
            warning off; mkdir(outpath); warning on
            outfile = sprintf( '%s%sTFRpooling_%s_%sfreq_%s_%s', outpath, filesep, TYP, analysistype, sens.leg{isoi}, diff_conds{idiff} ); % motor_conds{imotor}, pupilconds{ipup}
            disp(outfile)
            export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
        end;
end
cd(outpath)



%% Plot bars for the frontal gamma: 49-61 Hz, -0.4 0.2 s

SAV=1;
% close all
isoi=2;
frindoi = find(faxis <= 61 & faxis >= 49);
% taxis(tind) : -0.4 to 0.2
idrug = 1:2;
idiff = 1:2;
dum = squeeze(respavgstat(:,isoi,frindoi,tind,idrug, 3, idiff, 3, 3)); %average over subj
dum=mean(dum,2);
dum=squeeze(mean(dum,3)); % subj drug diff remain

sem = squeeze(std(dum) / sqrt(18));
pvaleasy = randtest1d(squeeze(dum(:,1,1)), squeeze(dum(:,2,1)), 0, 1000);
pvalhard = randtest1d(squeeze(dum(:,1,2)), squeeze(dum(:,2,2)), 0, 1000);

figure
barweb(squeeze(mean(dum))', sem', [], {'Easy' 'Hard'})
% legend(diff_conds(idiff))
legend(pharm_conds(idrug))
legend BOXOFF
title(sprintf('Easy: p = %g, Hard: p = %g', pvaleasy, pvalhard))
shading flat
ylabel('Modulation (%)')

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'bargraphs');
    
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%sGammadrugeffect_EasyVShard_%slocked_%s_%s_%sfreq_%s', outpath, filesep, trigger, TYP, sens.leg{isoi}, analysistype, motor_conds{imotor}  );
    disp(outfile)
    %             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
    %         print('-dpng', outfile) %, '-cmyk'
    export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
end;

%% plot spectrum of stimlocked activity: 0.25 to 0.5 s

cfg = [];
cfg.frequency        = 'all';
% cfg.latency          = [TIMLO TIMHI];
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
cfg.neighbours       = []; %in case no channel data present

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

freq2stats = freq;
freq2stats.dimord = 'subj_freq';
freq2stats.label = {'custompooling'};
freq2stats.time = 1;

SAV=1;
close all
if strcmp(trigger, 'stim')
    toi = [0.25 0.5];
else
    toi = [-0.5 0];
end
toiind = find(freq.time <= toi(2) & freq.time >= toi(1) );

figure; iplot=0;
set(gcf, 'Position', [0 0 375*2 210*4])
for isoi = 1:2
    iplot=iplot+1; subplot(2,1,iplot)
    dum = squeeze(mean(respavgstat(:,isoi,frind,toiind, 3, 3, 3, 3, 3), 4)); %average over tind
    
    freq2stats.powspctrm = dum; %DIMS:
    freq2statszero = freq2stats; %create zero freq to test against
    freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
    stat = ft_freqstatistics(cfg, freq2stats, freq2statszero);

    avg = squeeze(mean(dum));
    sem = std(dum) / sqrt(18);
    
    [~, minind]=min(avg);
    [~, maxind]=max(avg(5:end));
    
    h = shadedErrorBar(faxis,avg,sem);
    hold on
    plot([0 150], [0 0],'k')
    set(gca, 'Tickdir', 'out', 'XTick', 0:10:150)
    set(h.mainLine, 'Linewidth', 2)
    grid ON
    xlabel('Frequency (Hz)')
    ylabel('Modulation (%)')
    title(sprintf('Power spectrum %s across all conditions; interval %g-%g s %s-locked; min @ %d Hz, max @ %d Hz', sens.leg{isoi}, toi, trigger, faxis(minind), faxis(maxind) ))
    ylim([-0.5 0.2])
    
    pval = stat.prob';  % stat from ft_freqstatistics
    BARPOS= [-0.4];
    SIGCOL= ['k'];
    THR=0.025;
    %                 THR=1;
    BARWID= 5;
    sigint={};
    i = 1; sigint{i} = [];
    for ismp = 1:length(pval)
        if pval(ismp) < THR
            sigint{i} = [sigint{i} ismp];
        end % concatenate significant samples in sigint{iint}
        if (ismp < size(pval,1)) && (pval(ismp+1) >= THR) && ~isempty(sigint{i}) % jump to next interval if next sample is not sig
            i= i + 1;
            sigint{i} = [];
        end
    end
    % replot significant intervals in different color
    for iint = 1:length(sigint)
        if ~isempty(sigint{iint})
            begsmp = sigint{iint}(1); %*2
            endsmp = sigint{iint}(end); %*2
            
%             axylim = get(AX(1), 'YLim');
%             range = axylim(2) - axylim(1);
%             ylimsteps = axylim(1):range/20:axylim(2);
            
            %             plot( AX(1), stat{isoi, iband, idrug, idiff}.time(begsmp:endsmp),...   %begsmp-1:endsmp+1
            %                 ylimsteps(end-idrug)*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
            %                 'Color',SIGCOL,'LineWidth',BARWID)
            plot(stat.freq(begsmp:endsmp), ...   %begsmp-1:endsmp+1
                BARPOS*(ones(1,length(begsmp:endsmp))), ...  %begsmp-1:endsmp+1
                'Color',SIGCOL,'LineWidth',BARWID)
        end
    end

    
    
end
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    
    outfile = sprintf( '%s%sSpectrumAllconds_%slocked', outpath, filesep, trigger );
    disp(outfile)
    %             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
    %     print('-dpng', outfile) %, '-cmyk'
    %     orient landscape
    %     print('-depsc', outfile) %, '-cmyk'
    %     saveas(gcf, [outfile '.pdf'])
    export_fig(outfile, '-png', '-transparent') % '-png',
end;





