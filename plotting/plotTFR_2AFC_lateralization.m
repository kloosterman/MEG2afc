%% Lateralization wrt button press for motor and occ poolings
% first run plotTFR_2AFC_poolings

clear all
trigger = 'resp'

MEG2afc_load_respavg

fprintf(' . . . DONE \n')

%% freqstatistics on lateralization:
% test motor and occ lateralization vs 0
% motor and occ lateralization drug vs placebo

freq=[];
freq.time = taxis;
freq.freq = faxis;
freq.dimord = 'subj_freq_time';
freq.label = {'custompooling'};

FREQLO = 5;
FREQHI = 100;

XLIM = [TIMLO TIMHI];
YLIM = [FREQLO FREQHI];

loadstat = 0;
outfile = ['latstat_' trigger];
if loadstat;
    disp('Loading latstat . . .')
    load(fullfile(PREIN, 'respavg', outfile));
    return
end

cfg = [];
cfg.latency          = XLIM;
cfg.frequency        = YLIM;
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster'; % no
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
freq2stats.dimord = 'subj_freq_time';
freq2stats.label = {'custompooling'};

% test motor and occ lateralization vs 0
latstat=[];
for icorrect = 1:3
    for isoi = 1:3 %1:2 % occ and motor
        for idiff = 3 % 1:3
            for idrug = 1:4
                freq2stats.powspctrm =  squeeze(respavglat(:,isoi,:,:,idrug, 3, idiff, 3, icorrect));
                freq2statszero = freq2stats; %create zero freq to test against
                freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
                latstat{isoi, idrug, idiff, icorrect} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
            end
            % test motor and occ lateralization drug vs placebo
            freq2stats.powspctrm =  squeeze(respavglat(:,isoi,:,:,1, 3, idiff, 3, icorrect));
            freq2stats2 = freq2stats;
            freq2stats2.powspctrm =  squeeze(respavglat(:,isoi,:,:,2, 3, idiff, 3, icorrect));
            
            latstat{isoi, 4, idiff, icorrect} = ft_freqstatistics(cfg, freq2stats, freq2stats2);
        end
    end
end
save(fullfile(PREIN, 'respavg', outfile), 'latstat');

%% plot lateralization TFR's with stats
close all
SAV=1;
set(0,'DefaultAxesFontsize',8)

showstats = 1;

ZLIM = [-0.015 0.015];

% tind = find(taxis <= XLIM(2) & taxis >= XLIM(1));
frind = find(faxis <= YLIM(2) & faxis >= YLIM(1));
tind = find(taxis <= XLIM(2) & taxis >= XLIM(1));

imotor=3; idiff=3; istimorresp = 3; % icorrect=3; %irt= 4 ; icor=3;

TYP = '2afc';
for icorrect = 1:3
    for idiff= 3 %1:3 %1:3
        figure;    iplot=0; hold on
        % set(gcf, 'Position', [0 -600 210*4 297*4])
        set(gcf, 'Position', [0 -200 375*3 210*4])
        load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
        
        for isoi = 1:3%6 % occ or motor
            
            for idrug = 1:4 %[1,2,4] %
                iplot = iplot+1; subplot(3,4,iplot); hold on
                
                % %                 dum = squeeze(mean(respavglat(:,isoi,iwrt,5,frind,tind,idrug, imotor, idiff, 3))); %average over subj
                dum = squeeze(mean(respavglat(:,isoi,frind,tind,idrug, imotor, idiff, istimorresp, icorrect))); %average over subj
                
                if showstats
                    nsalpha=0.2;
                    pval = squeeze( latstat{isoi, idrug, idiff, icorrect}.prob);
                    pval = pval(frind,:);
                    thrcluster = zeros(size(pval)) + nsalpha; % set non significant alpha;
                    thrcluster(find(pval<0.05)) = 1;
                    tmpcdat = (dum + -ZLIM(1)) * (256 / (-ZLIM(1) + ZLIM(2)));
                    rgbcdat = ind2rgb(uint8(floor(tmpcdat)), colormap); % ind->rgb->hsv ||change saturation values||  hsv->rgb ->  plot
                    hsvcdat = rgb2hsv(rgbcdat);
                    hsvcdat(:,:,2) = hsvcdat(:,:,2) .* thrcluster;
                    dum = hsv2rgb(hsvcdat);
                    % make gray areas lighter to highlight significance borders
                    graypix = hsvcdat(:,:,2) < 0.1;
                    graypix=repmat(graypix, [1 1 3]);
                    dum(graypix) = 0.9;
                end
                
                %                 imagesc(taxis(tind),faxis(frind),dum,ZLIM);
                imagesc(taxis,faxis(frind),dum,ZLIM);
                
                %                 title({sprintf('Lateralization %s, %s', latrleg{isoi},  pharm_conds{idrug} ) ;
                %                     sprintf('%slocked [%g]', trigger, ZLIM(2))}); % 'FontWeight','bold'
                title(sprintf('Ltr %s %s, %s, psc=%g%%', correct_conds{icorrect}, latrleg{isoi},  pharm_conds{idrug}, ZLIM(2)*100)); % 'FontWeight','bold'
                
                hold on
                yaxis = [FREQLO,FREQHI];
                plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                
                xlim(XLIM)
                ylim(YLIM)
                YTICKS = [0:20:200];
                set(gca,'Box','off','XTick',-1:0.2:1,...    [-0.5,0,0.5,1,1.5,2,2.5]
                    'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                    'TickDir','out', 'FontSize', 8);
                LAB=1;
                xlabel(sprintf('Time from %s (s)', trigger));
                ylabel('Frequency (Hz)');
                zlabel('Response (%)');
                %             h=colorbar;
            end
            %         end %ises
        end
        if SAV
            outpath = fullfile(PREOUT, 'poolings');
            warning off; mkdir(outpath); warning on
            outfile = sprintf( '%s%sTFRlateralization_%slocked_%s_%s_%s_%s_%s', outpath, filesep, trigger, TYP, analysistype, motor_conds{imotor}, diff_conds{idiff}, correct_conds{icorrect}  );
            disp(outfile)
            export_fig(outfile, '-pdf', '-transparent',  '-depsc') %'-png',
        end;
    end
end

%% Define the two bands that show the choice effect, based on the placebo TFR: 1. lowest freqbin - 20 Hz; 2. 30 - 65 Hz.
% Just plot the time courses of these modulations for atomoxetine and placebo (again stim and response locked), and compare the drug conditions via cluster permutation test.  

% stats cfg for line plots
freq=[];
freq.time = taxis;
freq.freq = 1;
freq.dimord = 'subj_time';
freq.label = {'custompooling'};

cfg=[];
cfg.channel          = 'p';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
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

% plotting

close all
frind = find(faxis <= YLIM(2) & faxis >= YLIM(1));
tind = find(taxis <= XLIM(2) & taxis >= XLIM(1));

SAV=1;
set(0,'DefaultAxesFontsize',10)
showstats = 1;

bandoi = [5 20; 30 65]; % based on latr wrt choice
bandoileg = {'beta', 'gamma'};
imotor=3; istimorresp = 3; % icorrect=3; %irt= 4 ; icor=3;

% CLR =[0.8500    0.3250    0.0980;
%     0    0.4470    0.7410];
CLR = {'-r' '-b', '-y', '-k'};
YLIMS = [-0.025 0.02; -0.004 0.018]
figure;    iplot=0;
set(gcf, 'Position', [0 -200 375*3 210*4])
for iband = 1:2 % beta, gamma
    foi = find((faxis>=bandoi(iband,1)) & (faxis<=bandoi(iband,2)));
    for idiff= 3 %1:3 %1:3
        
        for isoi = 3 %1:3%6 % occ or motor
            for icorrect = 1:3
                dum = squeeze(mean(respavglat(:,isoi,foi,tind,:, imotor, idiff, istimorresp, icorrect), 3)); %average over foi
                
                iplot = iplot+1; subplot(2,3,iplot); hold on
                clear h
                for idrug = 1:2;
                    h(idrug) = shadedErrorBar(taxis, squeeze(mean(dum(:,:,idrug))), squeeze(std(dum(:,:,idrug)))/sqrt(nsub), CLR{idrug}, 1 ); %style{icorrect}
                end
                xlim(XLIM); plot(XLIM,[0 0],'--k');
                ylim(YLIMS(iband,:))
                if showstats; ctr=0;
                    for idrug = [1,2,4]
                        ctr = ctr+1;
                        freq.powspctrm = squeeze(dum(:,:,idrug));
                        freqzero = freq; %create zero freq to test against
                        freqzero.powspctrm = zeros(size(freq.powspctrm));
                        cfg.numrandomization = 10000;
                        stat = ft_freqstatistics(cfg, freq, freqzero);
                        if idrug < 3
                            plot_sig_intervals(stat, taxis, ctr, h(idrug).mainLine.Color, 5);
                        else
                            plot_sig_intervals(stat, taxis, ctr, [0 0 0] , 5);
                        end
                    end
                end
                yaxis = get(gca, 'ylim'); plot([0,0],yaxis,'k');
                set(gca, 'tickdir', 'out')
                title({sprintf('Lateralization %s', latrleg{isoi}) ...
                    sprintf('%s-band, %s trials', bandoileg{iband}, correct_conds{icorrect} )}); % 'FontWeight','bold'
                
                xlabel(sprintf('Time from %s (s)', trigger));
                ylabel('Response (%)');
                
            end
            legend([h.mainLine], pharm_conds{1:2})
            legend BOXOFF
        end
        
    end
    %         end %ises
end
TYP = '2afc';
if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%sTFRlateralization_bands_%slocked_%s_%s_%s_%s_%s', outpath, filesep, trigger, TYP, analysistype, motor_conds{imotor}, diff_conds{idiff}, correct_conds{icorrect}  );
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent',  '-depsc') %'-png',
    export_fig(outfile, '-png', '-transparent',  '-depsc') %'-png',
end;



