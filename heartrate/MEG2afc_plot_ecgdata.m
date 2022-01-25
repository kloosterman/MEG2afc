%% MEG2afc_plot_dva: plot dva Schurger
% function plotTFR_2AFC(~)
% plots SOI-TFRs
% average of selected sensors for multiple subjects
% sorted according to stimulus condition
% % % % % addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
% % % % % addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
% % % % % ft_defaults
cd('/Users/kloosterman/gridmaster2012/MATLAB/MEG_HH_analysis/plotting')

% MEG2afc_setup_paths
% trigger = 'resp'

MEG2afc_load_ecgdata

%% do stats
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
timelock=[];
timelock.time = taxis;
timelock.dimord = 'chan_subj_time';
timelock.label = {'custompooling'};

cfg = [];
% cfg.latency          = XLIM;
cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
cfg.statistic        = 'indepsamplesT'; % we compare groups here
cfg.correctm         = 'cluster';
% % % %     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
cfg.neighbours       = []; %in case no channel data present


% poolstat=[];
% for isoi=1%:2 % occ and motor
%     for idva=1:4
%         for idrug = [1,2,4] % 1:4
%             for idiff = 3 %1:4 %1:4
%                 timelock.trial = squeeze(respavg(:,idva,:,idrug, 3, idiff, 3, 3 )); %average over subj
%                 timelock.trial = shiftdim(timelock.trial,-1);
%                 timelockzero = timelock; %create zero timelock to test against
%                 timelockzero.trial = zeros(size(timelock.trial));
%                 poolstat{idva, idrug, idiff} = ft_timelockstatistics(cfg, timelock, timelockzero);
%             end
%         end
%     end
% end
% for respavg_ses:
poolstat=[];
istim = 3;
iresp = 3;
for idiff = 3 %1:4 %1:4
    for isoi = 1 %:3%:3%:2 % occ and motor
        for idva= [1:3,5]% 1:4
            for ises = 1:4
                %                 for idrug = 1:2
                subj_drug_ind = find(drug_by_order_mat(:,ises)==1);
                subj_plac_ind = find(drug_by_order_mat(:,ises)==2);
                %                         subj_motor_ind = find(motor_by_order_mat(:,ises)==idrug-1);
                
                dum_drug = squeeze(respavg_ses(subj_drug_ind, idva,isoi, tind, ises, idiff, istim, iresp ));
                dum_plac = squeeze(respavg_ses(subj_plac_ind, idva,isoi, tind, ises, idiff, istim, iresp ));
                
                timelock.trial = dum_drug;
                timelock.trial = shiftdim(timelock.trial,-1);
                
                timelocktwo = timelock; %create zero timelock to test against
                
                timelocktwo.trial = dum_plac;
                timelocktwo.trial = shiftdim(timelocktwo.trial,-1);
                %
                %                     design = zeros(2,2*nsub);
                %                     for i = 1:nsub
                %                         design(1,i) = i;
                %                     end
                %                     for i = 1:nsub
                %                         design(1,nsub+i) = i;
                %                     end
                %                     design(2,1:nsub)        = 1;
                %                     design(2,nsub+1:2*nsub) = 2;
                %
                design = zeros(1,size(timelock.trial,2) + size(timelocktwo.trial,2));
                design(1,1:size(timelock.trial,2)) = 1;
                design(1,(size(timelock.trial,2)+1):(size(timelock.trial,2) + size(timelocktwo.trial,2)))= 2;
                
                %                 design = zeros(2,size(timelock.trial,2) + size(timelocktwo.trial,2));
                %                 design(1,1:size(timelock.trial,2)) = 1:size(timelock.trial,2);
                %                 design(1,(size(timelock.trial,2)+1):(size(timelock.trial,2) + size(timelocktwo.trial,2)))= 1:size(timelocktwo.trial,2);
                %                 design(2,1:size(timelock.trial,2)) = 1;
                %                 design(2,(size(timelock.trial,2)+1):(size(timelock.trial,2) + size(timelocktwo.trial,2)))= 2;
                
                cfg.design   = design;
                %                     cfg.uvar     = 1;
                cfg.ivar     = 1;
                
                poolstat{ises, idva, idrug, isoi} = ft_timelockstatistics(cfg, timelock, timelocktwo);
                %                 end
            end
        end
    end
end
% save(fullfile(PREIN, 'dva', outfile), 'poolstat');

%% plotting: 2AFC heartrate drug vs placebo, collapsed over ssesions
close all
SAV=1;
set(0,'DefaultAxesFontsize',14)

for isoi = 1 %[1,3]
    for idiff = 3%1:4
        figure;    iplot=0; hold on
        set(gcf, 'Position', [0 -200 375*3 210*4])
        for ivar = 1:2 % [1,3:5] %1:4
            iplot=iplot+1; subplot(1,2,iplot); hold on
            ctr=0; clear h; clear legtxt
            ctr=ctr+1;
            if ivar == 1
                dum = squeeze(respavg(:,9, 1:2, 3 )); %average over subj
                dumstd = squeeze(nanstd(respavg(:,9, 1:2, 3  ))) / sqrt(nsub); %average over subj
                legtxt = pharm_conds(1:2);
            else
                dum = squeeze(respavg(:,9, 3, 1:2 )); %average over subj
                dumstd = squeeze(nanstd(respavg(:,9, 3, 1:2  ))) / sqrt(nsub); %average over subj
                legtxt = motor_conds(1:2);
            end
            
            barweb(mean(dum), dumstd', 0.75, [], [], [], 'Heartrate (bpm)', [1 0 0; 0 0 1]) %, pharm_conds(1:idrug)
            ylim([40 80])
            axis square
            
%             [~,pval] = ttest(dum(:,1), dum(:,2));
            pval = randtest1d(dum(:,1), dum(:,2), 0, 10000);
            
            title({sprintf('mean heartrate') ;
                sprintf('%g vs %g, p = %.2g', mean(dum), pval )}); % 'FontWeight','bold'
            legend(legtxt, 'Location', 'NorthEast'); legend BOXOFF
            
        end %
        
        if SAV
            outpath = fullfile(PREOUT, 'dva');
            warning off; mkdir(outpath); warning on
            outfile = sprintf( '%s%sheartrateoverall', outpath, filesep ); % motor_conds{imotor}, pupilconds{ipup}
            disp(outfile)
            export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
        end;
    end %
end

%% Barplots mean heartrate per session
% % TODO compute from fixation color change?
% if strcmp(trigger, 'stim')
%     TIMLO = 0.4;
%     TIMHI = 0.6;
% else
%     TIMLO = -0.6;
%     TIMHI = 0;
% end
% respavg_tim = squeeze(nanmean(respavg_ses(:,:,:, taxis >= TIMLO & taxis <= TIMHI, :,:,:,: ), 4 ));

ses_conds = {'Session 1', 'Session 2', 'Session 3', 'Session 4' };

close all
SAV = 1;
set(0,'DefaultAxesFontsize',14)
dvaleg = {'acrossdva', 'acrossnorm', 'withindva', 'withinnorm', 'withindva-baseline'};

showstats = 0;
imotor =3;
istim=3;
iresp=3;

plotdrug_or_motor = 1;
drug_motor_leg = {'drug' 'motor'};

% if strcmp(trigger, 'stim')
%     XLIM = [-0.5 1.5];
% else
%     XLIM = [-0.75 0.5];
% end
% % YLIM_norm = [   0.2800e-11    0.3500e-11];
% transparent = 1;
% linecols = {'r', 'b', 'g', 'k' ; '--r', '--b', '--g', '--k' ; 'r', 'b', 'g', 'k'};
% % linecols =  cool(4);

TYP = '2afc';
% for isoi = 1%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
for idiff = 3%1:2%1:4
    figure;    
    iplot=0; hold on
    set(gcf, 'Position', [0 -200 375*3 210*4])
    %     for idva = [3,1,2] % 1:3 %[1,2] %[1, 3:5] % 1:3%:4
    for idva = 1 % [1, 3, 5, 2]% [3,5,4] %[3,5] %[1, 3:5] % 1:3%:4

        for isoi= 1 % [1,3] %1:3 %3 %1:3
            dum = []; dumstd = [];
            iplot=iplot+1; subplot(2,2,iplot); hold on
            ctr=0; clear h; clear legtxt
            
            for ises = 1:4 %1:2%:2%:2
                if plotdrug_or_motor == 1
                    for idrug = 1:2
                        ctr=ctr+1;
                        subj_drug_ind = find(drug_by_order_mat(:,ises)==idrug);
                        dum(ises, idrug) = squeeze(nanmean(respavg_ses(subj_drug_ind, 9, ises ))); %average over subj
                        dumstd(ises, idrug) = squeeze(nanstd(respavg_ses(subj_drug_ind, 9, ises ))) / sqrt(length(subj_drug_ind)); %average over subj                        
                    end
                    
                else
                    for imotor = 1:2
                        ctr=ctr+1;                       
                        
                        subj_motor_ind = find(motor_by_order_mat(:,ises)==imotor);
%                         dum = squeeze(nanmean(respavg_ses(subj_motor_ind, idva,isoi, tind, ises, idiff, istim, iresp ))); %average over subj
%                         dumstd = squeeze(nanstd(respavg_ses(subj_motor_ind, idva,isoi, tind, ises, idiff, istim, iresp ))) / sqrt(length(subj_motor_ind)); %average over subj
                        
                        dum(ises, imotor) = squeeze(nanmean(respavg_ses(subj_motor_ind, 9, ises ))); %average over subj
                        dumstd(ises, imotor) = squeeze(nanstd(respavg_ses(subj_motor_ind, 9, ises ))) / sqrt(length(subj_motor_ind)); %average over subj

            
                    end
                end
            end
            
            barweb(dum, dumstd, 0.75, ses_conds, [], 'Session #', 'Heartrate (bpm)', [1 0 0; 0 0 1]) %, pharm_conds(1:idrug)
            
            for ises = 1:4
                if plotdrug_or_motor == 1
                    dat1 = squeeze(respavg_ses( find(drug_by_order_mat(:,ises)==1) , 9, ises ));
                    dat2 = squeeze(respavg_ses( find(drug_by_order_mat(:,ises)==2) , 9, ises ));
                else
                    dat1 = squeeze(respavg_ses( find(motor_by_order_mat(:,ises)==1), 9, ises ));
                    dat2 = squeeze(respavg_ses( find(motor_by_order_mat(:,ises)==2), 9, ises ));
                end
                                    
                [~,pvals(ises)] = ttest2(dat1, dat2);
                
%                 dat1 = dat1(~isnan(dat1));
%                 dat2 = dat2(~isnan(dat2));
%                 pvals(ises) = randtest1d(dat1, dat2, 0, 1000);
                
                
                disp(length(find(drug_by_order_mat(:,ises)==1)))
                disp(length(find(drug_by_order_mat(:,ises)==2)))
            end
            
            title({sprintf('mean heartrate') ;
                sprintf('sig: %.2g %.2g %.2g %.2g ', pvals' )}); % 'FontWeight','bold'
            hold on
            if iplot == 1
                if plotdrug_or_motor == 1
                    legend(pharm_conds{1:idrug}, 'Location', 'Southeast')
                else
                    legend(motor_conds{1:imotor})
                end
            end
%             if idva==1
                ylim([40 90])
%             elseif idva==2
% %                 ylim([-0.18 0])
%             elseif idva==5
%                 ylim([-0.18 0])
%             else
% %                 ylim([0.25 0.7])
%             end
            
            %             ylabel('dva');
            %             set(gca,'Box','on', ... 'XTick',-1:0.2:1,...    [-0.5,0,0.5,1,1.5,2,2.5]
            %                 'YDir','normal', ...'YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
            %                 'TickDir','out', 'FontSize', 12);
            % %             xlabel('Time (s)');
            % %             legend([h.mainLine], legtxt, 'Location', 'SouthEast'); legend BOXOFF
            % %             plot([0,0], get(gca, 'ylim'),'k' );
        end
    end %
    if SAV
        outpath = fullfile(PREOUT, 'dva');
        warning off; mkdir(outpath); warning on
        outfile = sprintf( '%s%sHeartrate%s_persession', outpath, filesep, drug_motor_leg{plotdrug_or_motor} ); % motor_conds{imotor}, pupilconds{ipup}
        disp(outfile)
        export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
    end;
    
end %




%% Barplots accumulated dva from stim onset until report
% TODO compute from fixation color change?
ses_conds = {'Session 1', 'Session 2', 'Session 3', 'Session 4' };

close all
SAV = 1;
set(0,'DefaultAxesFontsize',14)
dvaleg = {'rawwithindva', 'withindva-baseline', 'rawbaseline withindva'};

showstats = 0;
imotor =3;
istim=3;
iresp=3;

plotdrug_or_motor = 1;
drug_motor_leg = {'drug' 'motor'};

if strcmp(trigger, 'stim')
    XLIM = [-0.5 1.5];
else
    XLIM = [-0.75 0.5];
end
% YLIM_norm = [   0.2800e-11    0.3500e-11];
transparent = 1;
linecols = {'r', 'b', 'g', 'k' ; '--r', '--b', '--g', '--k' ; 'r', 'b', 'g', 'k'};
% linecols =  cool(4);

TYP = '2afc';
% for isoi = 1%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
for idiff = 3%1:2%1:4
    figure;    iplot=0; hold on
    set(gcf, 'Position', [0 -200 375*3 210*4])
    for idva = [3,1,2] % 1:3 %[1,2] %[1, 3:5] % 1:3%:4
        for isoi= [1,3] %1:3 %3 %1:3
            dum = []; dumstd = [];
            iplot=iplot+1; subplot(3,2,iplot); hold on
            ctr=0; clear h; clear legtxt
            
            for ises = 1:4 %1:2%:2%:2
                if plotdrug_or_motor == 1
                    for idrug = 1:2
                        ctr=ctr+1;
                        subj_drug_ind = find(drug_by_order_mat(:,ises)==idrug);
                        dum(ises, idrug) = squeeze(nanmean(respavg_ses_accum(subj_drug_ind, idva,isoi, ises, idiff, istim, iresp ))); %average over subj
                        dumstd(ises, idrug) = squeeze(nanstd(respavg_ses_accum(subj_drug_ind, idva,isoi, ises, idiff, istim, iresp ))) / sqrt(length(subj_drug_ind)); %average over subj                        
                    end
                    
                else
                    for imotor = 1:2
                        ctr=ctr+1;                       
                        
                        subj_motor_ind = find(motor_by_order_mat(:,ises)==imotor);
%                         dum = squeeze(nanmean(respavg_ses(subj_motor_ind, idva,isoi, tind, ises, idiff, istim, iresp ))); %average over subj
%                         dumstd = squeeze(nanstd(respavg_ses(subj_motor_ind, idva,isoi, tind, ises, idiff, istim, iresp ))) / sqrt(length(subj_motor_ind)); %average over subj
                        
                        dum(ises, imotor) = squeeze(nanmean(respavg_ses_accum(subj_motor_ind, idva,isoi, ises, idiff, istim, iresp ))); %average over subj
                        dumstd(ises, imotor) = squeeze(nanstd(respavg_ses_accum(subj_motor_ind, idva,isoi, ises, idiff, istim, iresp ))) / sqrt(length(subj_motor_ind)); %average over subj

            
                    end
                end
            end
            
            barweb(dum, dumstd, 0.75, ses_conds, [], 'Session #', 'dva', [1 0 0; 0 0 1]) %, pharm_conds(1:idrug)
            
            for ises = 1:4
                if plotdrug_or_motor == 1
                    dat1 = squeeze(respavg_ses_accum( find(drug_by_order_mat(:,ises)==1), idva,isoi, ises, idiff, istim, iresp ));
                    dat2 = squeeze(respavg_ses_accum( find(drug_by_order_mat(:,ises)==2), idva,isoi, ises, idiff, istim, iresp ));
                else
                    dat1 = squeeze(respavg_ses_accum( find(motor_by_order_mat(:,ises)==1), idva,isoi, ises, idiff, istim, iresp ));
                    dat2 = squeeze(respavg_ses_accum( find(motor_by_order_mat(:,ises)==2), idva,isoi, ises, idiff, istim, iresp ));
                end
                                    
                [~,pvals(ises)] = ttest2(dat1, dat2);
                disp(length(find(drug_by_order_mat(:,ises)==1)))
                disp(length(find(drug_by_order_mat(:,ises)==2)))
            end
            
            title({sprintf('mean dva from stimonset to resp' ) ;
                sprintf('%s %s, sig: %.2g %.2g %.2g %.2g ', dvaleg{idva}, sens.leg{isoi}, pvals' )}); % 'FontWeight','bold'
            hold on
            if plotdrug_or_motor == 1
                legend(pharm_conds{1:idrug})
            else
                legend(motor_conds{1:imotor})
            end
            if idva==1
%                 ylim([0.2 0.35])
            elseif idva==2
                ylim([-0.18 0])
            else
                ylim([0.25 0.7])
            end
            
            %             ylabel('dva');
            %             set(gca,'Box','on', ... 'XTick',-1:0.2:1,...    [-0.5,0,0.5,1,1.5,2,2.5]
            %                 'YDir','normal', ...'YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
            %                 'TickDir','out', 'FontSize', 12);
            % %             xlabel('Time (s)');
            % %             legend([h.mainLine], legtxt, 'Location', 'SouthEast'); legend BOXOFF
            % %             plot([0,0], get(gca, 'ylim'),'k' );
        end
    end %
    if SAV
        outpath = fullfile(PREOUT, 'dva');
        warning off; mkdir(outpath); warning on
        outfile = sprintf( '%s%saccum_dva_%slocked_%s_%s_%s_persession', outpath, filesep, trigger, TYP, sens.leg{isoi}, diff_conds{idiff},  drug_motor_leg{plotdrug_or_motor} ); % motor_conds{imotor}, pupilconds{ipup}
        disp(outfile)
        export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
    end;
    
end %





%
%
%
% %% Plot bars for the frontal gamma: 49-61 Hz, -0.4 0.2 s
%
% SAV=1;
% % close all
% isoi=2;
% frindoi = find(faxis <= 61 & faxis >= 49);
% % taxis(tind) : -0.4 to 0.2
% idrug = 1:2;
% idiff = 1:2;
% dum = squeeze(respavgstat(:,isoi,frindoi,tind,idrug, 3, idiff, 3, 3)); %average over subj
% dum=mean(dum,2);
% dum=squeeze(mean(dum,3)); % subj drug diff remain
%
% sem = squeeze(std(dum) / sqrt(18));
% pvaleasy = randtest1d(squeeze(dum(:,1,1)), squeeze(dum(:,2,1)), 0, 1000);
% pvalhard = randtest1d(squeeze(dum(:,1,2)), squeeze(dum(:,2,2)), 0, 1000);
%
% figure
% barweb(squeeze(mean(dum))', sem', [], {'Easy' 'Hard'})
% % legend(diff_conds(idiff))
% legend(pharm_conds(idrug))
% legend BOXOFF
% title(sprintf('Easy: p = %g, Hard: p = %g', pvaleasy, pvalhard))
% shading flat
% ylabel('Modulation (%)')
%
% if SAV
%     %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%     outpath = fullfile(PREOUT, 'bargraphs');
%
%     warning off; mkdir(outpath); warning on
%     outfile = sprintf( '%s%sGammadrugeffect_EasyVShard_%slocked_%s_%s_%sfreq_%s', outpath, filesep, trigger, TYP, sens.leg{isoi}, analysistype, motor_conds{imotor}  );
%     disp(outfile)
%     %             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
%     %         print('-dpng', outfile) %, '-cmyk'
%     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
% end;
%
% %% plot spectrum of stimlocked activity: 0.25 to 0.5 s
%
% cfg = [];
% cfg.frequency        = 'all';
% % cfg.latency          = [TIMLO TIMHI];
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
% cfg.correctm         = 'cluster';
% %     cfg.correctm         = 'no';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 1000;
% cfg.neighbours       = []; %in case no channel data present
%
% design = zeros(2,2*nsub);
% for i = 1:nsub
%     design(1,i) = i;
% end
% for i = 1:nsub
%     design(1,nsub+i) = i;
% end
% design(2,1:nsub)        = 1;
% design(2,nsub+1:2*nsub) = 2;
%
% cfg.design   = design;
% cfg.uvar     = 1;
% cfg.ivar     = 2;
%
% freq2stats = freq;
% freq2stats.dimord = 'subj_freq';
% freq2stats.label = {'custompooling'};
% freq2stats.time = 1;
%
% SAV=1;
% close all
% if strcmp(trigger, 'stim')
%     toi = [0.25 0.5];
% else
%     toi = [-0.5 0];
% end
% toiind = find(freq.time <= toi(2) & freq.time >= toi(1) );
%
% figure; iplot=0;
% set(gcf, 'Position', [0 0 375*2 210*4])
% for isoi = 1:2
%     iplot=iplot+1; subplot(2,1,iplot)
%     dum = squeeze(mean(respavgstat(:,isoi,frind,toiind, 3, 3, 3, 3, 3), 4)); %average over tind
%
%     freq2stats.powspctrm = dum; %DIMS:
%     freq2statszero = freq2stats; %create zero freq to test against
%     freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
%     stat = ft_freqstatistics(cfg, freq2stats, freq2statszero);
%
%     avg = squeeze(mean(dum));
%     sem = std(dum) / sqrt(18);
%
%     [~, minind]=min(avg);
%     [~, maxind]=max(avg(5:end));
%
%     h = shadedErrorBar(faxis,avg,sem);
%     hold on
%     plot([0 150], [0 0],'k')
%     set(gca, 'Tickdir', 'out', 'XTick', 0:10:150)
%     set(h.mainLine, 'Linewidth', 2)
%     grid ON
%     xlabel('Frequency (Hz)')
%     ylabel('Modulation (%)')
%     title(sprintf('Power spectrum %s across all conditions; interval %g-%g s %s-locked; min @ %d Hz, max @ %d Hz', sens.leg{isoi}, toi, trigger, faxis(minind), faxis(maxind) ))
%     ylim([-0.5 0.2])
%
%     pval = stat.prob';  % stat from ft_freqstatistics
%     BARPOS= [-0.4];
%     SIGCOL= ['k'];
%     THR=0.025;
%     %                 THR=1;
%     BARWID= 5;
%     sigint={};
%     i = 1; sigint{i} = [];
%     for ismp = 1:length(pval)
%         if pval(ismp) < THR
%             sigint{i} = [sigint{i} ismp];
%         end % concatenate significant samples in sigint{iint}
%         if (ismp < size(pval,1)) && (pval(ismp+1) >= THR) && ~isempty(sigint{i}) % jump to next interval if next sample is not sig
%             i= i + 1;
%             sigint{i} = [];
%         end
%     end
%     % replot significant intervals in different color
%     for iint = 1:length(sigint)
%         if ~isempty(sigint{iint})
%             begsmp = sigint{iint}(1); %*2
%             endsmp = sigint{iint}(end); %*2
%
%             %             axylim = get(AX(1), 'YLim');
%             %             range = axylim(2) - axylim(1);
%             %             ylimsteps = axylim(1):range/20:axylim(2);
%
%             %             plot( AX(1), stat{isoi, iband, idrug, idiff}.time(begsmp:endsmp),...   %begsmp-1:endsmp+1
%             %                 ylimsteps(end-idrug)*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
%             %                 'Color',SIGCOL,'LineWidth',BARWID)
%             plot(stat.freq(begsmp:endsmp), ...   %begsmp-1:endsmp+1
%                 BARPOS*(ones(1,length(begsmp:endsmp))), ...  %begsmp-1:endsmp+1
%                 'Color',SIGCOL,'LineWidth',BARWID)
%         end
%     end
%
%
%
% end
% if SAV
%     %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%     outpath = fullfile(PREOUT, 'poolings');
%     warning off; mkdir(outpath); warning on
%
%     outfile = sprintf( '%s%sSpectrumAllconds_%slocked', outpath, filesep, trigger );
%     disp(outfile)
%     %             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
%     %     print('-dpng', outfile) %, '-cmyk'
%     %     orient landscape
%     %     print('-depsc', outfile) %, '-cmyk'
%     %     saveas(gcf, [outfile '.pdf'])
%     export_fig(outfile, '-png', '-transparent') % '-png',
% end;
%
%
% %% Plot accumulated dva from stim onset to report
%
% showstats = 0;
% imotor =3;
% istim=3;
% iresp=3;
% isoi = 3;
% idiff = 3;
%
% figure;iplot=0;
% for inorm=1 %1:2  % noBC, BC
%     iplot = iplot+1;
%     subplot(2,2, iplot)
%     dat = squeeze(respavg_accum(:,inorm,isoi, 1:2, imotor, idiff, istim, iresp));
%     barweb(nanmean(dat), nanstd(dat)/sqrt(nsub))
%     ylabel('dva')
%     legend(pharm_conds(1:2))
% end
%
%
%
%
%
%
%
%
