%% plot trial to trial variance and mean across time for the lateralization over motor and occipital cortex
% drug vs nodrug
% dims latrvar:      subj, mean/var, isoi, iwrt, iband, tind, ipharm, imotor, 3
%  size(latrvar):        18     2     2     2    76    37     2     3     3

close all

SAV = 1;

leftind = {leftoccind leftmotorind};
rightind = {rightoccind rightmotorind};

latsoi = {'occipital', 'motor'};
latrWRT = {'wrttargetloc' 'wrtbp' };
latrleg = {'contra1', 'contra2', 'ipsi1', 'ipsi2', 'contra-ipsi' };
bandoi = [7 13; 15 30; 40 80];
bandleg = {'alpha' 'beta' 'gamma' };
latrvar(:,:,:,:, :,:,:,3,:) = mean(latrvar, 8); % avg across regime
% latrvar(:,3,:,:, :,:,:,:,:) = latrvar(:,1,:,:, :,:,:,:,:) ./ abs(latrvar(:,2,:,:, :,:,:,:,:)) ; %
measureleg = { 'variance' 'mean'};

if strcmp(trigger, 'stim')
    XLIM          = [-0.3 1.3];
else
    XLIM         = [-0.75 0.25];
end

% freqstatistics on variance and power:
% test mean and variance: drug vs placebo
% motor and occ lateralization wrt motor and visual
%
nsub=18;
% respavgstat = squeeze(mean(respavg(:,:,sensind,:,:,:,:),3));

cfg = [];
cfg.frequency        = 'all';
cfg.latency          = XLIM;
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
freq2stats.dimord = 'subj_time'; % freq_
freq2stats.label = {'custompooling'};

% dims latrvar:      subj, var/mean, isoi, iwrt, iband, tind, ipharm, imotor, idiff
%  size(latrvar):    18     2         2     2    76     37    2       3       3

stat=[];
for isoi= 1:2 %1:2 % occ and motor
    if isoi==1; iwrt = 2; end
    if isoi==2; iwrt = 1; end
    
    for idiff = 1:3
        for iband = 1:length(bandoi)
            for idrug=1:2
                freq2stats.powspctrm = squeeze(latrvar(:,1,isoi,iwrt,length(faxis)+iband,:, idrug,3,idiff)); % collapsed over [regime diff stimloc]
                freq2statszero = freq2stats; %create zero freq to test against
                freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
                stat{isoi, iband, idrug, idiff} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
                %                     stat{isoi, iband, idrug, idiff} = ft_timelockstatistics(cfg, freq2stats, freq2statszero);
            end
            % test motor and occ lateralization drug vs placebo
%             freq2stats.powspctrm =  squeeze(latrvar(:,1,isoi,iwrt,length(faxis)+iband,:, 1,3,idiff)); % drug, collapsed over [regime diff stimloc]
%             freq2stats2 = freq2stats;
%             freq2stats2.powspctrm = squeeze(latrvar(:,1,isoi,iwrt,length(faxis)+iband,:, 2,3,idiff)); % placebo, collapsed over [regime diff stimloc]
%             stat{isoi, iband, 4, idiff} = ft_freqstatistics(cfg, freq2stats, freq2stats2);
            
            %test difference between drug and plac against each other to
            %remove within subject variance
            freq2stats.powspctrm =  squeeze(latrvar(:,1,isoi,iwrt,length(faxis)+iband,:, 1,3,idiff)) - ...
                squeeze(latrvar(:,1,isoi,iwrt,length(faxis)+iband,:, 2,3,idiff)); % drug, collapsed over [regime diff stimloc]
            stat{isoi, iband, 4, idiff} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
            
            % % %             alldat(1,:,:) = freq2stats.powspctrm;
            % % %             alldat(2,:,:) = freq2stats2.powspctrm;

        end
    end
end



%%
% dims latrvar: subj var/mean isoi iwrt freq time drug regime
close all
figure; iplot = 0;
set(gcf, 'Position', [0 0 375*3 210*4])

showstats = 1;

for idiff = 1:2
    for isoi = 1:2
        %     for iwrt=1:2
        % if isoi==1; iwrt = 1; end
        % if isoi==2; iwrt = 2; end
        
        for iband = 1:length(bandoi)
            iplot = iplot+1; subplot(4,3,iplot); hold on
            %                 dum = squeeze(mean(latrvar(:,imeasure, isoi, isoi,   length(faxis)+iband,:,:,3,idiff) ));
            %                 plot(taxis, dum , 'Linewidth', 2)
            
            dum = squeeze(mean(latrvar(:,:, isoi, isoi,   length(faxis)+iband,:,:,3,idiff) ));
            
            [AX, H1, H2]=plotyy(taxis, squeeze(dum(2,:,:)), ...
                taxis, squeeze(dum(1,:,:)));
            set(H2, 'LineStyle', '--', 'Linewidth', 2)
            set(H2(1), 'Color', 'r', 'Linewidth', 2); set(H2(2), 'Color', 'b', 'Linewidth', 2) %; set(H2(3), 'Color', 'r'
            set(H1(1), 'Color', 'r', 'Linewidth', 2); set(H1(2), 'Color', 'b', 'Linewidth', 2) % ; set(H1(3), 'Color', 'r'
            %                 legend({'alphavar' 'betavar' 'gammavar' 'alphapow' 'betapow' 'gammapow'})
            ylabel(AX(1), 'Power (solid lines)')
            ylabel(AX(2), 'Variance (dotted lines)')
            plot([0 0], get(gca, 'YLim'), 'k')
            %                 plot([0 0], get(gca, 'YLim'), 'k')
            plot(AX(1), [taxis(1) taxis(end)], [0 0], 'k' )
            xlabel(sprintf('Time from %s (s)', trigger))
            %             set(gca(1), 'TickDir', 'out')
            xlim(AX(1), XLIM);
            xlim(AX(2), XLIM)
            if iplot ==1
                legend(pharm_conds{1:2})
                legend BOXOFF
            end
            set(AX, 'TickDir', 'out', 'XTick', [-2:0.5:2])
            title({sprintf('%s lateralization %s %s %s', bandleg{iband}, latsoi{isoi}, latrWRT{isoi}) ...
                sprintf('%s %slocked', diff_conds{idiff}, trigger) })
            
            if showstats
                BARPOS= [0 0 0 ];
                SIGCOL= ['r' 'b' 'k' 'k'];
                THR=0.025;
%                 THR=1;
                BARWID= 3;
                
                for idrug= [1,2,4]  % 4,2,1
                    pval = stat{isoi, iband, idrug, idiff}.prob';  % stat from ft_freqstatistics
                    
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
                            
                            axylim = get(AX(1), 'YLim');
                            range = axylim(2) - axylim(1);
                            ylimsteps = axylim(1):range/20:axylim(2);
                            
                            plot( AX(1), stat{isoi, iband, idrug, idiff}.time(begsmp:endsmp),...   %begsmp-1:endsmp+1
                                ylimsteps(end-idrug)*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
                                'Color',SIGCOL(idrug),'LineWidth',BARWID)
%                             plot( AX(1), stat{isoi, iband, idrug, idiff}.time(begsmp:endsmp),...   %begsmp-1:endsmp+1
%                                 BARPOS(idrug)*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
%                                 'Color',SIGCOL{idrug},'LineWidth',BARWID)
                        end
                    end
                end %istat            
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'variance');
    warning off; mkdir(outpath); warning on
    
    outfile = sprintf( '%s%sTFRlateralization_MeanvsVar_%slocked', outpath, filesep, trigger );
    disp(outfile)
    %             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
    %     print('-dpng', outfile) %, '-cmyk'
    %     orient landscape
    %     print('-depsc', outfile) %, '-cmyk'
    %     saveas(gcf, [outfile '.pdf'])
    export_fig(outfile, '-pdf', '-transparent') % '-png',
end;

% %% Divide var by the mean and plot
%
% figure; iplot = 0;
% set(gcf, 'Position', [0 0 297*4 210*3])
%
% for idiff = 1:2
%     for isoi = 1:2
%         %     for iwrt=1:2
%         % if isoi==1; iwrt = 1; end
%         % if isoi==2; iwrt = 2; end
%
%         for iband = 2:length(bandoi)
%             for imeasure=1 %  var or mean
%                 iplot = iplot+1; subplot(2,4,iplot); hold on
%                 %                 dum = squeeze(mean(latrvar(:,imeasure, isoi, isoi,   length(faxis)+iband,:,:,3,idiff) ));
%                 %                 plot(taxis, dum , 'Linewidth', 2)
%
%                 dum = squeeze(mean(latrvar(:,3, isoi, isoi,   length(faxis)+iband,:,:,3,idiff) ));
%                 plot(taxis, dum , 'Linewidth', 2)
%
% %                 [AX, H1, H2]=plotyy(taxis, squeeze(dum(2,:,:)), ...
% %                     taxis, squeeze(dum(1,:,:)));
% %                 set(H2, 'LineStyle', '--', 'Linewidth', 2)
% %                 set(H2(1), 'Color', 'r', 'Linewidth', 2); set(H2(2), 'Color', 'b', 'Linewidth', 2) %; set(H2(3), 'Color', 'r'
% %                 set(H1(1), 'Color', 'r', 'Linewidth', 2); set(H1(2), 'Color', 'b', 'Linewidth', 2) % ; set(H1(3), 'Color', 'r'
% %                 %                 legend({'alphavar' 'betavar' 'gammavar' 'alphapow' 'betapow' 'gammapow'})
% %                 ylabel(AX(1), 'Power (solid lines)')
% %                 ylabel(AX(2), 'Variance (dotted lines)')
% %                 plot([0 0], get(gca, 'YLim'), 'k')
% %                 %                 plot([0 0], get(gca, 'YLim'), 'k')
% %                 plot(AX(1), [taxis(1) taxis(end)], [0 0], 'k' )
% %                 %             xlabel('Time (s)')
% %                 %             set(gca(1), 'TickDir', 'out')
% %                 xlim(AX(1), XLIM);
% %                 xlim(AX(2), XLIM)
%                 if iplot ==1
%                     legend(pharm_conds{1:2})
%                     legend BOXOFF
%                 end
%
%
% %                 set(AX, 'TickDir', 'out', 'XTick', [-2:0.5:2])
%                 title({sprintf('%s lateralization %s %s %s', bandleg{iband}, latsoi{isoi}, latrWRT{isoi}) ...
%                     sprintf('%s %slocked', diff_conds{idiff}, trigger) })
%
%             end
%         end
%     end
% end
% if SAV
%     %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%     outpath = fullfile(PREOUT, 'variance');
%     warning off; mkdir(outpath); warning on
%
%     outfile = sprintf( '%s%sTFRlateralization_Fanofactor_%slocked', outpath, filesep, trigger );
%     disp(outfile)
%     %             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
%     %     print('-dpng', outfile) %, '-cmyk'
%     %     orient landscape
%     %     print('-depsc', outfile) %, '-cmyk'
%     %     saveas(gcf, [outfile '.pdf'])
%     export_fig(outfile, '-pdf', '-transparent') % '-png',
% end;

