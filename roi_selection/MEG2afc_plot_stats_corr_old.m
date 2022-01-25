% plot sig clusters from 5D freqanalysis of drift rate power correlation,
% integrating over different dims for plotting

% MEG2afc_load_respavg

%% load stats
statfolder = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/stats_5D_correlation';
cd(statfolder)

for idrug = 3:4
    for idiff = 3
        for itrig=1:2
            load( sprintf('stat_idrug%d_idiff%d_itrig%d', idrug, idiff, itrig) );
            allstat{idrug,idiff,itrig} = stat;
        end
    end
end

%% 1) TFR's for clusters p < 0.025, integrating over chans (of cluster)
%% 2) topo, integrating over freq and time of cluster
%% 3) time course of correlation values, integrating over chans and freq
%% 4) power spectrum of correlation, integrating over chans and time

%% plot integrated TFR's for stim and resp low and high freq

close all
THR = 0.025;
trap = 1;
foi = [faxis_all(1) faxis_all(end)];
corrfactor = {'driftrate'};
icorr = 1;
sign = {'Pos', 'Neg'};
isign = 1; %  positive or negagtive clusters
SAV = 1;
ZLIM = [-0.1 0.1];

for idrug = 3:4 %3:4
    figure; iplot = 0;
    set(gcf, 'Position', [0 -200 375*3 210*4])
    
    for idiff = 3
        for ifreq = 2:-1:1 % high, then low freq
            for itrig= 1:2 %1:2
                
                plotstat =  allstat{idrug,idiff,itrig};
                % select freq range
                frind_log = ismember(plotstat.freq, faxis{ifreq});
                plotstat.prob = plotstat.prob(:,frind_log,:);
                plotstat.posclusterslabelmat = plotstat.posclusterslabelmat(:,frind_log,:);
                plotstat.freq = faxis{ifreq};
                
                if isign == 1
                    clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters));
                    cluslabels = plotstat.posclusterslabelmat;
                    pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
                else
                    clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
                    cluslabels = plotstat.negclusterslabelmat;
                    pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
                end
                
                if ~isempty(clus)
                    for iclus = 1%:2 %clus % TODO all clusters at once
                        
                        % for TFR
                        rtemp = plotstat.rho(:,frind_log,:);
                        if trap
                            %                             rtemp(cluslabels~=clus(iclus)) = 0;
                            rtemp = squeeze(trapz(rtemp,1));%./numel(find(rtemp~=0));
                            scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                        else
                            rtemp(cluslabels~=clus(iclus)) = NaN;
                            rtemp = squeeze(nanmean(rtemp,1));%./numel(find(rtemp~=0));
                            rtemp(isnan(rtemp))=0;
                            scale = [-0.05 0.05];
                        end
                        cmap = cbrewer('div', 'RdBu',256);
                        cmap = flipud(cmap);
                        iplot = iplot + 1;
                        subplot(2,2, iplot); hold on
                        colormap(cmap)
                        %                                         rtemp = smooth2a(rtemp,1);
                        
                        %                         thrcluster = squeeze(any(plotstat.posclusterslabelmat == iclus, 1));
                        thrcluster = squeeze(any(ismember(plotstat.posclusterslabelmat, clus), 1));
                        %                         ft_plot_matrix(plotstat.time, plotstat.freq, rtemp, ... 'clim', scale);
                        %                             'highlight', thrcluster, 'highlightstyle', 'saturation', 'clim', scale)
                        
                        rtemp = showstatsTFR(rtemp, thrcluster, 1:length(faxis{ifreq}), scale, 1);
                        
                        imagesc(taxis{itrig}, faxis{ifreq}, rtemp,ZLIM);
                        %                         colorbar
                        for ic = clus
                            thrcluster = squeeze(any(plotstat.posclusterslabelmat == ic, 1));
                            thrcluster = thrcluster * (ic*2);
                            contour(taxis{itrig}, faxis{ifreq}, thrcluster, 1, 'ShowText', 'On');
                        end
                        
                        set(gca,'YDir','normal');
                        if ifreq == 2
                            set(gca, 'YTick', 40:20:150)
                        end
                        set(gca, 'XTick', -1:0.25:1)
                        plot([0,0],foi,'k',[0,0],foi,'k');
                        ylabel('Frequency (Hz)');
                        xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
                        ylim([faxis{ifreq}(1) faxis{ifreq}(end)])
                        if ifreq == 2
                            title(sprintf('%s %s %s-locked',  pharm_conds{idrug}, diff_conds{idiff}, trigger_leg{itrig} ))
                        end
                        
                    end
                end
            end
        end
        if SAV
            outpath = fullfile(PREOUT, 'corr5d');
            warning off; mkdir(outpath); warning on
            outfile = sprintf( '%s%sTFRstats_%sclusters_%s_%s', outpath, filesep, sign{isign},  pharm_conds{idrug}, diff_conds{idiff});
            display(outfile)
            export_fig(outfile, '-pdf') %'-png', ,  '-depsc'  '-transparent'
            
        end
        
    end
end
cd(outpath)

%% plot Hipp style figure for each cluster: TFR and topo

close all
THR = 0.025;
trap = 1;
foi = [faxis_all(1) faxis_all(end)];
corrfactor = {'driftrate'};
icorr = 1;
sign = {'Pos', 'Neg'};
SAV = 0;

subplotinds = [5, 3];

for idrug = 4 %3:4
    for idiff = 3
        for itrig= 2 %1:2 %1:2
            figure
            for ifreq = 2:-1:1 % high, then low freq
                
                plotstat =  allstat{idrug,idiff,itrig};
                frind_log = ismember(plotstat.freq, faxis{ifreq});
                plotstat.prob = plotstat.prob(:,frind_log,:);
                plotstat.posclusterslabelmat = plotstat.posclusterslabelmat(:,frind_log,:);
                plotstat.freq = faxis{ifreq};
                
                isign = 1; % only positive clusters for now
                if isign == 1
                    clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters));
                    cluslabels = plotstat.posclusterslabelmat;
                    pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
                else
                    clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
                    cluslabels = plotstat.negclusterslabelmat;
                    pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
                end
                
                if ~isempty(clus)
                    for iclus = clus
                        % for TFR
                        rtemp = plotstat.rho(:,frind_log,:);
                        if trap
                            %                             rtemp(cluslabels~=clus(iclus)) = 0;
                            rtemp = squeeze(trapz(rtemp,1));%./numel(find(rtemp~=0));
                            scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                        else
                            rtemp(cluslabels~=clus(iclus)) = NaN;
                            rtemp = squeeze(nanmean(rtemp,1));%./numel(find(rtemp~=0));
                            rtemp(isnan(rtemp))=0;
                            scale = [-0.05 0.05];
                        end
                        cmap = cbrewer('div', 'RdBu',256);
                        cmap = flipud(cmap);
                        subplot(3,2, subplotinds(ifreq)); colormap(cmap)
                        %                                         rtemp = smooth2a(rtemp,1);
                        thrcluster = squeeze(any(plotstat.posclusterslabelmat == iclus, 1));
                        %                         ft_plot_matrix(plotstat.time, plotstat.freq, rtemp, ... 'clim', scale);
                        %                             'highlight', thrcluster, 'highlightstyle', 'saturation', 'clim', scale)
                        
                        rtemp = showstatsTFR(rtemp, thrcluster, 1:length(faxis{ifreq}), scale, 1);
                        
                        imagesc(taxis{itrig}, faxis{ifreq}, rtemp,ZLIM);
                        %                         colorbar
                        
                        hold on
                        set(gca,'YDir','normal');
                        plot([0,0],foi,'k',[0,0],foi,'k');
                        ylabel('Frequency (Hz)');
                        if ifreq == 1
                            xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
                        end
                        
                        % for integration plots
                        rtemp = plotstat.rho(:,frind_log,:);
                        if trap
                            rtemp(cluslabels~=clus(iclus)) = 0;
                            rtemp = squeeze(trapz(rtemp,1));%./numel(find(rtemp~=0));
                            scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                        else
                            rtemp(cluslabels~=clus(iclus)) = NaN;
                            rtemp = squeeze(nanmean(rtemp,1));%./numel(find(rtemp~=0));
                            rtemp(isnan(rtemp))=0;
                            scale = [-0.05 0.05];
                        end
                        if unique(rtemp) == 0
                            disp('rtemp empty')
                            continue
                        end
                        rtemp_time = squeeze(trapz(rtemp,1));
                        rtemp_freq = squeeze(trapz(rtemp,2));
                        
                        subplot(3,2,1);
                        area(plotstat.time, rtemp_time); hold on;
                        plot([plotstat.time(1), plotstat.time(end)], [0 0],'k');
                        xlim([plotstat.time(1), plotstat.time(end)]);
                        %                         ylim([-(ceil(max(abs(rtemp_time)))).*1.2, ceil(max(abs(rtemp_time))).*1.2]);
                        ylim([0, ceil(max(abs(rtemp_time))).*1.2]);
                        %                         xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
                        ylabel('Integrated correlation');
                        title(sprintf('%s Corr %s %s-%d\n[%d %d] p=%1.4f', corrfactor{icorr}, pharm_conds{idrug}, sign{isign}, iclus, scale(1), scale(2), pvalues(iclus)));
                        
                        
                        subplot(3,2,subplotinds(ifreq)+1) 
                        area(plotstat.freq, rtemp_freq); hold on;
                        plot([plotstat.freq(1) plotstat.freq(end)],[0 0],'k');
                        xlim([plotstat.freq(1) plotstat.freq(end)]);
                        %                         ylim([-(ceil(max(abs(rtemp_freq)))).*1.2, ceil(max(abs(rtemp_freq))).*1.2]);
                        ylim([0, ceil(max(abs(rtemp_freq))).*1.2]);
                        set(gca,'XAxisLocation','top','YAxisLocation','right')
                        camroll(270)
                        camorbit(0,180)
                        xlabel('Frequency (Hz)');
                        ylabel('Integrated correlation');
                        
                        % for topo
                        % CFG
                        subplot(3,2,2); colormap(cmap);
                        cfg = [];
                        cfg.layout = 'CTF275.lay';
                        cfg.comment = 'no';
                        cfg.marker = 'off';
                        cfg.shading = 'flat';
                        cfg.style = 'straight'; %both
                        cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
                        cfg.markersize = 1;
                        cfg.highlight = 'on';
                        cfg.highlightchannel = chlabel(any(any(plotstat.posclusterslabelmat == iclus, 2),3));
                        
                        rtemp = plotstat.rho;
                        if trap
                            rtemp(cluslabels~=clus(iclus)) = 0;
                            rtemp = squeeze(trapz(rtemp,2));
                            rtemp = squeeze(trapz(rtemp,2));
                            cfg.zlim = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                            %                             cfg.zlim = [0, ceil(max(abs(rtemp(:))))];
                        else
                            rtemp(cluslabels~=clus(iclus)) = NaN;
                            rtemp = squeeze(nanmean(rtemp,2));
                            rtemp = squeeze(nanmean(rtemp,2));
                            rtemp(isnan(rtemp))=0;
                            cfg.zlim = [-0.25 0.25];
                        end
                        
                        freq2 = [];
                        freq2.label = chlabel;
                        freq2.dimord = 'chan';
                        freq2.powspctrm = rtemp;
                        freq2.time = 1;
                        freq2.freq = 1;
%                         warning off;
                        
                        ft_topoplotTFR(cfg,freq2);
                        colorbar
                        
                        if SAV
                            outdir = sprintf('~/data/MEG/freq/%s/plots/MEGsplit/',analysistype);
                            outfile = sprintf('%sHipp_corr_%s-locked_%s_%s_%s%d', outdir, trigger_leg{itrig}, corrfactor{icorr}, pharm_conds{idrug}, sign{isign}, clus(iclus));
                            display(outfile)
                            print('-dpdf',outfile)
                        end
                    end
                end
                
                
                
            end
        end
    end
end








% %% plot topo's correlation for different time points
% toi=[];
% toi(1,:) = [0 0.5 1]; % stimlocked
% % toi(2,:) = [-1 -0.5 0]; % resplocked
% toi(2,:) = [-1 -0.75 -0.5]; % resplocked
% toi_ind= [];
% for it = 1:length(toi)
%     for itrig= 1:2 %1:2
%         toi_ind(itrig, it) = find( toi(itrig, it) == taxis{itrig} );
%     end
% end
% %
% close all
% THR = 0.025;
% trap = 1;
% foi = [faxis_all(1) faxis_all(end)];
% 
% corrfactor = {'driftrate'};
% icorr = 1;
% 
% sign = {'Pos', 'Neg'};
% isign = 1; %  positive or negagtive clusters
% SAV = 1;
% 
% for idrug = 4 %3:4
%     for idiff = 3
%         %         for ifreq = 2:-1:1 % high, then low freq
%         for itrig= 1 %1:2 %1:2
%             
%             plotstat =  allstat{idrug,idiff,itrig};
%             %                 % select freq range
%             %                 frind_log = ismember(plotstat.freq, faxis{ifreq});
%             %                 plotstat.prob = plotstat.prob(:,frind_log,:);
%             %                 plotstat.posclusterslabelmat = plotstat.posclusterslabelmat(:,frind_log,:);
%             %                 plotstat.freq = faxis{ifreq};
% 
%             plotstat.posclusterslabelmat = plotstat.posclusterslabelmat(:,:,  toi_ind(itrig,:));
%             isign = 1; % only positive clusters for now
%             if isign == 1
%                 clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters ));
%                 cluslabels = plotstat.posclusterslabelmat;
% %                 pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
%             else
%                 clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
%                 cluslabels = plotstat.negclusterslabelmat;
%                 pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
%             end
%             
%             if ~isempty(clus)
%                 figure
%                 set(gcf, 'Position', [0 -200 375*3 210*4])
%                 ctr=0;
%                 for iclus = clus % [2,5,6] %
%                     
%                     %                         thrcluster = squeeze(any(ismember(plotstat.posclusterslabelmat, clus), 2)); % collapse freq
%                     
%                     % select time points oi
%                     rtemp = plotstat.rho(:,:,  toi_ind(itrig,:) );  % taxis{itrig}(toi_ind(itrig,:))
%                     if trap
%                         rtemp(cluslabels~=clus(iclus)) = 0;
%                         rtemp = squeeze(trapz(rtemp,2)); %integrate over freq
%                         %                             rtemp = squeeze(trapz(rtemp,2));
%                         cfg.zlim = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
%                     else
%                         rtemp(cluslabels~=clus(iclus)) = NaN;
%                         rtemp = squeeze(nanmean(rtemp,2));
%                         rtemp = squeeze(nanmean(rtemp,2));
%                         rtemp(isnan(rtemp))=0;
%                         cfg.zlim = [-0.25 0.25];
%                     end
%                     
% %                     plotstat.rho = temp;
% %                     plotstat.dimord = 'chan_freq_time';
% %                     plotstat.freq = 1;
% %                     plotstat.time = taxis{itrig}(toi_ind(itrig,:));
%                     plotstat.posclusterslabelmat = any(plotstat.posclusterslabelmat == iclus, 2); % select cluster 1, if any freq is part of the cluster, then include
%                     plotstat.negclusterslabelmat = zeros(size(plotstat.posclusterslabelmat));
%                     
%                     cmap_corr = cbrewer('div', 'RdBu',256);
%                     cmap_corr = flipud(cmap_corr);
%                     % subplot(223); colormap(cmap)
%  
%                     
%                     cfg = [];
%                     cfg.layout = 'CTF275.lay';
%                     cfg.comment = 'no';
%                     cfg.marker = 'off';
%                     cfg.shading = 'flat';
%                     cfg.style = 'straight'; %both
%                     cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
%                     cfg.markersize = 1;
%                     cfg.colormap =  cmap_corr;
%                     cfg.highlight = 'off';
%                     cfg.highlightchannel = chlabel(any(any(plotstat.posclusterslabelmat == iclus, 2),3));
%                     cfg.zlim = [-2 2];
%                     
%                     freq2 = [];
%                     freq2.label = chlabel;
%                     freq2.dimord = 'freq_chan_time';
%                     freq2.powspctrm = shiftdim(rtemp, -1);
%                     freq2.time = toi(itrig,:);
%                     freq2.freq = 1;
%                     warning off;
%                     
%                     for itoi = 1:3
%                         ctr = ctr+1;
%                         subplot(length(clus),3,ctr)
%                         cfg.xlim = [ toi(itrig,itoi)   toi(itrig,itoi) ];
%                         ft_topoplotTFR(cfg, freq2);
%                         title(sprintf('Clus %d, %g s wrt %s, %s', iclus, toi(itrig,itoi), trigger_leg{itrig}, pharm_conds{idrug}))
% %                         colorbar
%                     end
% 
% %                     close all
% %                     cfg = [];
% %                     cfg.alpha  = 0.025;
% %                     cfg.parameter = 'rho';
% %                     cfg.zlim   = [-15 15];
% %                     cfg.layout = 'CTF275.lay';
% %                     cfg.colorbar = 'no';
% %                     cfg.subplotsize    = [5 5];
% %                     cfg.highlightsizeseries = 1;
% %                     %         cfg.highlightsymbolseries = '.';
% %                     cfg.saveaspng = 'no';
% %                     cfg.colormap = cmap_corr;
% %                     cfg.maskparameter = 'posclusterslabelmat';
% %                     cfg.style = 'straight';
% %                     cd(PREOUT)
% %                     
% %                     ft_clusterplot(cfg, plotstat);
% %                     
% %                     colormap(cmap_corr);  %cmap = get(gcf, 'Colormap')
% %                     set(gcf, 'Position', [0 -200 375*3 210*4])
%                     
%                 end
%             end
%         end
%     end
% end
% 
% 
