% plot sig clusters from 5D freqanalysis of drift rate power correlation,
% integrating over different dims for plotting

% MEG2afc_load_respavg
% restoredefaultpath
if ismac
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman'; % on the cluster
    backend = 'none'; % local torque
else
    basepath = '/home/mpib/kloosterman'; % on the cluster
    addpath(fullfile(basepath, 'MATLAB',  'tools/qsub_tardis'))
    backend = 'torque'; % local torque
end
addpath(genpath(fullfile(basepath, 'MATLAB',   'MEG_HH_analysis')))
rmpath(genpath(fullfile(basepath, 'MATLAB', 'MEG_HH_analysis/.git/')))
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip-20161220')) %inc JJ edit ft_artifact_zvalue
ft_defaults
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip-20161220/external/spm8/')) %for spm_bwlabel

%%
addpath(genpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting'))
addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/cbrewer')

% Control panel
testdv = {'modulation'; 'correlation'}; % test modulation or power-drift correlation
itest = 1;

datatype = { 'powermod' 'latr_wrt_resp' 'latr_wrt_choice'  'latr_wrt_stim' }; % use powermodulation or powerlateralization
idat = 2;

freqrange = 'full'; % low high or full
% freqrange = 'lowhigh'; % low high or full
% freqrange = 'low'; % low high or full

corrtype = 'Pearson';

statfolder = fullfile(basepath, 'projectdata/MEG2afc/freq/stats_5D', testdv{itest}, datatype{idat}, freqrange, corrtype, 'stat');
cd(statfolder)
outpath = fullfile(statfolder, 'plots');
mkdir(outpath)

freqfolder = fullfile(basepath, 'projectdata/MEG2afc/freq/stats_5D', testdv{itest}, datatype{idat}, freqrange, corrtype, 'freq');

idiffs = 3 %1:2 %3:4
idrugs = 2


THR = 0.5;


allstat = {};
allfreq = {};
for idrug = idrugs
    for idiff = idiffs
        for itrig = 1:2
            load( fullfile(cd, sprintf('stat_idrug%d_idiff%d_itrig%d_itest%d_idat%d_regr.mat', idrug, idiff, itrig, itest, idat) ));
            load( fullfile(freqfolder, sprintf('freq_idrug%d_idiff%d_itrig%d_itest%d_idat%d_regr.mat', idrug, idiff, itrig, itest, idat) ));
            if itest == 1 % modulation
                stat.rho = stat.stat;
            end
            
            allstat{idrug,idiff,itrig} = stat;
            allfreq{idrug,idiff,itrig} = data;
        end
    end
end

faxis_all = allstat{idrug,idiff,itrig}.freq;
plotfaxis = {};
switch freqrange
    case 'full'
        %         plotfaxis{1} = faxis_all(1:34);
        %         plotfaxis{2} = faxis_all(35:end);
        %         ifreqs = [2 1];
        plotfaxis{1} = faxis_all;
        ifreqs = 1;
    case 'low'
        plotfaxis{1} = faxis_all;
        ifreqs = 1;
    case 'high'
        plotfaxis{2} = faxis_all;
        ifreqs = 2;
    case 'lowhigh'
        plotfaxis{1} = faxis_all;
        ifreqs = 1;
end

taxis = {};
taxis{1} = allstat{idrug,idiff,1}.time;
taxis{2} = allstat{idrug,idiff,2}.time;

% TODO fix plotting for mod
% 1) TFR's for clusters p < 0.025, integrating over chans (of cluster)
% 2) topo, integrating over freq and time of cluster
% 3) time course of correlation values, integrating over chans and freq
% 4) power spectrum of correlation, integrating over chans and time

% plot integrated TFR's for stim and resp low and high freq
close all


% Condition labels
%--------------------------------------------------------------------------
pharm_conds = {'atomox' 'placebo' 'drugscomb' 'atomox - placebo'};
motor_conds = {'pressipsi' 'presscontra' 'regimecomb' 'pressipsi - presscontra'};
diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};
pupilconds = {'pupilhigh' 'pupillow' 'pupilcomb' 'pupilhi-lo'};
stim_conds = {'stimleft' 'stimright' 'allstim' 'stimleft-right'};
choice_conds = {'choiceleft' 'choiceright' 'allchoice' 'choiceleft-right'};
rt_conds = {'slow' 'medium' 'fast' 'allrts'};
correct_conds = {'correct' 'error' 'corr+err'};
sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
trigger_leg = {'stim', 'resp'};
freq_leg = {'low', 'full'};
sign = {'Pos', 'Neg'}; % neg or pos clusters

% Plot Hipp style figure for each cluster: TFR and topo

% close all
MEG2afc_sensorselection

trap = 1;
foi = [faxis_all(1) faxis_all(end)];
corrfactor = {'driftrate'};
icorr = 1;
SAV = 1;
XTICKS = -2:0.1:2;
% YTICKS = {0:5:150, 40:20:150};
YTICKS = {0:10:150, 40:20:150};

% for low and high freq plots
% nrows = 3;
% subplotinds = [5, 3];

% for full freq
nrows = 2;
subplotinds = [3, 1];

ZLIM = [-0.1 0.1]; % not really needed

for idrug = idrugs
    for idiff = idiffs
        for itrig = 1:2 %1:2
            plotstat =  allstat{idrug,idiff,itrig};
            
            for isign = 1:2
                if isign == 1
                    clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters));
                    cluslabels = plotstat.posclusterslabelmat;
                    pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
                else
                    clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
                    cluslabels = plotstat.negclusterslabelmat;
                    pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
                end
                
                if isempty(clus); continue; end
                
                for iclus = clus
                    figh = figure;    iplot=0; hold on
                    %                     set(gcf, 'Position', [0 -200 375*3 210*4])
                    set(gcf, 'Position', [100 150 1000 750])
                    chansel = any(any(cluslabels == iclus, 2),3); % get chans before high and low are separated
                    
                    for ifreq = ifreqs % high, then low freq
                        
                        plotstat =  allstat{idrug,idiff,itrig};
                        frind_log = ismember(plotstat.freq, plotfaxis{ifreq});
                        plotstat.prob = plotstat.prob(:,frind_log,:);
                        plotstat.posclusterslabelmat = plotstat.posclusterslabelmat(:,frind_log,:);
                        plotstat.negclusterslabelmat = plotstat.negclusterslabelmat(:,frind_log,:);
                        plotstat.freq = plotfaxis{ifreq};
                        
                        if isign == 1
                            clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters));
                            cluslabels = plotstat.posclusterslabelmat;
                            pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
                        else
                            clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
                            cluslabels = plotstat.negclusterslabelmat;
                            pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
                        end
                        
                        % for TFR
                        rtemp = plotstat.rho(:,frind_log,:);
                        %                         chansel = any(any(cluslabels == iclus, 2),3);
                        if trap
                            %                             rtemp(cluslabels~=clus(iclus)) = 0;
                            %                             rtemp = squeeze(trapz(rtemp,1));%./numel(find(rtemp~=0));
                            rtemp = squeeze(trapz(rtemp(chansel,:,:),1));%./numel(find(rtemp~=0));
                            scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                        else
                            rtemp(cluslabels~=clus(iclus)) = NaN;
                            rtemp = squeeze(nanmean(rtemp,1));%./numel(find(rtemp~=0));
                            rtemp(isnan(rtemp))=0;
                            scale = [-0.05 0.05];
                        end
                        cmap = cbrewer('div', 'RdBu',256);
                        cmap = flipud(cmap);
                        subplot(nrows,2, subplotinds(ifreq)); colormap(cmap)
                        %                                         rtemp = smooth2a(rtemp,1);
                        if isign == 1
                            thrcluster = squeeze(any(plotstat.posclusterslabelmat == iclus, 1));
                        else
                            thrcluster = squeeze(any(plotstat.negclusterslabelmat == iclus, 1));
                        end
                        %                         ft_plot_matrix(plotstat.time, plotstat.freq, rtemp, ... 'clim', scale);
                        %                             'highlight', thrcluster, 'highlightstyle', 'saturation', 'clim', scale)
                        
                        rtemp = showstatsTFR(rtemp, thrcluster, 1:length(plotfaxis{ifreq}), scale, 1);
                        
                        imagesc(taxis{itrig}, plotfaxis{ifreq}, rtemp, ZLIM);
                        title(sprintf('scale [%g %g]', scale))
                        hold on
                        %                         colorbar
                        contour(taxis{itrig}, plotfaxis{ifreq},  thrcluster * (iclus*2), 1, 'ShowText', 'On');
                        
                        set(gca,'YDir','normal', 'XTick', XTICKS, 'YTick', YTICKS{ifreq});
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
                        
                        %                         if isign == 2
                        %                             rtemp_time = abs(rtemp_time);
                        %                             rtemp_freq = abs(rtemp_freq);
                        %                         end
                        
                        subplot(nrows,2,1);
                        area(plotstat.time, rtemp_time); hold on;
                        xlim([plotstat.time(1), plotstat.time(end)]);
                        %                         ylim([-(ceil(max(abs(rtemp_time)))).*1.2, ceil(max(abs(rtemp_time))).*1.2]);
                        if isign == 1
                            YRANGE = [0, ceil(max(abs(rtemp_time))).*1.2];
                        else
                            YRANGE = [floor(min(rtemp_time)).*1.2, 0];
                        end
                        
                        ylim(YRANGE);
                        plot([0,0],YRANGE,'k');
                        
                        %                             ylim([0, ceil(max(abs(rtemp_time))).*1.2]);
                        %                             plot([0,0],[0, ceil(max(abs(rtemp_time))).*1.2],'k');
                        
                        %                         xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
                        ylabel('Integrated correlation');
                        if itest == 2
                            title(sprintf('%s\n%s Corr %s %s %s-%d\n p=%1.4f', datatype{idat}, corrfactor{icorr}, pharm_conds{idrug}, diff_conds{idiff}, sign{isign}, iclus, pvalues(iclus)));
                        else
                            title(sprintf('%s %s %s %s-%d\n p=%1.4f', datatype{idat}, pharm_conds{idrug}, diff_conds{idiff}, sign{isign}, iclus, pvalues(iclus)));
                        end
                        set(gca, 'XTick', XTICKS)
                        
                        % integrated corr plot 2
                        %                         subplot(nrows,2,subplotinds(ifreq)+1)
                        %                         area(plotstat.freq, rtemp_freq); hold on;
                        %
                        %                         if isign == 1
                        %                             YRANGE = [0, ceil(max(abs(rtemp_freq))).*1.2];
                        %                         else
                        %                             YRANGE = [floor(min(rtemp_freq)).*1.2, 0];
                        %                         end
                        %                         ylim(YRANGE);
                        %                         xlim([plotstat.freq(1), plotstat.freq(end)]);
                        %                         set(gca,'XAxisLocation','top','YAxisLocation','right', 'XTick', YTICKS{ifreq})
                        %                         camroll(270);  camorbit(0,180)
                        %                         ylabel('Integrated correlation');
                        %
                        %                         % for topo
                        %                         plotstat =  allstat{idrug,idiff,itrig};
                        %                         if isign == 1
                        %                             clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters));
                        %                             cluslabels = plotstat.posclusterslabelmat;
                        %                             pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
                        %                         else
                        %                             clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
                        %                             cluslabels = plotstat.negclusterslabelmat;
                        %                             pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
                        %                         end
                        
                        subplot(nrows,2,subplotinds(ifreq)+1)
                        % plot scatter of correlation for cluster
                        if itest == 2
                            if isign == 1
                                clusind = plotstat.posclusterslabelmat == iclus;
                            else
                                clusind = plotstat.negclusterslabelmat == iclus;
                            end
                            plotfreq = allfreq{idrug,idiff,itrig};
                            cluspow = plotfreq.powspctrm(:,clusind);                        % binary ROI:
                            
                            weighted_ROI = 0;
                            if weighted_ROI % TODO fix this
                                cluscorr = stat.rho(clusind);
                                for isub=1:nsub
                                    cluspow(isub,:) = cluspow(isub,:) .* cluscorr';
                                end
                            end
                            cluspow = mean(cluspow,2);
                            
                            set(0, 'DefaultAxesFontSize', 14)
                            %figure;
                            scatter(cluspow, plotfreq.drifts', 150, 'k', 'filled')
                            [r,p] = corr(cluspow, plotfreq.drifts');
                            title(sprintf('r = %.2f, p = %.3g', r, p))
                            axis square; box on;
                            if p < 0.05
                                h = lsline;
                                if p < 0.01
                                    set(h, 'Linewidth', 3)
                                end
                            end
                            xlabel('MEG power modulation (psc)')
                            ylabel('Drift rate')
                            
                            %                         xlim([-0.11 0.1])
                            %                         ylim([-0.3 0.11])
                            
                            xL = xlim;
                            yL = ylim;
                            line([0 0], yL, 'LineStyle', '--', 'color', [0.5 0.5 0.5 ])
                            line(xL, [0 0], 'LineStyle', '--', 'color', [0.5 0.5 0.5 ])
                        end
                        
                        % plot topo
                        subplot(nrows,2,2); colormap(cmap);
                        cfg = [];
                        cfg.layout = 'CTF275.lay';
                        cfg.marker = 'off';
                        cfg.markersize = 3;
                        cfg.markersymbol = '.'; %'o'
                        cfg.comment = 'no';                        
                        cfg.shading = 'flat';
                        cfg.style = 'both'; %both 'imsat' straight
                        cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
                        cfg.highlight = 'on';
                        %cfg.highlightstyle = 'opacity';
                        if isign == 1
%                             cfg.highlightchannel = chlabel(any(any(plotstat.posclusterslabelmat == iclus, 2),3));
                            cfg.highlightchannel = find(any(any(plotstat.posclusterslabelmat == iclus, 2),3));
                        else
%                             cfg.highlightchannel = chlabel(any(any(plotstat.negclusterslabelmat == iclus, 2),3));
                            cfg.highlightchannel = find(any(any(plotstat.negclusterslabelmat == iclus, 2),3));
                        end
                        
%                         cfg.highlightchannel = sens.ind{2};
                        
                        cfg.highlightsymbol = '.';
                        cfg.parameter = 'powspctrm';
                        cfg.interactive = 'no';
                        % cfg.maskparameter = 'mask';
                                                                        
                        rtemp = plotstat.rho;
                        if trap
                            rtemp(cluslabels~=clus(iclus)) = 0;
                            rtemp = squeeze(trapz(rtemp,2));
                            rtemp = squeeze(trapz(rtemp,2));
                            
%                             cluslog = cluslabels == iclus;
%                             cluslog = squeeze(any(cluslog,1));
%                             rtemp = squeeze(trapz(rtemp(:,cluslog),2));%./numel(find(rtemp~=0));

                            cfg.zlim = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                            %                             cfg.zlim = [0, ceil(max(abs(rtemp(:))))];
                        else
                            rtemp(cluslabels~=clus(iclus)) = NaN;
                            rtemp = squeeze(nanmean(rtemp,2));
                            rtemp = squeeze(nanmean(rtemp,2));
                            rtemp(isnan(rtemp))=0;
                            cfg.zlim = [-0.25 0.25];
                        end
                        
                        freq = [];
                        freq.label = plotstat.label; % chlabel(LR_subtract_mat(:,2));
                        freq.dimord = 'chan';
                        freq.powspctrm = rtemp;
                        
%                         freq.mask = zeros(size(chlabel)); %+ 0.25 
%                         freq.mask(ismember(chlabel, cfg.highlightchannel)) = 1;
                        
                        ft_topoplotTFR(cfg, freq);
                        
                        colorbar
                    end
                    
                    if SAV && ~isempty(clus) && pvalues(iclus) < 1
                        pstr = num2str(pvalues(iclus));
                        pstr = pstr(3:6);
                        outfile = sprintf('%s_corr_p=%s_%s_%slocked_%s_%s_%sclus%d', datatype{idat}, pstr, corrtype, trigger_leg{itrig}, pharm_conds{idrug}, diff_conds{idiff}, sign{isign}, clus(iclus));
                        display(outfile)
                        %                         print('-dpdf', fullfile(outpath, outfile))
                        %                         export_fig( fullfile(outpath, outfile), '-pdf')
                        export_fig( fullfile(outpath, outfile), '-png')
                        cd(outpath)
                    end
                end
            end
        end
    end
end




%% plot integrated TFR's, handy if several clusters present

% close all
% trap = 1;
% foi = [faxis_all(1) faxis_all(end)];
% % corrfactor = {'driftrate'};
% % icorr = 1;
% SAV = 0;
% ZLIM = [-0.1 0.1];

% for idrug = idrugs
%     for idiff = idiffs
%         figure; iplot = 0;
%         set(gcf, 'Position', [0 -200 375*3 210*4])
%
%         for ifreq = ifreqs % high, then low freq
%             for itrig= 1:2 %1:2
%
%                 plotstat =  allstat{idrug,idiff,itrig};
%                 % select freq range
%                 frind_log = ismember(plotstat.freq, plotfaxis{ifreq});
%                 plotstat.prob = plotstat.prob(:,frind_log,:);
%                 plotstat.posclusterslabelmat = plotstat.posclusterslabelmat(:,frind_log,:);
%                 plotstat.freq = plotfaxis{ifreq};
%
%                 if isign == 1
%                     clus = find(arrayfun(@(x) x.prob<=THR, plotstat.posclusters));
%                     cluslabels = plotstat.posclusterslabelmat;
%                     pvalues = arrayfun(@(x) x.prob, plotstat.posclusters(clus));
%                 else
%                     clus = find(arrayfun(@(x) x.prob<=THR, plotstat.negclusters));
%                     cluslabels = plotstat.negclusterslabelmat;
%                     pvalues = arrayfun(@(x) x.prob, plotstat.negclusters(clus));
%                 end
%
%                 % for TFR
%                 rtemp = plotstat.rho(:,frind_log,:);
%                 if trap
%                     %                             rtemp(cluslabels~=clus(iclus)) = 0;
%                     rtemp = squeeze(trapz(rtemp,1));%./numel(find(rtemp~=0));
%                     scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
%                     %                         else
%                     %                             rtemp(cluslabels~=clus(iclus)) = NaN;
%                     %                             rtemp = squeeze(nanmean(rtemp,1));%./numel(find(rtemp~=0));
%                     %                             rtemp(isnan(rtemp))=0;
%                     %                             scale = [-0.05 0.05];
%                 end
%                 cmap = cbrewer('div', 'RdBu',256);
%                 cmap = flipud(cmap);
%                 iplot = iplot + 1;
%                 subplot(2,2, iplot); hold on
%                 colormap(cmap)
%                 %                                         rtemp = smooth2a(rtemp,1);
%
%                 %                         thrcluster = squeeze(any(plotstat.posclusterslabelmat == iclus, 1));
%                 thrcluster = squeeze(any(ismember(plotstat.posclusterslabelmat, clus), 1));
%                 %                         ft_plot_matrix(plotstat.time, plotstat.freq, rtemp, ... 'clim', scale);
%                 %                             'highlight', thrcluster, 'highlightstyle', 'saturation', 'clim', scale)
%
%                 rtemp = showstatsTFR(rtemp, thrcluster, 1:length(plotfaxis{ifreq}), scale, 1);
%
%                 imagesc(taxis{itrig}, plotfaxis{ifreq}, rtemp,ZLIM);
%                 %                         colorbar
%                 if ~isempty(clus)
%                     for ic = clus
%                         thrcluster = squeeze(any(plotstat.posclusterslabelmat == ic, 1));
%                         thrcluster = thrcluster * (ic*2);
%                         contour(taxis{itrig}, plotfaxis{ifreq}, thrcluster, 1, 'ShowText', 'On');
%                     end
%                 end
%
%                 set(gca,'YDir','normal');
%                 if ifreq == 2
%                     set(gca, 'YTick', 40:20:150)
%                 end
%                 set(gca, 'XTick', -1:0.25:1)
%                 plot([0,0],foi,'k',[0,0],foi,'k');
%                 ylabel('Frequency (Hz)');
%                 xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
%                 xlim([plotstat.time(1), plotstat.time(end)]);
%                 ylim([plotfaxis{ifreq}(1) plotfaxis{ifreq}(end)])
%                 if ifreq == 2
%                     title(sprintf('%s %s %s-locked, thr = %g',  pharm_conds{idrug}, diff_conds{idiff}, trigger_leg{itrig}, THR ))
%                 end
%
%                 %                     end
%             end
%         end
%         if SAV && ~isempty(clus)
%             outfile = sprintf( '%s%sTFRstats_%s_%sclusters_%s_%s', outpath, filesep, corrtype, sign{isign},  pharm_conds{idrug}, diff_conds{idiff});
%             display(outfile)
%             export_fig(outfile, '-png') %'-png', ,  '-depsc'  '-transparent'
%
%         end
%
%     end
% end
% cd(outpath)

