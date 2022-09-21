% function MEG2afc_mergefreq_corrstats_plot(megdat)
%%  plot 3D integrated cluster

SAV = 1;
close all
load colormap_jetlightgray.mat
subplotind = [2 1; 3 4];
clussign = {'pos', 'neg'};

ifig = 1; % counter for saving
f = figure;   f.Position = [ 680          75         612        792 ]; % A4 formaat
irow = 0;
% imod = 2; idrug = 4; idiff = 3;
for imod = 1 %1:4
  for ibehav = 1:size(megdat.corrstat,2)
    
    cfg=[];
    cfg.clus2plot = 1; 
    cfg.parameter = 'rho';
    cfg.integratetype = 'trapz'; % mean or trapz
    cfg.colormap = cmap;
    cfg.subplotsize = [4 4];
    
    for isign = 1%:2
      for ifreq = 1 %2:-1:1 % 1:2 %
        cfg.clussign = clussign{isign};
        
        % check if stim or resp locked is significant, if one, plot both
        pval.stim = megdat.corrstat{imod, ibehav,ifreq,1}.([cfg.clussign 'clusters'])(cfg.clus2plot).prob;
        pval.resp = megdat.corrstat{imod, ibehav,ifreq,2}.([cfg.clussign 'clusters'])(cfg.clus2plot).prob;
        disp(pval)
        if pval.stim && pval.resp > 0.9
          disp('Stim and Resp-locked not significant')
          continue
        end
        
        for itrig = 1:2
          curstat = megdat.corrstat{imod, ibehav,ifreq,itrig};          
          cfg.titleTFR = sprintf('%s\n',curstat.behavname{:}, curstat.megtype);;
          cfg.subplotind = irow*4 + subplotind(itrig,:);
          if length(curstat.label) == 131 % only 131 sensors for latr
            cfg.layout = 'CTF275_helmet_latr.mat';
          else
            cfg.layout = 'CTF275_helmet.mat';
          end
          
          plotsuccess = ft_clusterplot3D(cfg, curstat);
          
        end
        irow = irow+1;
      end
        
%       if plotsuccess
%         if irow == 3 || ibehav == size(megdat.corrstat,1) % only 4 fit, starting from 0
          if SAV
            %               saveas(gcf, fullfile(megdat.PREOUT, sprintf('corr_%svs%s.pdf',  megdat.corrstat{ibehav,1,1}.megtype,  [megdat.corrstat{ibehav,1,1}.behavname{:}] )))
            saveas(gcf, fullfile(megdat.PREOUT, sprintf('clus%d_%s_%s_corr.pdf', ifig,  curstat.megtype, curstat.behavname{:} )))
          end
          f = figure;      f.Position = [ 680          75         612        792 ]; % A4 formaat
          ifig = ifig+1;
          irow = 0;
%         else
%           irow = irow+1;
%         end
        %       end
        plotscatter = 1;
        if plotscatter
          bcorr = ft_findcfg(curstat.cfg, 'design')';
          acorr = curstat.scatterdat;
          f = figure; f.Position = [744   950   147   100];
          scatter(acorr, bcorr(:,1), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 30);
          box on; axis square;  axis tight
%           set(gca, 'XLim', [-3.3840 3.3840])
%           set(gca, 'YLim', [-0.36 0.36])
          title(sprintf('r = %1.2f, rho = %1.2f', partialcorr(acorr, bcorr(:,1), bcorr(:,2:end), 'type', 'Pearson'), partialcorr(acorr, bcorr(:,1), bcorr(:,2:end), 'type', 'Spearman' )))
          lsline
          saveas(gcf, fullfile(megdat.PREOUT, sprintf('clus%d_vs_%s_%s_scatter.pdf', ifig,  curstat.megtype, curstat.behavname{:} ))) %
        end
    end
  end
end
cd(megdat.PREOUT)



%% plot scatter of cluster OLD
for ibehav = 1:size(megdat.corrstat,1)
  f = figure; iplot=0;
  f.Position = [ 680          75         612        792 ]; % A4 formaat
  titstr = sprintf('%s\n', megdat.corrstat{ibehav,1,1}.behavname{:});
  
  irow = 0;
  for ifreq = 2:-1:1
    for isign = 1:2
      for itrig = 1:2
        curstat = megdat.corrstat{ibehav,ifreq,itrig};
        mask = curstat.([clussign{isign} 'clusterslabelmat']) == 1;
        seldat = mean(curstat.powspctrm_subj(:, mask),2);        
        design = ft_findcfg( curstat.cfg, 'design');

        iplot = iplot+1;
        subplot(4,2,iplot)
        scatter(design(:), seldat(:))
        title(titstr)
        text(design(:), seldat(:), num2cell(megdat.SUBJ(ismember(megdat.SUBJ, curstat.SUBJ_inc))));
      end
    end
  end
  if SAV
    saveas(gcf, fullfile(megdat.PREOUT, sprintf('corrscatter_megvs%s.pdf', titstr )))
  end
end
length(seldat)
cd(megdat.PREOUT)


% % drugdat = megdat.freq(1,imod,itrig,ifreq,3);
% % placdat = megdat.freq(2,imod,itrig,ifreq,3);
% % 
% % cfg2=[];
% % cfg2.parameter = 'powspctrm';
% % cfg2.operation = 'subtract';
% % freq = ft_math(cfg2, drugdat, placdat); % idrug, itrig, ifreq, idiff
% % seldat = mean(freq.powspctrm(:, curstat.mask),2);
% 
% mask = curstat.([clussign 'clusterslabelmat']) == 1;
% % seldat = mean(freq.powspctrm(:, mask),2);
% 
% design = ft_findcfg( curstat.cfg, 'design');
% % seldat = nanmean(curstat.powspctrm_subj(:, curstat.mask),2); 
% 
% figure; scatter(design(:), seldat(:))
% title(sprintf('r = %g', corr(design(:), seldat(:))))
% text(design(:), seldat(:), num2str(megdat.SUBJ'))
% % xlabel('atx-plac rep prob')
% % ylabel('atx-plac power')
xlabel('atx-plac behavior')
ylabel('atx-plac power')


%% ft_multiplotTFR
% stat = megdat.corrstat;

% multiplots
% close all
% % beta corr, p=0.02
% ifreq = 1; itrig = 2; imod = 3; 
% clussign = 'neg';

ifreq = 2; itrig = 1; imod = 1; ibehav=1; clussign = 'pos';
curstat = megdat.corrstat{imod, ibehav,ifreq,itrig};
% % % 
% % gamma corr with z
% ifreq = 2; itrig = 2; imod = 3; clussign = 'pos';

load colormap_jetlightgray.mat
cfg=[];
if imod == 1
  cfg.layout = 'CTF275_helmet.mat';
else
  cfg.layout = 'CTF275_helmet_latr.mat';
end
%   cfg.baseline = [-0.25 0];
%   cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.colormap = cmap;
cfg.zlim = 'maxabs';
cfg.hotkeys = 'yes';

% cfg.maskparameter = 'mask';
% cfg.maskalpha = 0.25;
curstat.mask = curstat.([clussign 'clusterslabelmat']) == 1; %  curstat
% curstat.mask = curstat.prob < 0.05;

cfg.parameter = 'rho';
% ctr = ctr + 1;
% subplot(2,2,ctr)
%   ft_multiplotTFR(cfg, freq(3,2,1,1)) % idrug, itrig, ifreq, idiff
f = figure;
f.Position = [   680   444   781   654];
ft_multiplotTFR(cfg, curstat) % idrug, itrig, ifreq, idiff
%   ft_multiplotTFR(cfg, freqall(3,1,2,1)) % idrug, itrig, ifreq, idiff
% title(curstat.posclusters(1).prob)
title(sprintf('corr atx - plac power vs ..., cluster p = %1.3f', curstat.([clussign 'clusters'])(1).prob))
% title(sprintf('corr atx - plac power vs rep prob, cluster p = %1.3f', curstat.negclusters(1).prob))


% %% Plot single subjects baseline spectrum
% % close all
% SAV = 1;
% 
% ifreq = 1;
% itrig = 1;
% 
% cfg=[];
% 
% cfg.channel = 'M*O*'; %
% % cfg.channel = {'MRF25', 'MRT21'} % right frontal stuff
% % cfg.channel = {'MLO11', 'MLO12', 'MLO13', 'MLO14', 'MLO21', 'MLO22', 'MLO23', 'MLO24', 'MLO31', 'MLO32', 'MLO33', 'MLO34', 'MLO44', 'MLP31', 'MLP41', 'MLP42', 'MLP51', 'MLP52', 'MLP53', 'MLP54', 'MLT16', 'MLT26', 'MLT27', 'MLT47', 'MRO11', 'MRO12', 'MRO21', 'MRO22', 'MRO23', 'MRO31', 'MRO32', 'MRP31', 'MRP41', 'MRP42', 'MRP51', 'MRP52', 'MRP53', 'MZO01', 'MZO02', 'MZP01'};
% % cfg.channel = {'MLF24', 'MLF25', 'MRF25' , 'MRT21', 'MRT32'};
% cfg.channel = {'MLO23'}; % left occ chan with high alpha in dr-pl
% cfg.latency = [-0.5 0];
% cfg.avgovertime = 'yes';
% cfg.avgoverchan = 'yes';
% stat{1}.dimord = 'subj_chan_freq_time';
% freqsel = ft_selectdata(cfg, stat{ifreq,itrig})
% 
% % figure; plot( freqsel.freq, squeeze(freqsel.powspctrm)' )
% 
% datsel = squeeze(freqsel.powspctrm(:,:, 1:end))';
% nsub = size(datsel,2);
% % datsel = squeeze(freqsel.powspctrm([1:10, 12:end],:, 1:end))';
% % datsel = squeeze(freqsel.powspctrm([3:10, 12:end],:, 1:end))';
% % datsel = squeeze(freqsel.powspctrm(11,:, 1:end))';
% figure; hold on
% cmap = jet(nsub);
% for isub = [1, 3:nsub]
%   pl = plot( freqsel.freq(1:end), datsel(:,isub), 'LineWidth', 3, 'Color', cmap(isub,:) );
% end
% legend(cellstr(num2str([1:nsub]')))
% title(['chans:' cfg.channel])
% set(gca, 'Xtick', 0:2:35); grid on
% ylabel('PSD')
% xlabel('freq (Hz)')
% 
% if SAV
%   saveas(gcf, fullfile(megdat.PREOUT, 'singlesub_drug-plac.png'))
%   cd(megdat.PREOUT)
% end
% 
% %% plot subj drug plac 1 by 1 subplots
% SAV = 1
% close all
% cfg=[];
% cfg.latency = [-0.5 0];
% % cfg.channel = 'M*O*'; %
% cfg.channel = {'MRF25', 'MRT21'} % right frontal stuff
% cfg.channel = 'MRF25' % right frontal stuff
% 
% cfg.avgovertime = 'yes';
% cfg.avgoverchan = 'yes';
% clear freqsel
% freqsel(1) = ft_selectdata(cfg, megdat.freq(1,1,1,4)); % idrug, itrig, ifreq, idiff
% freqsel(2) = ft_selectdata(cfg, megdat.freq(2,1,1,4)); % idrug, itrig, ifreq, idiff
% 
% f = figure;
% f.Position = [         357          53        1307        1052];
% for isub = 1:nsub
%   subplot(4,5,isub); hold on
%   %   % log
%   pl = plot( freqsel(1).freq(1:end), log(freqsel(1).powspctrm(isub,:)), 'LineWidth', 2 );
%   pl = plot( freqsel(2).freq(1:end), log(freqsel(2).powspctrm(isub,:)), 'LineWidth', 2 );
%   
%   %   pl = plot( freqsel(1).freq(1:end), freqsel(1).powspctrm(isub,:), 'LineWidth', 2 );
%   %   pl = plot( freqsel(2).freq(1:end), freqsel(2).powspctrm(isub,:), 'LineWidth', 2 );
%   ylabel('PSD')
%   xlabel('freq (Hz)')
%   title(sprintf('subj %d\nchan %s', megdat.SUBJ(isub), cfg.channel))
%   xlim([0 freqsel(1).freq(end)])
%   %   ylim([0 1e-28])
%   %   if SUBJ(isub) == 14
%   %     ylim([0 1e-26])
%   %   end
% end
% legend({'drug', 'plac'});
% 
% if SAV
%   saveas(gcf,'singlesub_drugplac.png')
% end
% 
% %% plot CTF275 lay
% cd /Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/layout
% cfg = [];
% cfg.layout = 'CTF275_helmet.mat';
% layout = ft_prepare_layout(cfg);
% 
% figure
% ft_plot_layout(layout);
% h = title(cfg.layout);
% set(h, 'Interpreter', 'none');
