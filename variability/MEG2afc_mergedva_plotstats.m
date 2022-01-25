function MEG2afc_mergedva_plotstats(megdat)

%% plot drug and plac dva
close all
f = figure; iplot = 0;
for idva = 1:2
  for ichan = 1:2%:2
    iplot=iplot+1;
    subplot(2,2,iplot); hold on
    %   for idrug = 3%:2
    %     dvadat = megdat.dvalock(idrug,3);
    %     plot(dvadat.avg())
    
    cfg = [];
    cfg.channel = ichan;
    cfg.title = megdat.dvalock(1,idva,3).label{ichan}(1:10);
    cfg.fontsize = 12;
    ft_singleplotER(cfg, megdat.dvalock(1,idva,3), megdat.dvalock(2,idva,3))
    legend drug plac
    legend boxoff
    
  end
end

%% plot drug and plac dva error bars = 1;
SAV=1;
close all
f = figure;
linecols = 'rb';
ylab = {'within dva', 'within norm';'across dva', 'across norm' };
f = figure; iplot = 0;
for idva = 1:2
  for ichan = 1:2%:2
    iplot=iplot+1;
    subplot(2,2,iplot); hold on
    %   for idrug = 3%:2
    %     dvadat = megdat.dvalock(idrug,3);
    %     plot(dvadat.avg())
    clear H;
    for idrug = 1:2
      x = megdat.dvalock(idrug,idva,3).time;
      y = squeeze(mean(megdat.dvalock(idrug,idva,3).avg(:,ichan,:)));
      errBar = squeeze(std(megdat.dvalock(idrug,idva,3).avg(:,ichan,:)))/ sqrt(16);
      lineProps = linecols(idrug);
      H(idrug)=shadedErrorBar(x,y,errBar,lineProps,1);
      ax=gca; plot([0,0], ax.YLim, 'k')
    end
    legend([H.mainLine], {'drug' 'plac'})
    legend boxoff
    axis tight
    xlabel('Time from stim onset (s)')
    ylabel(ylab{idva,ichan})
    title(ylab{idva,ichan})
    
    plot_sig_bar(x, megdat.stat{idrug, idva}.mask) % vpos, height, Color
    
  end
end

if SAV 
  mkdir(megdat.PREOUT)
  saveas(gcf, fullfile(megdat.PREOUT, 'withindva_drugvsplac.png'))
  cd(megdat.PREOUT)
end


%%
% stat = megdat.stat;
% 
% %% multiplots
% close all
% ifreq = 1;
% itrig = 2;
% 
% load colormap_jetlightgray.mat
% cfg=[];
% cfg.layout = 'CTF275_helmet.mat';
% %   cfg.layout = 'CTF275.lay';
% %   cfg.baseline = [-0.25 0];
% %   cfg.baselinetype = 'relchange';
% cfg.colorbar = 'yes';
% cfg.colormap = cmap;
% cfg.zlim = 'maxabs';
% cfg.hotkeys = 'yes';
% 
% cfg.maskparameter = 'mask';
% cfg.maskalpha = 0.25;
% stat{ifreq,itrig}.mask = stat{ifreq,itrig}.posclusterslabelmat == 1; %  stat{ifreq,itrig}
% 
% % ctr = ctr + 1;
% % subplot(2,2,ctr)
% %   ft_multiplotTFR(cfg, freq(3,2,1,1)) % idrug, itrig, ifreq, idiff
% f = figure; 
% f.Position = [   680   444   781   654];
% ft_multiplotTFR(cfg, stat{ifreq,itrig}) % idrug, itrig, ifreq, idiff
% %   ft_multiplotTFR(cfg, freqall(3,1,2,1)) % idrug, itrig, ifreq, idiff
% % title(stat{ifreq,itrig}.posclusters(1).prob)
% title(sprintf('Drug - placebo, cluster p = %1.3f', stat{ifreq,itrig}.posclusters(1).prob))
% 
% %%
% 
% % ft_topoplotTFR(cfg, stat{ifreq,itrig}) % idrug, itrig, ifreq, idiff
% 
% 
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
% %   % log
%   pl = plot( freqsel(1).freq(1:end), log(freqsel(1).powspctrm(isub,:)), 'LineWidth', 2 );
%   pl = plot( freqsel(2).freq(1:end), log(freqsel(2).powspctrm(isub,:)), 'LineWidth', 2 );
%   
% %   pl = plot( freqsel(1).freq(1:end), freqsel(1).powspctrm(isub,:), 'LineWidth', 2 );
% %   pl = plot( freqsel(2).freq(1:end), freqsel(2).powspctrm(isub,:), 'LineWidth', 2 );
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
