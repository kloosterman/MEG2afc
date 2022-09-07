%%  plot 3D integrated cluster

SAV = 0;
close all
load colormap_jetlightgray.mat
subplotind = [2 1; 3 4];
clussign = {'pos', 'neg'};

ifig = 1; % counter for saving
% imod = 2; idrug = 4; idiff = 3;
for idrug =4  %[2, 4]
  for imod = 1 %[1 3]     
    cfg=[];
    cfg.clus2plot = 1; 
    cfg.parameter = 'powspctrm';
    cfg.integratetype = 'trapz'; % mean or trapz
    cfg.colormap = cmap;
    cfg.subplotsize = [4 4];    
    for isign = 1%:2
      f = figure;   f.Position = [ 680          75         612        792 ]; % A4 formaat
      irow = 0;
      for ifreq = 2:-1:1 % 1:2 %
        cfg.clussign = clussign{isign};
        
        % check if stim or resp locked is significant, if one, plot both
        pval.stim = megdat.stat{imod, idrug, ifreq,1}.([cfg.clussign 'clusters'])(cfg.clus2plot).prob;
        pval.resp = megdat.stat{imod, idrug, ifreq,2}.([cfg.clussign 'clusters'])(cfg.clus2plot).prob;
        disp(pval)
        if pval.stim && pval.resp > 0.9
          disp('Stim and Resp-locked not significant')
          continue
        end
        
        for itrig = 1:2
          curstat = megdat.stat{imod, idrug, ifreq,itrig};
          cfg.titleTFR = sprintf('%s\nidrug%d', curstat.megtype, idrug);
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
      
      if irow == 2 %|| ibehav == size(megdat.stat,1) % only 4 fit, starting from 0
        if SAV
          saveas(gcf, fullfile(megdat.PREOUT, sprintf('%sclus%d_%s_idrug%d_stat.pdf', cfg.clussign, cfg.clus2plot,  curstat.megtype, idrug )))
        end
        f = figure;      f.Position = [ 680          75         612        792 ]; % A4 formaat
        ifig = ifig+1;
        irow = 0;
      else
        irow = irow+1;
      end
    end
  end
end
cd(megdat.PREOUT)


%% plot multiplot
% function MEG2afc_mergefreq_plotstats(megdatraw)

stat = megdat.stat;

% multiplots
close all
idrug = 4; idiff = 3; imod = 3; itrig = 1; ifreq = 2;

load colormap_jetlightgray.mat
cfg=[];
% cfg.layout = 'CTF275_helmet.mat';
cfg.layout = 'CTF275_helmet_latr.mat';
%   cfg.layout = 'CTF275.lay';
%   cfg.baseline = [-0.25 0];
%   cfg.baselinetype = 'relchange';
cfg.colorbar = 'yes';
cfg.colormap = cmap;
if ifreq == 2
  cfg.ylim = [35 100];
end
cfg.zlim = 'maxabs';
% cfg.zlim = [-20 20];
cfg.hotkeys = 'yes';
if itrig == 1; cfg.xlim = [-0.2 0.6]; else cfg.xlim = [-0.8 0.2]; end

cfg.maskparameter = 'mask';
cfg.maskalpha = 0.25;
stat{idrug, imod, itrig, ifreq, idiff}.mask = stat{idrug, imod, itrig, ifreq, idiff}.posclusterslabelmat == 1; %  stat{idrug, imod, itrig, ifreq, idiff}
% stat{idrug, imod, itrig, ifreq, idiff}.mask = stat{idrug, imod, itrig, ifreq, idiff}.negclusterslabelmat == 1; %  stat{idrug, imod, itrig, ifreq, idiff}

% stat{idrug, imod, itrig, ifreq, idiff}.mask = stat{idrug, imod, itrig, ifreq, idiff}.prob < 0.05;

% ctr = ctr + 1;
% subplot(2,2,ctr)
%   ft_multiplotTFR(cfg, freq(3,2,1,1)) % idrug, itrig, ifreq, idiff
f = figure; 
f.Position = [   680   444   781   654];
ft_multiplotTFR(cfg, stat{idrug, imod, itrig, ifreq, idiff}) % idrug, itrig, ifreq, idiff
%   ft_multiplotTFR(cfg, freqall(3,1,2,1)) % idrug, itrig, ifreq, idiff
% title(stat{idrug, imod, itrig, ifreq, idiff}.posclusters(1).prob)
title(sprintf('Drug - placebo, cluster p = %1.3f', stat{idrug, imod, itrig, ifreq, idiff}.posclusters(1).prob))

%%  plot 3D integrated cluster

SAV = 1;
close all
load colormap_jetlightgray.mat
subplotind = [2 1; 3 4];
clussign = {'pos', 'neg'};

for idrug = 4 1:4
  for idiff = 3 3:4
    for imod = 1 1:3     
      cfg=[];
      if imod == 1
        cfg.layout = 'CTF275_helmet.mat';
      else
        cfg.layout = 'CTF275_helmet_latr.mat';
      end
      cfg.clus2plot = 1; % to do report pvals
      cfg.integratetype = 'trapz'; % mean or trapz
      cfg.colormap = cmap;
      cfg.subplotsize = [4 4];
      f = figure;
      f.Position = [ 680          75         612        792 ]; % A4 formaat
      
      irow = 0;
      for ifreq = 2:-1:1
        for isign = 1:2
          cfg.clussign = clussign{isign};
          for itrig = 1:2
            cfg.subplotind = irow*4 + subplotind(itrig,:);
            ft_clusterplot3D(cfg, megdat.stat{idrug, imod, itrig, ifreq, idiff})
          end
          irow = irow+1;
        end
      end
      if SAV
        saveas(gcf, fullfile(megdat.PREOUT, sprintf('mod%d_drug%d_diff%d_clus%d.pdf', imod, idrug, idiff, cfg.clus2plot)))
      end
    end
  end
end
cd(megdat.PREOUT)

%% Plot single subjects baseline spectrum
% close all
SAV = 1;

ifreq = 2;
itrig = 1;

cfg=[];

cfg.channel = 'M*O*'; % 
% cfg.channel = {'MRF25', 'MRT21'} % right frontal stuff
% cfg.channel = {'MLO11', 'MLO12', 'MLO13', 'MLO14', 'MLO21', 'MLO22', 'MLO23', 'MLO24', 'MLO31', 'MLO32', 'MLO33', 'MLO34', 'MLO44', 'MLP31', 'MLP41', 'MLP42', 'MLP51', 'MLP52', 'MLP53', 'MLP54', 'MLT16', 'MLT26', 'MLT27', 'MLT47', 'MRO11', 'MRO12', 'MRO21', 'MRO22', 'MRO23', 'MRO31', 'MRO32', 'MRP31', 'MRP41', 'MRP42', 'MRP51', 'MRP52', 'MRP53', 'MZO01', 'MZO02', 'MZP01'};
% cfg.channel = {'MLF24', 'MLF25', 'MRF25' , 'MRT21', 'MRT32'};
cfg.channel = {'MLO23'}; % left occ chan with high alpha in dr-pl
cfg.latency = [-0.5 0];
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
stat{1}.dimord = 'subj_chan_freq_time';
freqsel = ft_selectdata(cfg, stat{idrug, imod, itrig, ifreq, idiff})

% figure; plot( freqsel.freq, squeeze(freqsel.powspctrm)' )

datsel = squeeze(freqsel.powspctrm(:,:, 1:end))';
nsub = size(datsel,2);
% datsel = squeeze(freqsel.powspctrm([1:10, 12:end],:, 1:end))';
% datsel = squeeze(freqsel.powspctrm([3:10, 12:end],:, 1:end))';
% datsel = squeeze(freqsel.powspctrm(11,:, 1:end))';
figure; hold on
cmap = jet(nsub);
for isub = [1, 3:nsub]
  pl = plot( freqsel.freq(1:end), datsel(:,isub), 'LineWidth', 3, 'Color', cmap(isub,:) );
end
legend(cellstr(num2str([1:nsub]')))
title(['chans:' cfg.channel])
% set(gca, 'Xtick', 0:2:35); 
grid on
ylabel('PSD')
xlabel('freq (Hz)')

if SAV  
  saveas(gcf, fullfile(megdatraw.PREOUT, 'singlesub_drug-plac.png'))
  cd(megdatraw.PREOUT)
end

%% plot subj drug plac 1 by 1 subplots
SAV = 1
close all
ifreq = 1; 
cfg=[];
% cfg.latency = [-0.5 0];
cfg.latency = [ 0.25 0.75];
cfg.channel = 'M*O*'; % 
% cfg.channel = 'all'; % 
% cfg.channel = {'MRF25', 'MRT21'} % right frontal stuff
% cfg.channel = 'MRF25' % right frontal stuff

cfg.avgovertime = 'yes';
cfg.avgoverchan = 'yes';
clear freqsel
freqsel(1) = ft_selectdata(cfg, megdatraw.freq(1,1,ifreq,3)); % idrug, itrig, ifreq, idiff
freqsel(2) = ft_selectdata(cfg, megdatraw.freq(2,1,ifreq,3)); % idrug, itrig, ifreq, idiff

f = figure;
f.Position = [         357          53        1307        1052];
for isub = 1:nsub
  subplot(4,5,isub); hold on
% %   % log
%   pl = plot( freqsel(1).freq(1:end), log(freqsel(1).powspctrm(isub,:)), 'LineWidth', 2 );
%   pl = plot( freqsel(2).freq(1:end), log(freqsel(2).powspctrm(isub,:)), 'LineWidth', 2 );
  
  pl = plot( freqsel(1).freq(1:end), freqsel(1).powspctrm(isub,:), 'LineWidth', 2 );
  pl = plot( freqsel(2).freq(1:end), freqsel(2).powspctrm(isub,:), 'LineWidth', 2 );
  ylabel('PSD')
  xlabel('freq (Hz)')
  title(sprintf('subj %d\nchan %s', megdatraw.SUBJ(isub), cfg.channel))
  xlim([0 freqsel(1).freq(end)])
  %   ylim([0 1e-28])
  %   if SUBJ(isub) == 14
  %     ylim([0 1e-26])
  %   end
end
legend({'drug', 'plac'});

if SAV
  saveas(gcf,'singlesub_drugplac.png') 
end

%% plot CTF275 lay
cd /Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/fieldtrip/template/layout
cfg = [];
cfg.layout = 'CTF275_helmet.mat';
layout = ft_prepare_layout(cfg);

figure
ft_plot_layout(layout);
h = title(cfg.layout);
set(h, 'Interpreter', 'none');
