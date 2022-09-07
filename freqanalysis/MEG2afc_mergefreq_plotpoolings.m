% %%  plot TFR for prespecified CTF O and C poolings, creating mask field
% mask: [274×34×15 logical]
% TODO select chan, freq and time bins, set to 1 in mask
% call ft_clusterplot3D :)


SAV = 1;
close all
load colormap_jetlightgray.mat
subplotind = [2 1; 3 4];
clussign = {'pos', 'neg'};
poolingnames = {'C' 'O' }; % occipital and Central
freqoi = [12 25; 40 80];

ifig = 1; % counter for saving
% imod = 2; idrug = 4; idiff = 3;
for itrig = 1
  f = figure;   f.Position = [ 680          75         612        792 ]; % A4 formaat
  irow = 0;
  
  for ip = 2:-1:1
    choi = ft_channelselection(sprintf('M*%s*', poolingnames{ip}), megdat.freq(1).label);
    [~, chanind] = ismember(choi, megdat.freq(1).label);
    for idrug = 2  %[2, 4]
      for imod = 1 % [1 3]
        for isign = 1%:2
%           freq = megdat.freq(idrug, imod, itrig, ip , 3);% ip instead of ifreq to get correct freq range
          cfg=[];
          cfg.latency = [-0.2 0.8];
          freq = ft_selectdata(cfg,  megdat.freq(idrug, imod, itrig, ip , 3));
          freq.powspctrm = squeeze(mean(freq.powspctrm));
          freq.dimord = 'chan_freq_time';
          %make mask
          freq.mask = false(size(freq.powspctrm));
          toi = freq.time > 0.2 & freq.time < 0.8;
          %             foi = freq.freq > 40 & freq.freq < 80;
          foi = freq.freq > freqoi(ip,1) & freq.freq < freqoi(ip,2) ;
          freq.mask(chanind, foi, toi) = true;
          
          %               curstat = megdat.stat{imod, idrug, ifreq,itrig};
          %               cfg.titleTFR = sprintf('%s\nidrug%d', curstat.megtype, idrug);
          cfg=[];
          cfg.clus2plot = 1;
          cfg.parameter = 'powspctrm';
          cfg.integratetype = 'mean'; % mean or trapz
          cfg.colormap = cmap;
          cfg.subplotsize = [4 4];
          cfg.subplotind = irow*4 + subplotind(itrig,:);
          if length(curstat.label) == 131 % only 131 sensors for latr
            cfg.layout = 'CTF275_helmet_latr.mat';
          else
            cfg.layout = 'CTF275_helmet.mat';
          end
          
          plotsuccess = ft_clusterplot3D(cfg, freq);
          
          irow = irow+1;
          
          irow = irow+1;
        end
      end
    end
    cd(megdat.PREOUT)
  end
  if SAV
    saveas(gcf, fullfile(megdat.PREOUT, sprintf('poolings%d_%s_idrug%d_stat.pdf', cfg.clus2plot,  curstat.megtype, idrug )))
  end
%   f = figure;      f.Position = [ 680          75         612        792 ]; % A4 formaat
%   ifig = ifig+1;
%   irow = 0;
end





