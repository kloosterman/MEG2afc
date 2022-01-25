%% plot TFR's for selected poolings
% function plotTFR_2AFC(~)
% plots SOI-TFRs
% average of selected sensors for multiple subjects
% sorted according to stimulus condition
% % % % % addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
% % % % % addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
% % % % % ft_defaults
% cd('/Users/kloosterman/gridmaster2012/kloosterman')
% MEG2afc_setup_paths

% MEG2afc_setup_paths
% trigger = 'stim' % both are now loaded!!

% MEG2afc_load_respavg
%
% freqrange = 2;

%% freqstatistics on TFR:
% test motor and occ vs 0
% motor and occ drug vs placebo
nsub=length(respavg.SUBJ);

loadstat = 0;
outfile = 'poolstat';
if loadstat;
  disp('Loading poolstat . . .')
  load(fullfile(PREIN, 'respavg', outfile));
  return
end
freq=[];
% freq.dimord = 'chan_subj_freq_time';
freq.dimord = 'subj_chan_freq_time';
freq.label = {'custompooling'};

cfg = [];
% cfg.frequency        = YLIM(3,:);
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
% cfg.alpha            = 0.1;
cfg.numrandomization = 100;
cfg.neighbours       = []; %in case no channel data present
cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = 'ctf275_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);

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
for ilatr = 1%:4 % latr_leg = {'modulation', 'latr wrt resp', 'latr wrt choice', 'latr wrt stim'};
  for isoi = 1 %1%:2 %:2 % occ and motor
    for idiff = 3% 3:4 %1:4
      for idrug = 4 %3:4 %1:4 %1:4% 1:4
        for itrig =  1:2 %1:4 %1:4% 1:4
          for ifreq = 1
            %                     cfg.latency          = respavg.time{itrig};
            freq.time = respavg.time{itrig};
            freq.freq = respavg.freq{ifreq};
            
            %                         %average over SENSORS
            %                         freq.powspctrm = squeeze(mean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr), 2)); %tind{itrig}
            
            % in 3D
            freq.powspctrm = squeeze(respavg.dat(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr));
%             freq.powspctrm = squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr));
            freq.label = respavg.label;
            
            %                         freq.powspctrm = squeeze(respavg(:,isoi,:,tind{itrig}, idrug,3,idiff, itrig)); %DIMS:
            
            %                         freq.powspctrm = shiftdim(freq.powspctrm, -1);
            freqzero = freq; %create zero freq to test against
            freqzero.powspctrm = zeros(size(freq.powspctrm));
            poolstat{ilatr, isoi, idrug, idiff, itrig, ifreq} = ft_freqstatistics(cfg, freq, freqzero);
          end
        end
      end
    end
  end
end
% save(fullfile(PREIN, outfile), 'poolstat');

%% plotting: 2AFC TFR drug vs placebo: motor and visual cortex + stats
close all
SAV=1;
set(0,'DefaultAxesFontsize',10)

addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/plotting');
showstats = 1;

iband = 3;

TYP = '2afc';
for ilatr = 1%:4 % latr_leg = {'modulation', 'latr wrt resp', 'latr wrt choice', 'latr wrt stim'};
  
  figh = figure;    iplot=0; hold on
  set(gcf, 'Position', [200 200 1200 800])
  %     load('colormap170613.mat');
  load colormap_jetlightgray.mat
  colormap(cmap);  %cmap = get(gcf, 'Colormap')
  % colormap(jet(256));  %cmap = get(gcf, 'Colormap')
  %         for ipup= [3 1 2 4] %1:4
  for isoi = 1%:2 %:4 %1:3%  1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
    for idiff = 3 %:4 % 1:4
      for idrug = 4 %3:4 %[1:2,4] % 1:2 1:4 %1:4 %
        for ifreq = 1% 2:-1:1
          for itrig = 1:2
            ZLIM = [-0.05 0.05];
            %                     ZLIM = [-0.1 0.1];
            iplot = iplot+1;
            subplot(2,2,iplot); hold on
            if iplot > 12;  ZLIM = [-0.1 0.1];  end
            if idiff == 4;  ZLIM = [-0.05 0.05];  end
            if idrug == 4;  ZLIM = [-0.05 0.05];  end
            %                     if idrug == 4;  ZLIM = [-0.025 0.025];  end
            
            %                     % subj chan freq time ipharm, idiff, itrig, ifreq, ilatr
            %                     freq.powspctrm = squeeze(nanmean(respavg(:,sens.ind{isoi}, frind{3}, tind{itrig} ,idrug, idiff, itrig, 3, ilatr))); %average over subj
            %                     freq.powspctrm = shiftdim(freq.powspctrm, -1);
            %                     cfg_p = [];
            %                     cfg_p.maskparameter = poolstat{ilatr, isoi, idrug, idiff, itrig}.mask;
            %                     cfg_p.maskstyle = 'saturation';
            % %                     cfg_p.zlim = ZLIM;
            %                     ft_singleplotTFR(cfg, freq);
            
            % old
            % subj chan freq time ipharm, idiff, itrig, ifreq, ilatr
            dum = squeeze(nanmean(respavg.dat(:,respavg.sens.ind{isoi}, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr))); %average over subj tind{itrig}
            dum = squeeze(nanmean(dum)); % average over sensors, nanmean for lateralizations (half sensors are nan)
            %                     if showstats
            %                         dum = showstatsTFR(dum, poolstat{ilatr, isoi, idrug, idiff, itrig}.prob, frind{3}, ZLIM, showstats);
            %                     end
            %                     imagesc(taxis{itrig}, faxis{3}, dum, ZLIM);
            
            if showstats
              mask = double(squeeze(poolstat{ilatr, isoi, idrug, idiff, itrig, ifreq}.mask));
              mask(mask==0) = 0.25;
              ft_plot_matrix(respavg.time{itrig}, respavg.freq{ifreq}, dum, 'clim', ZLIM, 'box', 'no', ...
                'highlight', mask, 'highlightstyle', 'opacity'); % opacity
            else
              ft_plot_matrix(respavg.time{itrig}, respavg.freq{ifreq}, dum, 'clim', ZLIM, 'box', 'no' ); % opacity
            end
            
            
            title({ sprintf('%s %s ', respavg.latr_leg{ilatr}, respavg.sens.leg{isoi})
              sprintf('%s %s [%g]', respavg.pharm_conds{idrug}, respavg.diff_conds{idiff}, ZLIM(2)) }); % 'FontWeight','bold'
            hold on
            yaxis = [respavg.freq{ifreq}(1), respavg.freq{ifreq}(end)];
            plot([0,0],yaxis,'k',[0,0],yaxis,'k');
            
            %                     xlim(XLIM(itrig,:))
            xlim([respavg.time{itrig}(1) respavg.time{itrig}(end)])
            ylim(yaxis)
            
            %                     ylim([FREQLO(3) 35])
            YTICKS = [0:10:200];
            set(gca,'Box','off','XTick',-2:0.1:2,...    [-0.5,0,0.5,1,1.5,2,2.5]
              'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
              'TickDir','out', 'FontSize', 12);
            if itrig == 1
              ylabel('Frequency (Hz)');
            else
              %                         set(gca, 'YTickLabel', []);
            end
            xlabel(sprintf('Time from %s (s)', respavg.trigger_leg{itrig}));
            %                                         h=colorbar;
            box on
            
          end
        end
      end
    end
  end
  
  if SAV
    outpath = fullfile(respavg.PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%sTFRpooling_%sfreq_%s_%s', outpath, filesep, TYP,  respavg.sens.leg{isoi}, respavg.latr_leg{ilatr} ); % motor_conds{imotor}, pupilconds{ipup}
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
    cd(outpath)
  end;
end

%% multiplot

ilatr = 1; idrug=4; idiff=3; itrig=2; ifreq=1;

freq=[];
freq.label = respavg.label;
freq.time = respavg.time{itrig};
freq.freq = respavg.freq{ifreq};
freq.powspctrm = squeeze(respavg.dat(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr));
% freq.powspctrm = squeeze(respavg.pow(:,:, 1:length(respavg.freq{ifreq}), 1:length(respavg.time{itrig}) ,idrug, idiff, itrig, ifreq, ilatr));
freq.powspctrm = freq.powspctrm *1e21;
freq.dimord = 'subj_chan_freq_time';


freq.mask = poolstat{ilatr, isoi, idrug, idiff, itrig, ifreq}.mask;
freq.powspctrm(~freq.mask) = nan;
% freq.mask = poolstat{ ilatr, isoi, idrug, idiff, itrig, ifreq }.posclusterslabelmat == 1 ; % | ...
cfg=[];
%         %                 freq.mask = poolstat{ itrig, icond, istim, iresp }.negclusterslabelmat == 1 ; % | ...
%         %                                      freq.mask = poolstat{ itrig, icond, istim, iresp }.posclusterslabelmat == 1;
%         % % % %                 % %                 freq.mask = poolstat{ itrig, icond, istim, iresp }.negclusterslabelmat == 3;
%         % % % % %                 freq.mask = poolstat{ itrig, icond, istim, iresp }.posclusterslabelmat == 1;
%         %
cfg.maskparameter = 'mask';
cfg.maskalpha = 0; 0.35;


cfg.layout = 'CTF275.lay';  %CTF275_helmet
cfg.layout = ft_prepare_layout(cfg);
% cfg.layout.width(:) = 0.075;
% cfg.layout.height(:) = 0.075;

% cfg.xlim = [-0.75 2];
%                 cfg.zlim = [-0.1 0.1];
cfg.zlim = 'maxabs';
%         cfg.zlim = 'maxmin';
% cfg.zlim = [-1 1];
% load( 'colormap_jetlightgray.mat')
cfg.colormap = cmap;
cfg.hotkeys = 'yes';
cfg.fontsize = 18;
cfg.colorbar = 'yes';

f = figure;
f.Position = [ 680   678   1200   1200];

ft_multiplotTFR(cfg, freq);
