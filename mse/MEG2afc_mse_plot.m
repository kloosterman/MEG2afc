% function critEEG_mse_plot(mseavg)



try isstruct(mseavg) % nargin == 0
catch
  mseavg = MEG2afc_mse_merge_runs();
end

%% timelockstatistics on chan X scales:
% put scales in time dim
nsub = length(mseavg.SUBJ);

cfg = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
% cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
% cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;
% prepare_neighbours determines what sensors may form clusters
cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/ctf275_neighb.mat';
% cfg0_neighb.elecfile  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/electrode/standard_1020.elc';
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

% cfg.frequency = [1 21]
cfg.latency = [-0.5 1]
cfg.minnbchan = 1

poolstat=cell( 2,3,4,2,2 );
ises = 4; % 4 = collapsed over session
for icond = 3%:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
  for istim = 1 %:2%:3
    for iresp = 1%:2
      for imse = 2 % 1:2
        % 5D stats:cc
        freq = [];
        freq.dimord = 'subj_chan_freq_time';
        freq.label = mseavg.label;
        freq.freq = mseavg.timescales;
        freq.time = mseavg.time;
        freq.powspctrm = squeeze( mseavg.dat(:, imse, :,:,:,  icond, istim, iresp) );
        %                 freq.powspctrm = shiftdim(freq.powspctrm, -1);
        
        freqzero = freq; %create zero freq to test against
        freqzero.powspctrm = zeros(size(freq.powspctrm));
        poolstat{ imse, icond, istim, iresp } = ft_freqstatistics(cfg, freq, freqzero);
        
      end
    end
  end
end

%% plot multiplot TFR for significant clusters

close all
for idrug = 1
  for imotor = 1
    for idiff = 1%:2
      for imse = 1 % raw or psc or r
        cfg = [];
        freq = [];
        %                 freq.mask = poolstat{ imse, icond, istim, iresp }.mask;
        %                 freq.mask = poolstat{ imse, icond, istim, iresp }.negclusterslabelmat == 2 ; % | ...
        % %                       poolstat{ imse, icond, istim, iresp }.posclusterslabelmat == 1;
        % %                 freq.mask = poolstat{ imse, icond, istim, iresp }.negclusterslabelmat == 3;
        % %                 freq.mask = poolstat{ imse, icond, istim, iresp }.posclusterslabelmat == 1;
        %                 cfg.maskparameter = 'mask';
        %                 if ~any(freq.mask)
        %                     continue;
        %                 end
        %             freq.time = respavg.time{itrig};
        %             freq.freq = respavg.freq{iband};
        %                 freq.powspctrm = squeeze( mseavg.datraw(:,4, imse, :,:,:,  icond, istim, iresp) );
        
        freq.powspctrm = squeeze(  mseavg.dat(:, imse, :,:,:, idrug, imotor, idiff) );
        
        freq.label = mseavg.label;
        freq.freq = mseavg.timescales;
        freq.time = mseavg.time;
        %                 freq.dimord = 'subj_chan_freq_time';
%         freq.dimord = 'subj_chan_freq_time';
        freq.dimord = 'chan_freq_time';
        
        %                 cfg.layout = 'CTF275.lay';  % biosemi64  elec1010
        cfg.layout = 'CTF275_helmet.mat';  % biosemi64  elec1010
        cfg.layout = ft_prepare_layout(cfg);
        %                 cfg.layout.width(:) = 0.075;
        %                 cfg.layout.height(:) = 0.075;
        
        %                 cfg.xlim = [-0.75 0.75];
        %                 cfg.zlim = [-0.1 0.1];
        cfg.zlim = 'maxabs';
        %                 cfg.zlim = 'maxmin';
        % cfg.zlim = [-1 1];
        load( 'colormap_jetlightgray.mat')
        cfg.colormap = cmap;
        %                 cfg.colormap = cmap(129:end,:);
        cfg.hotkeys = 'yes';
        cfg.fontsize = 18;
        cfg.colorbar = 'yes';
        
        f = figure;
        f.Position = [ 680   678   1200   1200];
        
        ft_multiplotTFR(cfg, freq);
        
        title(sprintf('%s, %s', mseavg.mse_leg{imse}, mseavg.pharm_conds{idrug}))
      end
    end
  end
end

%% line plot multiplot 
%     dimord: 'chan_time'       defines how the numeric data should be interpreted
%        avg: [151x600 double]  the numeric data (in this example it contains the average values of the activity for 151 channels x 600 timepoints)
%      label: {151x1 cell}      the channel labels (e.g. 'MRC13')
%       time: [1x600 double]    the timepoints in seconds
%        var: [151x600 double]  the variance of the activity for 151 channels x 600 timepoints
%       grad: [1x1 struct]      information about the sensor array (for EEG-data it is called elec)
%        cfg: [1x1 struct]      configuration structure used by the invoking FieldTrip function

scoi = 1:10;
% scoi = 11:42;

close all
for idrug = 4
  for imotor = 3
    for idiff = 3%:2
      for imse = 1 % raw or psc or r
        cfg = [];
        timelock = [];       
        timelock.avg = squeeze(  mseavg.dat(:, imse, :,:,:,  idrug, imotor, idiff) );
        timelock.avg = squeeze(mean(timelock.avg(:,:,scoi,:),3));
        
        timelock.label = mseavg.label;        
        timelock.time = mseavg.time;
        %                 timelock.dimord = 'subj_chan_timelock_time';
        timelock.dimord = 'subj_chan_time';
        
        %                 cfg.layout = 'CTF275.lay';  % biosemi64  elec1010
        cfg.layout = 'CTF275_helmet.mat';  % biosemi64  elec1010
        cfg.layout = ft_prepare_layout(cfg);
        %                 cfg.layout.width(:) = 0.075;
        %                 cfg.layout.height(:) = 0.075;
        
        %                 cfg.xlim = [-0.75 0.75];
        %                 cfg.zlim = [-0.1 0.1];
        cfg.zlim = 'maxabs';
        %                 cfg.zlim = 'maxmin';
        % cfg.zlim = [-1 1];
        load( 'colormap_jetlightgray.mat')
        cfg.colormap = cmap;
        %                 cfg.colormap = cmap(129:end,:);
        cfg.hotkeys = 'yes';
        cfg.fontsize = 18;
        cfg.colorbar = 'yes';
        
        f = figure;
        f.Position = [ 680   678   1200   1200];
        
        ft_multiplotER(cfg, timelock);
        
        title(sprintf('%s, %s', mseavg.mse_leg{imse}, mseavg.pharm_conds{idrug}))
      end
    end
  end
end




%%
cfg2 = [];
cfg2.toilim = [-1 1.5];
data = ft_redefinetrial(cfg2, data)

cfg2 = [];
cfg2.vartrllength = 2;
timelock = ft_timelockanalysis(cfg2, data);

cfg2 = [];
%                 cfg.layout = 'CTF275.lay';  % biosemi64  elec1010
cfg.layout = 'CTF275_helmet.mat';  % biosemi64  elec1010
cfg2.hotkeys = 'yes';

cfg2.parameter = 'var' % var avg
cfg2.xlim =   [-1 1.5];
cfg2.zlim = 'maxabs';
cfg2.colorbar = 'yes';
ft_multiplotER(cfg2, timelock)


%% plot multiplot single mse struct
% drug = load('NK3_ses1_20141121_drug_ipsi_run1_stim_data.mat');
% plac = load('NK3_ses4_20141206_plac_ipsi_run1_stim_data.mat');
drug = load('NK3_ses2_20141125_drug_contra_run1_stim_data.mat');
plac = load('NK3_ses3_20141205_plac_contra_run1_stim_data.mat');

close all
cfg = [];
% freq.powspctrm = plac.mse.sampen;
% freq.powspctrm = mse.r;
freq.powspctrm = drug.mse.sampen - plac.mse.sampen;
% freq.powspctrm = drug.mse.r - plac.mse.r;

freq.label = mse.label;
freq.freq = mse.timescales;
freq.time = mse.time;
freq.dimord = 'chan_freq_time';

% cfg.layout = 'CTF275.lay';  % biosemi64  elec1010
cfg.layout = 'CTF275_helmet.mat';  % biosemi64  elec1010
% cfg.layout = ft_prepare_layout(cfg);
%                 cfg.layout.width(:) = 0.075;
%                 cfg.layout.height(:) = 0.075;

cfg.xlim = [-0.75 0.75];
cfg.zlim = 'maxabs';
% load( 'colormap_jetlightgray.mat')
cfg.colormap = cmap;
cfg.hotkeys = 'yes';
cfg.fontsize = 18;
cfg.colorbar = 'yes';
% cfg.baseline = [-0.75 -0.5];
% cfg.baseline = [-0.5 -0.25];
% cfg.baselinetype = 'relchange' ;

f = figure;
f.Position = [ 680   678   1200   1200];

ft_multiplotTFR(cfg, freq);



%% plot baseline MSEn across scales
close all

BASELO = -0.5;
BASEHI = -0.25;
basetind = mseavg.time > BASELO & mseavg.time < BASEHI

lincols = ['r', 'b'];

figure; hold on
clear s
for icond = 1:2 %:4%:4% 1:4 %[1,2,4] %3:4 %[1,3, 4]  % %1:4 %% 1:4
  for istim = 1 %:2%:3
    for iresp = 1%:2
      for imse = 2 % 1:2
        cfg = [];
        %                 freq.mask = poolstat{ imse, icond, istim, iresp }.mask;
        %                                 freq.mask = poolstat{ imse, icond, istim, iresp }.negclusterslabelmat == 1;
        %                                 freq.mask = poolstat{ imse, icond, istim, iresp }.posclusterslabelmat == 1;
        %                 cfg.maskparameter = 'mask';
        %                 if ~any(freq.mask)
        %                     continue;
        %                 end
        %             freq.time = respavg.time{itrig};
        %             freq.freq = respavg.freq{iband};
        freq.powspctrm = squeeze( mseavg.datraw(:,4, imse, :,:,:,  icond, istim, iresp) );
        %                 freq.powspctrm = squeeze( mseavg.dat(:,4, imse, :,:,:,  icond, istim, iresp) );
        freq.label = mseavg.label;
        freq.freq = mseavg.scales;
        freq.time = mseavg.time;
        freq.dimord = 'subj_chan_freq_time';
        
        basedat = mean(freq.powspctrm(:,:,:,basetind),4);
        basedat = squeeze(mean(basedat,2));
        
        s(icond) = shadedErrorBar(mseavg.scales, nanmean(basedat), nanstd(basedat) / sqrt(size(basedat,1)), lincols(icond), 1)
        
        %                 %             cfg.layout = 'biosemi64incI1I2.lay';  % biosemi64  elec1010
        %                 % cfg.layout = 'biosemi64.lay';  % biosemi64  elec1010
        %                 cfg.layout = 'elec1010.lay';  % biosemi64  elec1010
        %                 cfg.layout = ft_prepare_layout(cfg);
        %                 cfg.layout.width(:) = 0.075;
        %                 cfg.layout.height(:) = 0.075;
        %
        %                 cfg.xlim = [-0.75 2];
        % %                 cfg.zlim = [-0.1 0.1];
        %                 cfg.zlim = 'maxabs';
        %                 % cfg.zlim = [-1 1];
        %                 load( 'colormap_jetlightgray.mat')
        %                 cfg.colormap = cmap;
        %                 cfg.hotkeys = 'yes';
        %                 cfg.fontsize = 18;
        %                 cfg.colorbar = 'yes';
        %
        %                 f = figure;
        %                 f.Position = [ 680   678   1200   1200];
        %
        %                 ft_multiplotTFR(cfg, freq);
        
        title(sprintf('%s %s %s', mseavg.mse_leg{imse}, mseavg.sdt_conds{istim,iresp}, mseavg.behav_conds{icond}))
      end
    end
  end
end
legend([s.mainLine],mseavg.behav_conds(1:2), 'Location', 'East')
legend boxoff

%
%






%% collapse over relevant scales











%%
% critEEG_conditionlabels

load('critEEG_chlabel.mat')

% critEEG_sensorselection
% 'Pz'   'Fz'   'Cz'    missing
left_front = {'Fp1' 'AF3' 'F7'  'F3'  };
left_centrfron = { 'FC5'  'FC1' 'C3'};
left_centrpar = { 'T7' 'TP7'  'CP5'  'CP3'  'CP1' };
left_par = {'P9' 'P7'  'P5' 'P3' 'P1' 'PO7'   };
left_occ = { 'PO3'  'O1' 'I1' };

mid_occ = { 'Iz'   'Oz'   'POz' };

right_occ = { 'I2' 'O2'  'PO4' };
right_par = { 'PO8' 'P2' 'P4' 'P6' 'P8' 'P10'  };
right_centrpar = { 'CP2'  'CP4'  'CP6'  'TP8' 'T8'  };
right_centrfron = { 'C4' 'FC2' 'FC6'};
right_front = {'F4' 'F8' 'AF4'  'Fp2'};

lay_1D = [left_front left_centrfron left_centrpar  left_par left_occ  mid_occ  right_occ  right_par right_centrpar right_centrfron right_front];
[~, lay_1Dind]= match_str(lay_1D, mseavg.label);

lay_1D_leg =  {'frontal' 'centrofr.' 'centropar.' 'parietal'  'oc.parietal'  'occipital' 'oc.parietal'  'parietal'  'centropar.' 'centrofr.' 'frontal'}



%% Plot chan X scale heatmaps
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/plots/mse';
mkdir(PREOUT)

close all
SAV = 1;
cmap_RuBd = flipud(cbrewer('div', 'RdBu',256));
cmap_Reds = cbrewer('seq', 'Reds',256);

load colormap_jetlightgray.mat

showstats = 1;

set(0, 'defaultaxesfontsize', 14)
for istim = 1%:2 % 1 = fig, 2 = hom
  
  for iresp = 1%:2 % resp or noresp
    %         if iresp > 1 && strcmp(triggers{itrg}, 'resp')
    %             continue % no button presses for resp = 2
    %         end
    %         colormap(cmap_RuBd)
    
    for icond = 3:4%:4 %[1,2,4] % crit
      %             for itoi = 1:3%:3 %3:4 %[1,2,4]
      
      figure; iplot = 0; ax = [];
      set(gcf, 'Position', [100 150 1000 1000])
      temp = [];
      for imse = 2 %1:2
        for isoi = 1:6
          iplot = iplot+1;
          ax(iplot) = subplot(3, 2, iplot);
          hold on
          %                 if itoi < 3 || icond < 3
          %                 SCALE = [-0.25 0.25];
          %                     if itoi == 3 %|| icond == 3
          %                         SCALE = [-0.25 0.25];
          %                     elseif itoi == 4;
          %                         SCALE = [-0.5 0.5];
          %                     end
          SCALE = [-0.1 0.1];
          if icond == 4
            SCALE = [-0.05 0.05];
          end
          %                     SCALE = [1.1 1.25];
          %                     SCALE = [1 1.25];
          
          dat = squeeze( mean(mseavg.dat(:,4, imse, mseavg.sens.ind{isoi},:,:,  icond, istim, iresp), 4 ));
          dat = squeeze(mean(dat));
          
          %                     % select msen at 0.5 s, as before
          %                     dat = squeeze( mean(mseavg.dat(:,4, imse, :,:,17,  icond, istim, iresp)));
          %                     figure; plot(dat'); colorbar
          %                     temp(:,:,imse) = squeeze( nanmean(mseavg.dat(:,4, imse, itoi,:,:,  icond, istim, iresp)) )';
          %                     if min(min(dat)) < 0
          %                         colormap(ax(iplot), cmap_RuBd);
          %                     else
          %                                             colormap(ax(iplot), cmap_Reds);
          %                     end
          %                     colormap(cmap_Reds)
          colormap(cmap)
          
          %                     if itoi == 3
          % %                         showstats = 1;
          %                         showstats = 0;
          %                     end
          if showstats
            mask = double(squeeze( poolstat{ imse, isoi, icond, istim, iresp }.mask));
            mask(mask==0) = 0.25;
            % mask = mask(lay_1Dind,:)';
            
            ft_plot_matrix(mseavg.time, mseavg.scales,dat, 'clim', SCALE, 'box', 'no', ...
              'highlight', mask, 'highlightstyle', 'opacity'); % opacity
            contour(mask,1)
          else
            ft_plot_matrix(mseavg.time, mseavg.scales, dat, 'clim', SCALE, 'box', 'yes' ); % opacity
          end
          
          %figure; imagesc(mseavg.time, mseavg.scales, dat, [-0.07 0.07]); colorbar
          
          ylim([1 42]);
          plot([0.16 0.16], get(gca, 'ylim'), '--k');
          plot([0 0], get(gca, 'ylim'), '-k');
          plot([1 1], get(gca, 'ylim'), '-k');
          
          xlim( [ mseavg.time(1)  mseavg.time(end)] )
          colorbar
          ax = gca;
          ax.XTick = -1:0.25:2;
          % ax.XTickLabel = mseavg.label(lay_1Dind);
          ax.YDir = 'reverse';
          ax.FontSize = 14;
          %                     xtickangle(45)
          
          
          title(sprintf('%s %s %s %s', mseavg.mse_leg{imse}, mseavg.sdt_conds{istim,iresp}, mseavg.behav_conds{icond}, mseavg.sens.leg{isoi}))
          %                     ylabel('Scale No.')
          %
          %                     %plot scale X mse plot
          %                     iplot = iplot+1;
          % %                     ax(iplot) = subplot(2, 2, iplot);
          %                     subplot(2, 2, iplot);
          %                     chans = 1:45 %28:35 %20:27
          %                     ax = plot(dat(:,chans));
          %                     ax(:).LineWidth = 2;
          %                     ylim([0 1.2])
          %                     if imse==2
          %                         % avg across scales, rank order channels by sorting, check if rank
          %                         % stability changes
          %                         temp2 = squeeze(mean(temp));
          %                         [~, ind] = sort(temp2);
          %                         [r_rank, p_rank] = corr(ind(:,1), ind(:,2));
          %
          %                         [r,p] = corr(temp(:,:,1), temp(:,:,2));
          %                         title(sprintf('mean r = %g across 48 chan\nchan rank r = %g', mean(diag(r)), r_rank))
          %                     end
          xlabel('Time from stim')
          ylabel('Scaleno.')
          %                     legend(lay_1D(chans))
          box on
        end
        outfile = sprintf('%s_%s_%s_%s', mseavg.mse_leg{imse}, mseavg.behav_conds{icond}, mseavg.sdt_conds{istim,iresp} )
        if SAV
          %                     export_fig( fullfile(PREOUT, outfile), '-pdf')
          export_fig( fullfile(PREOUT, outfile), '-png')
          cd(PREOUT)
        end
      end
    end
  end
end


%% plot topo's for each scale

SAV = 1;
close all

cfg = [];
% cfg.layout = 'elec1010.lay';
cfg.comment = 'no';
cfg.marker = 'off';
% cfg.shading = 'flat';
% cfg.style = 'straight'; %both
% cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfg.markersize = 1;
cfg.highlight = 'off';
%                         cfg.highlightchannel = mseavg.label(any(any(plotstat.posclusterslabelmat == iclus, 2),3));
cfg.parameter = 'powspctrm';
cfg.interactive = 'no';

PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/critEEG/plots/mse';
mkdir(PREOUT)

cfg2=[];
cfg2.layout= 'elec1010.lay';
cfg.layout = ft_prepare_layout(cfg2);

cmap_RuBd = flipud(cbrewer('div', 'RdBu',256));
cmap_Reds = cbrewer('seq', 'Reds',256);

for icond = 3:4 %3:4 % [1,2,4] % crit
  figure; iplot = 0; ax = [];
  set(gcf, 'Position', [100 150 1920 1050])
  for istim = 1%:2 % 1 = fig, 2 = hom
    for iresp = 1%:2 % resp or noresp
      
      for itoi = 1:3 %4 %[1,2,4]
        %                 if itoi==1 || itoi==3
        %                 end
        
        for iscale = [1,2,4,8,16,32]  % 1:42 %  1:7:42
          iplot = iplot+1;
          %                     subplot(4, 6, iplot)
          ax(iplot) = subplot(4, 6, iplot);
          %                     ZLIM = [0.8 1.2];
          if itoi > 2 || icond == 4
            colormap(ax(iplot), cmap_RuBd);
            ZLIM = [-0.15 0.15];
            if icond == 4
              ZLIM = [-0.1 0.1];  %[-0.1 0.1];
            end
            cfg.zlim = ZLIM;
          else
            colormap(ax(iplot), cmap_Reds);
            try cfg = rmfield(cfg, 'zlim'); catch end
          end
          
          freq = [];
          freq.label = mseavg.label; % mseavg.label(LR_subtract_mat(:,2));
          freq.dimord = 'chan';
          freq.powspctrm = squeeze( nanmean(mseavg.dat(:,4, imse, itoi,:,iscale,  icond, istim, iresp)) );  %mseavg(isub,ises,3, 1:2,:,1:40,  icond, istim, iresp)
          ft_topoplotTFR(cfg, freq);
          %                     colorbar
          
          t = title(sprintf('%s, %s, %s, sc %d, %g Hz\n%s', mseavg.mse_leg{imse}, mseavg.sdt_conds{istim, iresp}, behav_conds{icond}, iscale, 256/iscale, mseavg.toi_leg{itoi}));
          
        end
      end
      %             outfile = sprintf('topo_%s_%s_%s', mse_leg{imse}, behav_conds{icond}, toi_leg{itoi} )
      outfile = sprintf('topo_%s_%s_%s_%s_%s', mseavg.fun2run, mseavg.mse_leg{imse}, sdt_conds{istim,iresp}, behav_conds{icond}, mseavg.toi_leg{itoi} )
      if SAV
        %             export_fig( fullfile(PREOUT, outfile), '-pdf')
        export_fig( fullfile(PREOUT, outfile), '-png')
      end
      
    end
  end
end
cd(PREOUT)


%% Plot figure for HBM poster

close all
SAV = 1;
cmap_RuBd = flipud(cbrewer('div', 'RdBu',256));
cmap_Reds = cbrewer('seq', 'Reds',256);
for istim = 1%:2 % 1 = fig, 2 = hom
  
  for iresp = 1%:2 % resp or noresp
    %         if iresp > 1 && strcmp(triggers{itrg}, 'resp')
    %             continue % no button presses for resp = 2
    %         end
    %         colormap(cmap_RuBd)
    
    for icond = 3%:4%:4 %[1,2,4] % crit
      for itoi = 1:2%:3 %3:4 %[1,2,4]
        figure; iplot = 0; ax = [];
        %                 set(gcf, 'Position', [0 -200 1400 1050])  % 2143          33        1400        1050
        set(gcf, 'Position', [2143          33        1400        1050])  %
        %                 set(gcf, 'Position', [100 150 1000 500])
        set(0, 'defaultaxesfontsize', 20)
        
        temp = [];
        for imse = 2 %1:2
          
          iplot = iplot+1;
          ax(iplot) = subplot(2, 2, iplot);
          hold on
          %                 if itoi < 3 || icond < 3
          SCALE = [1 1.25];
          %                 SCALE = [-0.25 0.25];
          if itoi == 3 %|| icond == 3
            SCALE = [-0.15 0.15];
          elseif itoi == 4;
            SCALE = [-0.5 0.5];
          end
          if icond == 4
            SCALE = [-0.1 0.1];
          end
          
          dat = squeeze( nanmean(mseavg.dat(:,4, imse, itoi,lay_1Dind,:,  icond, istim, iresp)) )';
          temp(:,:,imse) = squeeze( nanmean(mseavg.dat(:,4, imse, itoi,:,:,  icond, istim, iresp)) )';
          if min(min(dat)) < 0
            colormap(ax(iplot), cmap_RuBd);
          else
            colormap(ax(iplot), cmap_Reds);
          end
          
          showstats = 0;
          if itoi == 3
            %                         showstats = 1;
            showstats = 1;
          end
          if showstats
            mask = double(squeeze( poolstat{ imse, itoi, icond, istim, iresp }.mask));
            mask(mask==0) = 0.25;
            mask = mask(lay_1Dind,:)';
            
            ft_plot_matrix(dat, 'clim', SCALE, 'box', 'no', ...
              'highlight', mask, 'highlightstyle', 'opacity'); % opacity
            contour(mask,1)
          else
            ft_plot_matrix(dat, 'clim', SCALE, 'box', 'no' ); % opacity
          end
          
          ylim([1 42]);
          
          ax = gca;
          ax.XTick = 1:length(lay_1Dind);
          ax.XTickLabel = mseavg.label(lay_1Dind);
          ax.YDir = 'reverse';
          %ax.FontSize = 9;
          %                     xtickangle(45)
          
          ax.XTick = [1  6  11 15 19 ...
            23 ...
            27 31 35 40 45];
          %                     ax.XTick = [1:5:45];
          ax.XTickLabel = lay_1D_leg;
          ax.XTickLabelRotation = 45;
          ax.XLim = [1 45];
          
          %                     title(sprintf('%s %s %s %s', mseavg.mse_leg{imse}, sdt_conds{istim,iresp}, behav_conds{icond}, mseavg.toi_leg{itoi}))
          ylabel('SLOW        Time scale         FAST')
          th = text(-2, 56, 'ANTERIOR LEFT'); th.FontSize = 20;
          th = text(17, 56, 'POSTERIOR'); th.FontSize = 20;
          th = text(32, 56, 'ANTERIOR RIGHT'); th.FontSize = 20;
          
          colh =  colorbar
          %                   colh.Position = [0.4365 0.5933 0.0143 0.3314]
          colh.Position = [0.48 0.71 0.0143 0.1]
          colh.Box = 'off';
          %                     colh.Label.String = 'Entropy (Post - pre stimulus onset)';
          colh.Label.String = 'Entropy';
          colh.Label.FontSize = 20;
          colh.Ticks = SCALE;
          box on
        end
        outfile = sprintf('%s_%s_%s_%s', mseavg.mse_leg{imse}, behav_conds{icond}, mseavg.sdt_conds{istim,iresp}, mseavg.toi_leg{itoi} )
        if SAV
          %                     export_fig( fullfile(PREOUT, outfile), '-pdf')
          export_fig( fullfile(PREOUT, outfile), '-png')
          %             export_fig( fullfile(PREOUT, outfile), '-png')
          cd(PREOUT)
        end
      end
    end
  end
end

