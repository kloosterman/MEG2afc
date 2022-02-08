function [freq_subj] = MEG2afc_mergefreq_runs(  )
%UNTITLED6 Summary of this function goes here
%   BLC = baseline correction: raw or BLC

if nargin == 0
  baseline = 'alreadynormalized'; % BLC
end

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/'; %yesno or 2afc kloosterman
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
end

PREIN = fullfile(basepath, 'freqzapline-plus'); % freqbandstop freqzapline-plus
cd(PREIN)

disp 'load behavior'
load(fullfile(basepath, 'behav/behavstruct.mat'));

drugleg = {'drug', 'plac'};
motorleg = {'contra' 'ipsi'};

% SUBJ= [1:5, 7:9, 11:21]; % all 19 subj
SUBJ= [1, 3:5, 7:9, 11:21]; % 2 incomplete

megdat = []; % output
megdat.SUBJ = SUBJ;
megdat.PREOUT = fullfile(PREIN, 'plots');  mkdir(megdat.PREOUT);
megdat.dimord = 'drug_mod_trig_freq_diff';

nsub = length(SUBJ);
imod = 1; % modulation
freqall = {};
freq_subj = {};
for isub = 1:nsub
  clear freq_ses
  for is = 1:2
    clear freq_drug
    for idrug = 1:2
      sesfile = dir(sprintf('NK%d_%s_%s_freq.mat', SUBJ(isub), drugleg{idrug}, motorleg{is}));
      fprintf(sesfile.name)
      load(sesfile.name, 'freq') % freq{iltr+1, itrig, ifreq, idiff}
      if ~isfield(freq, 'powspctrm')
        disp('File empty')
        continue
      end
      
      %take difficulty avg
      cfg=[];
      cfg.parameter = 'powspctrm';
      cfg.operation = 'add';
      freq = arrayfun(@(x,y) ft_math(cfg, x,y), freq(imod,:,:,1), freq(imod,:,:,2)); % modulation only
      cfg.operation = 'divide';
      cfg.scalar = 2;
      freq = squeeze(arrayfun(@(x) ft_math(cfg, x), freq));
      %       ft_databrowser([], freq(1))
      
      % put freq together
      freq = squeeze(arrayfun(@(x,y) ft_appendfreq([], x,y), freq(:,1), freq(:,2)));
      
      % select toi stim and resp
      latencies = [-0.1 0.15; -0.45 0.25];
      for itrig = 1:2
        cfg = [];
        cfg.latency = latencies(itrig,:);
        freq(itrig).time = round(freq(itrig).time*100)/100; % to get rid of slight time difs
        freq(itrig) = ft_selectdata(cfg, freq(itrig)); % for latency cut
      end
      resptimeshift = freq(1).time(end) + abs(freq(2).time(1)) + 0.05;
      freq(2).time = freq(2).time + resptimeshift;
      cfg = [];
      cfg.appenddim = 'time';
      freq_drug(idrug) = ft_appendfreq(cfg, freq(1), freq(2));
      clear freq
    end
    % match nruns atx and plac
    nruns = min([size(freq_drug(1).powspctrm,1) size(freq_drug(2).powspctrm,1)]);
    freq_drug(1).powspctrm = freq_drug(1).powspctrm(1:nruns,:,:,:);
    freq_drug(2).powspctrm = freq_drug(2).powspctrm(1:nruns,:,:,:);
    
    %     atx-plac per run
    cfg=[];
    cfg.parameter = 'powspctrm';
    cfg.operation = 'subtract';
    freq_ses{is} = ft_math(cfg, freq_drug(1), freq_drug(2));
    
    % add behavior: DDM pars
    b=behavior.ddm_histbias_perrun;
    isub2 = find(SUBJ(isub) == behavior.SUBJ); % index of subject in ddm file
    freq_ses{is}.trialinfo(:,1) = squeeze(b.v(isub2, 1:nruns, 1, is) - b.v(isub2, 1:nruns, 2, is, 3)); % dimord: 'subj_runs_drug_motor_diffodatselrevresp'
    freq_ses{is}.trialinfo(:,2) = squeeze(b.a(isub2, 1:nruns, 1, is) - b.a(isub2, 1:nruns, 2, is)); % dimord: 'subj_runs_drug_motor_diffodatselrevresp'
    freq_ses{is}.trialinfo(:,3) = squeeze(b.t(isub2, 1:nruns, 1, is) - b.t(isub2, 1:nruns, 2, is)); % dimord: 'subj_runs_drug_motor_diffodatselrevresp'
    freq_ses{is}.trialinfo(:,4) = squeeze(b.histshift_dc(isub2, 1:nruns, 1, is) - b.histshift_dc(isub2, 1:nruns, 2, is)); % dimord: 'subj_runs_drug_motor_diffodatselrevresp'
    freq_ses{is}.trialinfo(:,5) = squeeze(b.histshift_z(isub2, 1:nruns, 1, is) - b.histshift_z(isub2, 1:nruns, 2, is)); % dimord: 'subj_runs_drug_motor_diffodatselrevresp'
    if any(any(isnan(freq_ses{is}.trialinfo)))
      disp('nan found')
    end
  end
  % append
  cfg=[];
  cfg.appenddim = 'rpt';
  freq_ses = ft_appendfreq(cfg,freq_ses{:});
  
  % zscore per subject to put on equal footing
  freq_ses.powspctrm = (freq_ses.powspctrm - mean(freq_ses.powspctrm)) ./ ...
    std(freq_ses.powspctrm);
  freq_ses.trialinfo = (freq_ses.trialinfo - mean(freq_ses.trialinfo)) ./ ...
    std(freq_ses.trialinfo);
  
  %put in big mat
  freq_subj{end+1} = freq_ses;
end

freq_subj = ft_appendfreq(cfg,freq_subj{:});

save freq_subj freq_subj

% cfg=[];
% cfg.layout = 'CTF275.lay'
% ft_multiplotTFR(cfg, freq_subj)

%% run stats


% behavnames = {...
%   {'ddm_acc_perrun' 'v'};
%   %   %     {'ddm_acc_perrun' 'a'};
%   %   %     {'ddm_acc_perrun' 't'};
%   %       {'dprime'};
%   %     {'RT'};
%   %     {'criterion'};
%   %     {'p_repeatbalanced'};
%   %     {'p_repeatunbalanced'};
%   %         {'ddm_histbias_perrun' 'histshift_dc'};
%   %       {'ddm_histbias_perrun' 'histshift_z'};
%   };
%
% behavnames_ctrl = {... % control for these vars
%   {'ddm_acc_perrun' 'a'}
%   {'ddm_acc_perrun' 't'}
%   };
% % behavnames_ctrl = {... % control for these vars
% % %   {'ddm_histbias_perrun' 'histshift_z'}
% % %   {'ddm_histbias_perrun' 'a'}
% % %   {'ddm_histbias_perrun' 't'}
% % %   {'ddm_histbias_perrun' 'v'}
% % };

cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = 'ctf275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg0_neighb);

ibehav = 1;
cfg = [];
cfg.design           = freq_subj.trialinfo(:,ibehav);
cfg.uvar     = [];
cfg.ivar     = 1;
cfg.method           = 'montecarlo';
%         cfg.statistic        = 'ft_statfun_correlationT';  %depsamplesT ft_statfun_correlationT_corrcol
cfg.statistic        = 'ft_statfun_partialcorrelationT';  %depsamplesT ft_statfun_correlationT_corrcol
cfg.type             = 'Pearson'; % Spearman Pearson
cfg.correctm         = 'cluster';  %'no'
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 100;
cfg.neighbours       = neighbours;
cfg.minnbchan        = 0;
cfg.spmversion = 'spm12';

tempstat = ft_freqstatistics(cfg, freq_subj);
% disp 'get single subj values for scatter plot correlation'
% posmask1 = tempstat.posclusterslabelmat == 1;
% scatterdat = mean(curfreq.powspctrm(:,posmask1(:)),2);
%         figure; scatter(scatterdat, cfg.design(1,:));
%         title(sprintf('r = %1.2f, rho = %1.2f', corr(scatterdat, cfg.design(1,:)', 'type', 'Pearson'), corr(scatterdat, cfg.design(1,:)', 'type', 'Spearman' )))

save rmcorr_stat tempstat 

%% plot
clussign = 'pos'

load colormap_jetlightgray.mat
cfg=[];
imod=1;
if imod == 1
  cfg.layout = 'CTF275_helmet.mat';
else
  cfg.layout = 'CTF275_helmet_latr.mat';
end
cfg.colorbar = 'yes';
cfg.colormap = cmap;
cfg.zlim = 'maxabs';
cfg.hotkeys = 'yes';

cfg.maskparameter = 'mask';
cfg.maskalpha = 0.25;
tempstat.mask = tempstat.([clussign 'clusterslabelmat']) == 1; %  tempstat
% tempstat.mask = tempstat.prob < 0.05;

cfg.parameter = 'rho';
% ctr = ctr + 1;
% subplot(2,2,ctr)
%   ft_multiplotTFR(cfg, freq(3,2,1,1)) % idrug, itrig, ifreq, idiff
f = figure;
f.Position = [   680   444   781   654];
ft_multiplotTFR(cfg, tempstat) % idrug, itrig, ifreq, idiff

%%
plotscatter = 1;
if plotscatter
  bcorr = ft_findcfg(tempstat.cfg, 'design')';
  posmask1 = tempstat.posclusterslabelmat == 1;
  scatterdat = mean(freq_subj.powspctrm(:,posmask1(:)),2);
  acorr = scatterdat;
  f = figure; f.Position = [744   950   147   100];
  scatter(acorr, bcorr(:,1), 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 30);
  box on; axis square;  axis tight
  %           set(gca, 'XLim', [-3.3840 3.3840])
  %           set(gca, 'YLim', [-0.36 0.36])
  title(sprintf('r = %1.2f, rho = %1.2f', partialcorr(acorr, bcorr(:,1), bcorr(:,2:end), 'type', 'Pearson'), partialcorr(acorr, bcorr(:,1), bcorr(:,2:end), 'type', 'Spearman' )))
  lsline
  %   saveas(gcf, fullfile(megdat.PREOUT, sprintf('clus%d_%s_vs_%s_%s_scatter.pdf', ifig,  tempstat.megtype, tempstat.behavname{:} ))) %
end

%%
disp('break up parts again')
for itrig=1:2
  for ifreq = 1:2
    cfg=[];
    cfg.latency = latencies(itrig,:);
    if itrig==2
      cfg.latency = cfg.latency + resptimeshift;
    end
    if ifreq==1
      cfg.frequency = [2 35];
    else
      cfg.frequency = [36 100];
    end
    corrstat{imod, ibehav,ifreq,itrig} = ft_selectdata(cfg, tempstat);
    if itrig==2; corrstat{imod, ibehav,ifreq,itrig}.time = corrstat{imod, ibehav,ifreq,itrig}.time - resptimeshift; end
    
    % keep track of data etc
    tempfreq = ft_selectdata(cfg, freq_subj);
    corrstat{imod, ibehav,ifreq,itrig}.powspctrm_subj = tempfreq.powspctrm;
    corrstat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(tempfreq.powspctrm));
    corrstat{imod, ibehav,ifreq,itrig}.scatterdat = scatterdat;
    clear tempfreq
    corrstat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ;
    corrstat{imod, ibehav,ifreq,itrig}.behavname = {'drift'}; % behavnames{ibehav}; % behav name
    corrstat{imod, ibehav,ifreq,itrig}.megtype = 'modulation';% megdat.megleg{imod}; % meg measure name
  end
end

megdat.corrstat = corrstat;
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
      for ifreq = 2:-1:1 % 1:2 %
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
%           cfg.titleTFR = sprintf('%s\n',curstat.behavname{:}, curstat.megtype);;
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
        if irow == 3 || ibehav == size(megdat.corrstat,1) % only 4 fit, starting from 0
          if SAV
            %               saveas(gcf, fullfile(megdat.PREOUT, sprintf('corr_%svs%s.pdf',  megdat.corrstat{ibehav,1,1}.megtype,  [megdat.corrstat{ibehav,1,1}.behavname{:}] )))
            saveas(gcf, fullfile(megdat.PREOUT, sprintf('clus%d_%s_vs_%s_corr.pdf', ifig,  curstat.megtype, curstat.behavname{:} )))
            saveas(gcf, fullfile(megdat.PREOUT, sprintf('clus%d_%s_vs_%s_corr.png', ifig,  curstat.megtype, curstat.behavname{:} )))
          end
          f = figure;      f.Position = [ 680          75         612        792 ]; % A4 formaat
          ifig = ifig+1;
          irow = 0;
        else
          irow = irow+1;
        end
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
          saveas(gcf, fullfile(megdat.PREOUT, sprintf('clus%d_%s_vs_%s_scatter.pdf', ifig,  curstat.megtype, curstat.behavname{:} ))) %
        end
    end
  end
end
cd(megdat.PREOUT)