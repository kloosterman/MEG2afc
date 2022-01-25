function [megdat] = MEG2afc_mergefreq_corrstats(megdat)
% Run stats, add them to the megdat struct
% SUBJ_inc = [1:5, 7:9, 11:21]; % all 19 subj

stattype = 'cluster'; % cluster PLS
SUBJ_inc = [1:5, 7:9, 11:21]; % NK nrs! all subj, 19 total
% SUBJ_inc = [1, 3:5, 7:9, 11:21]; % NK nrs! NK2 out, incomplete, 18 total
incsubj = ismember(megdat.SUBJ, SUBJ_inc);

idrug = 4 % change change correlations!
ises = 3
idiff = 3

tiletype = 'timefreq_together'; % for corr with drift/dprime
latencies = [-0.1 0.15; -0.45 0.25]; % p = 0.059 1000 perms

% tiletype = 'separate'; % for resplatr vs p repeat
% latencies = [-0.2 0.3; -0.3 0.2]; % for latr correlation

% tiletype = 'freq_together';

% latencies = [-0.1 0.3; -0.25 0.25]; % rising edge RT is at 0.3, mode RT at 0.65
% latencies = [-0.2 0.4; -0.4 0.2];
% latencies = [-0.2 0.5; -0.5 0.2];
% latencies = [-0.1 0.15; -0.45 0.15]; 

behavnames = {...
  {'ddm_acc_perrun' 'v'};
  %   %     {'ddm_acc_perrun' 'a'};
  %   %     {'ddm_acc_perrun' 't'};
%       {'dprime'};
  %     {'RT'};
  %     {'criterion'};
%     {'p_repeatbalanced'};
  %     {'p_repeatunbalanced'};
%         {'ddm_histbias_perrun' 'histshift_dc'};
  %       {'ddm_histbias_perrun' 'histshift_z'};
  
  
  % old or less interesting:
  %   %   {'button_bias'};
  %   {'chi_accuracy_basic_runs' 'v'};
  %     {'chi_prevresp_z_dc_runs' 'histshift_dc'};
  %       {'chi_prevresp_z_dc' 'histshift_dc'};
  %       {'chi_prevresp_z_dc' 'histshift_z'};
  %   {'chi_prevresp_z_dc_nomotor' 'histshift_z'}; % ses = 1 for this
  %   {'chi_prevresp_z_dc_nomotor' 'histshift_dc'};
  %   {'chi_prevresp_z_dc_runs' 'histshift_dc'}; % OLD, has dc and v
  
  %   {'ddm_histbias_perrun_ol' 'histshift_dc'};
  %   {'ddm_histbias_perrun_ol' 'histshift_z'};
  
  };

behavnames_ctrl = {... % control for these vars
  {'ddm_acc_perrun' 'a'}
  {'ddm_acc_perrun' 't'}
  };
% behavnames_ctrl = {... % control for these vars
% %   {'ddm_histbias_perrun' 'histshift_z'}
% %   {'ddm_histbias_perrun' 'a'}
% %   {'ddm_histbias_perrun' 't'}
% %   {'ddm_histbias_perrun' 'v'}
% };

cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = 'ctf275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg0_neighb);

cfg = [];
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
cfg.spmversion = 'spm12'

switch tiletype
  case 'separate'
    for imod = 3%1:4 % 2 is latr %  2=stim, 3=button, 4 = prevbutton
      for ibehav = 1:length(behavnames)
        for ifreq = 1:2 %2:-1:1
          for itrig = 1:2% :2
            cfg.latency = latencies(itrig,:);
            design = getfield(megdat.behavior, behavnames{ibehav}{:});
            if length(size(design)) == 4
              design = squeeze(design(:,end,idrug,ises)');
            elseif length(size(design)) == 5
              design = squeeze(design(:,end,idrug,ises,idiff)'); % TODO ok for all?
              %       elseif length(size(design)) == 3
              %         cfg.design = squeeze(design(:,idrug,ises)');
            end
            if ~isempty(behavnames_ctrl)
              disp('Adding control variables')
              for ictrl = 1:size(behavnames_ctrl,1)
                ctrl = getfield(megdat.behavior, behavnames_ctrl{ictrl}{:});
                if length(size(ctrl)) == 4
                  ctrl = squeeze(ctrl(incsubj,end,idrug,ises)');
                elseif length(size(ctrl)) == 5
                  ctrl = squeeze(ctrl(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
                end
                design = [design; ctrl]; % append ctrl variable
              end
            end
            
            curfreq = ft_selectdata(cfg, megdat.freq(idrug,imod,itrig,ifreq,idiff)); % for latency cut
            curfreq.powspctrm = curfreq.powspctrm(incsubj,:,:,:);
            cfg.design = design(incsubj);
            
            corrstat{imod, ibehav,ifreq,itrig} = ft_freqstatistics(cfg, curfreq);
            
            disp 'get single subj values for scatter plot correlation'
            posmask1 = corrstat{imod, ibehav,ifreq,itrig}.posclusterslabelmat == 1;
            corrstat{imod, ibehav,ifreq,itrig}.scatterdat = mean(curfreq.powspctrm(:,posmask1(:)),2);

            % keep track of data etc
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm_subj = curfreq.powspctrm;
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(curfreq.powspctrm));
            corrstat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            corrstat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
            corrstat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.corrstat = corrstat;
    
  case 'timefreq_together'
    %% stim, resp and low high together
    disp 'stats'
    corrstat = {};
    for imod = 1 %1:4 % 2 is latr %  2=stim, 3=button, 4 = prevbutton
      for ibehav = 1:length(behavnames)
        % prepare design
        design = getfield(megdat.behavior, behavnames{ibehav}{:});
        if length(size(design)) == 4
          design = squeeze(design(incsubj,end,idrug,ises)');
        elseif length(size(design)) == 5
          design = squeeze(design(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
          %       elseif length(size(design)) == 3
          %         cfg.design = squeeze(design(:,idrug,ises)');
        end
        if ~isempty(behavnames_ctrl)
          disp('Adding control variables')
          for ictrl = 1:size(behavnames_ctrl,1)
            ctrl = getfield(megdat.behavior, behavnames_ctrl{ictrl}{:});
            if length(size(ctrl)) == 4
              ctrl = squeeze(ctrl(incsubj,end,idrug,ises)');
            elseif length(size(ctrl)) == 5
              ctrl = squeeze(ctrl(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
            end
            design = [design; ctrl]; % append ctrl variable
          end
        end
        
        % prepare data
        curfreq = {};
        for ifreq = 1:2
          cfg = [];
          cfg.latency = latencies(1,:);
          freq_sr={};
          freq_sr{1} = ft_selectdata(cfg, megdat.freq(idrug,imod,1,ifreq,idiff)); % for latency cut
          cfg.latency = latencies(2,:);
          freq_sr{2} = ft_selectdata(cfg, megdat.freq(idrug,imod,2,ifreq,idiff)); % for latency cut
          resptimeshift = freq_sr{1}.time(end) + abs(freq_sr{2}.time(1)) + 0.05;
          freq_sr{2}.time = freq_sr{2}.time + resptimeshift;
          cfg = [];
          cfg.appenddim = 'time';
          curfreq{ifreq} = ft_appendfreq(cfg, freq_sr{:});
        end
        cfg = [];
        cfg.appenddim = 'freq';
        curfreq = ft_appendfreq(cfg, curfreq{:});
        curfreq.powspctrm = curfreq.powspctrm(incsubj,:,:,:);
    
        switch stattype
          case 'cluster'
            cfg = [];
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
            cfg.design           = design;
            cfg.minnbchan        = 0;
            cfg.spmversion = 'spm12'

            tempstat = ft_freqstatistics(cfg, curfreq);
          case 'PLS'
            cfg = [];
            cfg.type             = 'Pearson';
            cfg.design           = design;
            cfg.method           = 3;  % Regular Behavior PLS
            cfg.LVthreshold      = [-3 3];
            tempstat = ft_freqstatistics_PLS(cfg, curfreq);
        end
        disp 'get single subj values for scatter plot correlation'
        posmask1 = tempstat.posclusterslabelmat == 1;
        scatterdat = mean(curfreq.powspctrm(:,posmask1(:)),2);
%         figure; scatter(scatterdat, cfg.design(1,:));
%         title(sprintf('r = %1.2f, rho = %1.2f', corr(scatterdat, cfg.design(1,:)', 'type', 'Pearson'), corr(scatterdat, cfg.design(1,:)', 'type', 'Spearman' )))
        
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
            tempfreq = ft_selectdata(cfg, curfreq);
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm_subj = tempfreq.powspctrm;
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(tempfreq.powspctrm));
            corrstat{imod, ibehav,ifreq,itrig}.scatterdat = scatterdat;
            clear tempfreq
            corrstat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            corrstat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
            corrstat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.corrstat = corrstat;
    % megdat.corrstat{1}.posclusters(1).prob
    
  case 'time_together'
    %% take stim and resp together, separate for low and high
    disp 'stats'
    corrstat = {};
    for imod = 1%1:4 % 2 is latr %  2=stim, 3=button, 4 = prevbutton
      for ibehav = 1:length(behavnames)
        for ifreq = 1:2% 1:2 %2:-1:1
          cfg = [];
          cfg.uvar     = [];
          cfg.ivar     = 1;
          cfg.method           = 'montecarlo';
          cfg.statistic        = 'ft_statfun_correlationT';  %depsamplesT
          cfg.type             = 'Pearson'; % Spearman Pearson
          % cfg.type             = 'Spearman'; %  Pearson
          cfg.correctm         = 'cluster';
          %     cfg.correctm         = 'no';
          cfg.clusteralpha     = 0.05;
          cfg.clusterstatistic = 'maxsum';
          cfg.tail             = 0;
          cfg.clustertail      = 0;
          cfg.alpha            = 0.025;
          cfg.numrandomization = 100;
          cfg.neighbours       = neighbours;
    
          design = getfield(megdat.behavior, behavnames{ibehav}{:});
          if length(size(design)) == 4
            cfg.design = squeeze(design(:,end,idrug,ises)');
          elseif length(size(design)) == 5
            cfg.design = squeeze(design(:,end,idrug,ises,idiff)'); % TODO ok for all?
            %       elseif length(size(design)) == 3
            %         cfg.design = squeeze(design(:,idrug,ises)');
          end
          %         curfreq = ft_selectdata(cfg, megdat.freq(idrug,imod,itrig,ifreq,idiff)); % for latency cut
    
          cfg.latency = [-0.2 0.4];
          freq_sr={};
          freq_sr{1} = ft_selectdata(cfg, megdat.freq(idrug,imod,1,ifreq,idiff)); % for latency cut
          cfg.latency = [-0.4 0.2];
          freq_sr{2} = ft_selectdata(cfg, megdat.freq(idrug,imod,2,ifreq,idiff)); % for latency cut
          freq_sr{2}.time = freq_sr{2}.time + 0.85;
          c=[];
          c.appenddim = 'time';
          curfreq=ft_appendfreq(c, freq_sr{:});
          cfg = rmfield(cfg, 'latency');
    
          curfreq.powspctrm = curfreq.powspctrm(incsubj,:,:,:);
          cfg.design = cfg.design(incsubj);
    
          %         corrstat{imod, ibehav,ifreq,itrig} = ft_freqstatistics(cfg, curfreq);
          tempstat = ft_freqstatistics(cfg, curfreq);
    
          for itrig=1:2
            cfg=[];
            if itrig==1
              cfg.latency = [-0.2 0.4];
            else
              cfg.latency = [-0.4 0.2] + 0.85;
            end
            corrstat{imod, ibehav,ifreq,itrig} = ft_selectdata(cfg, tempstat);
            if itrig==2; corrstat{imod, ibehav,ifreq,itrig}.time = corrstat{imod, ibehav,ifreq,itrig}.time - 0.85; end
    
            % keep track of data etc
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm_subj = freq_sr{itrig}.powspctrm;
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(freq_sr{itrig}.powspctrm));
            corrstat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            corrstat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
            corrstat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.corrstat = corrstat;

  case 'freq_together'
    %% low high freq together, separately for stim and resp
    disp 'stats'
    corrstat = {};
    for imod = 5 %1:4 % 2 is latr %  2=stim, 3=button, 4 = prevbutton
      for ibehav = 1:length(behavnames)
        % prepare design
        design = getfield(megdat.behavior, behavnames{ibehav}{:});
        if length(size(design)) == 4
          design = squeeze(design(incsubj,end,idrug,ises)');
        elseif length(size(design)) == 5
          design = squeeze(design(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
          %       elseif length(size(design)) == 3
          %         cfg.design = squeeze(design(:,idrug,ises)');
        end
        if ~isempty(behavnames_ctrl)
          disp('Adding control variables')
          for ictrl = 1:size(behavnames_ctrl,1)
            ctrl = getfield(megdat.behavior, behavnames_ctrl{ictrl}{:});
            if length(size(ctrl)) == 4
              ctrl = squeeze(ctrl(incsubj,end,idrug,ises)');
            elseif length(size(ctrl)) == 5
              ctrl = squeeze(ctrl(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
            end
            design = [design; ctrl]; % append ctrl variable
          end
        end
    
        curfreq = {};
        for ifreq = 1:2
          curfreq{ifreq} = megdat.freq(idrug,imod,1,ifreq,idiff);
        end
        cfg = [];
        cfg.appenddim = 'freq';
        curfreq = ft_appendfreq(cfg, curfreq{:});
        curfreq.powspctrm = curfreq.powspctrm(incsubj,:,:,:);
    
        switch stattype
          case 'cluster'
            cfg = [];
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
            cfg.design           = design;
            cfg.minnbchan        = 0;
            tempstat = ft_freqstatistics(cfg, curfreq);
          case 'PLS'
            cfg = [];
            cfg.type             = 'Pearson';
            cfg.design           = design;
            cfg.method           = 3;  % Regular Behavior PLS
            cfg.LVthreshold      = [-3 3];
            tempstat = ft_freqstatistics_PLS(cfg, curfreq);
        end
    
        disp('break up parts again')
        for itrig=1:2
          for ifreq = 1:2
            cfg=[];
%             cfg.latency = latencies(itrig,:);
%             if itrig==2
%               cfg.latency = cfg.latency + resptimeshift;
%             end
            if ifreq==1
              cfg.frequency = [2 35];
            else
              cfg.frequency = [36 100];
            end
            corrstat{imod, ibehav,ifreq,itrig} = ft_selectdata(cfg, tempstat);
%             if itrig==2; corrstat{imod, ibehav,ifreq,itrig}.time = corrstat{imod, ibehav,ifreq,itrig}.time - resptimeshift; end
    
            % keep track of data etc
            tempfreq = ft_selectdata(cfg, curfreq);
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm_subj = tempfreq.powspctrm;
            corrstat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(tempfreq.powspctrm));
            clear tempfreq
            corrstat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            corrstat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
            corrstat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.corrstat = corrstat;
    % megdat.corrstat{1}.posclusters(1).prob
    
    
end



%% old hacked code for starters, remove soon
%     cfg.design   = rep_probs(:,4)';
%     cfg.design   = db(:,4)';
%     cfg.design   = z(:,4)';
%     cfg.design   = design;

% behavpath = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/MEG_HH_analysis/behavior/historybias';
% % % load Anke history bias data
% % rep_probs = csvread(fullfile(behavpath, 'rep_probs.csv'), 1,1) % atx, plac
% % NK's debugged:
% load(fullfile(behavpath, 'rep_probsNK.mat'))
%
% rep_probs(:,4) = rep_probs(:,2) - rep_probs(:,1);
%
% % load history ddm
% ddmpars = MEG2afc_load_ddm_paras('/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_histbias.csv')
% db(:,1) = ddmpars.dc_1_1 - ddmpars.dc_1_0; %drug
% db(:,2) = ddmpars.dc_0_1 - ddmpars.dc_0_0 % plac
% db(:,4) = db(:,1) - db(:,2);
% mean(db)
% [h,p] = ttest(db(:,1), db(:,2))
% [h,p] = ttest(rep_probs(:,1), rep_probs(:,2))
% r = corr(db(:,1), rep_probs(:,1))
% figure; subplot(131); scatter(db(:,1), rep_probs(:,1)); title(r)
% r = corr(db(:,2), rep_probs(:,2))
% subplot(132); scatter(db(:,2), rep_probs(:,2)); title(r)
% r = corr(db(:,4), rep_probs(:,4))
% subplot(133); scatter(db(:,4), rep_probs(:,4)); title(r)
%
% %starting point
% z(:,1) = ddmpars.z_1_1 - ddmpars.z_1_0; %drug
% z(:,2) = ddmpars.z_0_1 - ddmpars.z_0_0 % plac
% z(:,4) = z(:,1) - z(:,2);
% mean(z)
% [h,p] = ttest(z(:,1), z(:,2))
% [h,p] = ttest(rep_probs(:,1), rep_probs(:,2))
% r = corr(z(:,1), rep_probs(:,1))
% figure; subplot(131); scatter(z(:,1), rep_probs(:,1)); title(r)
% r = corr(z(:,2), rep_probs(:,2))
% subplot(132); scatter(z(:,2), rep_probs(:,2)); title(r)
% r = corr(z(:,4), rep_probs(:,4))
% subplot(133); scatter(z(:,4), rep_probs(:,4)); title(r)
%
% % [r,p] = corr(db(:,4), z(:,4));
% % figure; scatter(db(:,4), z(:,4));
%
% if length(megdat.SUBJ) ~= 19
%   % SUBJ= [1, 3, 5, 7:9, 11:13, 15:21]; % NK2 incomplete, NK4 artifact drug-plac?, NK14 magnitude bigger responses
%   % 6 and 10 do not exist
%   disp 'remove dropped subj'
%   rep_probs = rep_probs([1,3:end],:); % so 2 out
%   db = db([1,3:end],:); % so 2 out
%   z = z([1,3:end],:); % so 2 out
% end
%
