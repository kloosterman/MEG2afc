function [megdat] = MEG2afc_mergefreq_stats(megdat)
% Run stats, add them to the megdat struct
% SUBJ_inc = [1:5, 7:9, 11:21]; % all 19 subj

stattype = 'cluster'; % cluster PLS
SUBJ_inc = [1:5, 7:9, 11:21]; % NK nrs! all subj, 19 total
% SUBJ_inc = [1, 3:5, 7:9, 11:21]; % NK nrs! NK2 out, incomplete, 18 total
incsubj = ismember(megdat.SUBJ, SUBJ_inc);

ises = 3; idiff = 3;

% latencies = [-0.1 0.3; -0.25 0.25]; % rising edge RT is at 0.3, mode RT at 0.65
tiletype = 'timefreq_together';

% latencies = [-0.2 0.3; -0.3 0.2];
% latencies = [-0.2 0.4; -0.4 0.2];
% latencies = [-0.2 0.5; -0.5 0.2];
% latencies = [-0.1 0.15; -0.45 0.15]; % p = 0.039604
latencies = [-0.1 0.15; -0.45 0.25]; % p = 0.059 1000 perms

cfg0_neighb = [];
cfg0_neighb.method    = 'template';
cfg0_neighb.template  = 'ctf275_neighb.mat';
neighbours       = ft_prepare_neighbours(cfg0_neighb);


switch tiletype
%   case 'separate'
%     for imod = 1%1:4 % 2 is latr %  2=stim, 3=button, 4 = prevbutton
%       for ibehav = 1:length(behavnames)
%         for ifreq = 1:2 %2:-1:1
%           for itrig = 1:2% :2
%             cfg = [];
%             cfg.uvar     = [];
%             cfg.ivar     = 1;
%             cfg.method           = 'montecarlo';
%             %         cfg.statistic        = 'ft_statfun_correlationT';  %depsamplesT ft_statfun_correlationT_corrcol
%             cfg.statistic        = 'ft_statfun_partialcorrelationT';  %depsamplesT ft_statfun_correlationT_corrcol
%             cfg.type             = 'Pearson'; % Spearman Pearson
%             cfg.correctm         = 'cluster';  %'no'
%             cfg.clusteralpha     = 0.05;
%             cfg.clusterstatistic = 'maxsum';
%             cfg.tail             = 0;
%             cfg.clustertail      = 0;
%             cfg.alpha            = 0.025;
%             cfg.numrandomization = 1000;
%             cfg.neighbours       = neighbours;
%             cfg.minnbchan        = 0;
%             cfg.latency = latencies(itrig,:);
%             design = getfield(megdat.behavior, behavnames{ibehav}{:});
%             if length(size(design)) == 4
%               cfg.design = squeeze(design(:,end,idrug,ises)');
%             elseif length(size(design)) == 5
%               cfg.design = squeeze(design(:,end,idrug,ises,idiff)'); % TODO ok for all?
%               %       elseif length(size(design)) == 3
%               %         cfg.design = squeeze(design(:,idrug,ises)');
%             end
%             if ~isempty(behavnames_ctrl)
%               disp('Adding control variables')
%               ctrl = getfield(megdat.behavior, behavnames_ctrl{ibehav}{:});
%               if length(size(ctrl)) == 4
%                 ctrl = squeeze(ctrl(incsubj,end,idrug,ises)');
%               elseif length(size(ctrl)) == 5
%                 ctrl = squeeze(ctrl(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
%               end
%               cfg.design = [cfg.design; ctrl]; % append ctrl variable
%             end
%             
%             curfreq = ft_selectdata(cfg, megdat.freq(idrug,imod,itrig,ifreq,idiff)); % for latency cut
%             curfreq.powspctrm = curfreq.powspctrm(incsubj,:,:,:);
%             cfg.design = cfg.design(incsubj);
%             
%             stat{imod, ibehav,ifreq,itrig} = ft_freqstatistics(cfg, curfreq);
%
%             % keep track of data etc
%             stat{imod, ibehav,ifreq,itrig}.powspctrm_subj = curfreq.powspctrm;
%             stat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(curfreq.powspctrm));
%             stat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
%             stat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
%             stat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
%           end
%         end
%       end
%     end
%     megdat.stat = stat;

  case 'timefreq_together'
    %% stim, resp and low high together
    disp 'stats'
    stat = {};
    for idrug = 4 [2 4] % 2 plac, 4 atx-plac
      for imod = [1] %1:4 % 2 is latr %  2=stim, 3=button, 4 = prevbutton
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
        
        cfg = [];
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'depsamplesT';  %depsamplesT ft_statfun_correlationT ft_statfun_correlationT_corrcol
        cfg.correctm         = 'cluster';  %'no'
        cfg.clusteralpha     = 0.05;
        cfg.clusterstatistic = 'maxsum';
        cfg.tail             = 0;
        cfg.clustertail      = 0;
        cfg.alpha            = 0.025;
        cfg.numrandomization = 100;
        cfg.neighbours       = neighbours;
        cfg.minnbchan        = 0;
        nsub = size(megdat.freq(1).powspctrm,1);
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
        cfg.spmversion = 'spm12'

        freqzero = curfreq; %create zero freq to test against
        freqzero.powspctrm = zeros(size(curfreq.powspctrm));
        tempstat = ft_freqstatistics(cfg, curfreq, freqzero);
        
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
            stat{imod, idrug, ifreq,itrig} = ft_selectdata(cfg, tempstat);
            if itrig==2; stat{imod, idrug, ifreq,itrig}.time = stat{imod, idrug, ifreq,itrig}.time - resptimeshift; end
            
            % keep track of data etc
            tempfreq = ft_selectdata(cfg, curfreq);
            stat{imod, idrug, ifreq,itrig}.powspctrm_subj = tempfreq.powspctrm;
            stat{imod, idrug, ifreq,itrig}.powspctrm = squeeze(mean(tempfreq.powspctrm));
            clear tempfreq
            stat{imod, idrug, ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            stat{imod, idrug, ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.stat = stat;
    % megdat.stat{1}.posclusters(1).prob
    
  case 'time_together'
    %% take stim and resp together, separate for low and high
    disp 'stats'
    stat = {};
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
          
          %         stat{imod, ibehav,ifreq,itrig} = ft_freqstatistics(cfg, curfreq);
          tempstat = ft_freqstatistics(cfg, curfreq);
          
          for itrig=1:2
            cfg=[];
            if itrig==1
              cfg.latency = [-0.2 0.4];
            else
              cfg.latency = [-0.4 0.2] + 0.85;
            end
            stat{imod, ibehav,ifreq,itrig} = ft_selectdata(cfg, tempstat);
            if itrig==2; stat{imod, ibehav,ifreq,itrig}.time = stat{imod, ibehav,ifreq,itrig}.time - 0.85; end
            
            % keep track of data etc
            stat{imod, ibehav,ifreq,itrig}.powspctrm_subj = freq_sr{itrig}.powspctrm;
            stat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(freq_sr{itrig}.powspctrm));
            stat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            stat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
            stat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.stat = stat;
    
  case 'freq_together'
    %% low high freq together, separately for stim and resp
    disp 'stats'
    stat = {};
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
          ctrl = getfield(megdat.behavior, behavnames_ctrl{ibehav}{:});
          if length(size(ctrl)) == 4
            ctrl = squeeze(ctrl(incsubj,end,idrug,ises)');
          elseif length(size(ctrl)) == 5
            ctrl = squeeze(ctrl(incsubj,end,idrug,ises,idiff)'); % TODO ok for all?
          end
          design = [design; ctrl]; % append ctrl variable
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
            stat{imod, ibehav,ifreq,itrig} = ft_selectdata(cfg, tempstat);
            %             if itrig==2; stat{imod, ibehav,ifreq,itrig}.time = stat{imod, ibehav,ifreq,itrig}.time - resptimeshift; end
            
            % keep track of data etc
            tempfreq = ft_selectdata(cfg, curfreq);
            stat{imod, ibehav,ifreq,itrig}.powspctrm_subj = tempfreq.powspctrm;
            stat{imod, ibehav,ifreq,itrig}.powspctrm = squeeze(mean(tempfreq.powspctrm));
            clear tempfreq
            stat{imod, ibehav,ifreq,itrig}.SUBJ_inc = SUBJ_inc;
            stat{imod, ibehav,ifreq,itrig}.behavname = behavnames{ibehav}; % behav name
            stat{imod, ibehav,ifreq,itrig}.megtype = megdat.megleg{imod}; % meg measure name
          end
        end
      end
    end
    megdat.stat = stat;
    % megdat.stat{1}.posclusters(1).prob
    
    
end




%% ori, trig and freq separate
% function [megdat] = MEG2afc_mergefreq_stats(megdat)
% % % % % Run stats, add them to the megdat struct
%
% cfg = [];
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
% cfg.correctm         = 'cluster';
% %     cfg.correctm         = 'no';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% % cfg.alpha            = 0.1;
% cfg.numrandomization = 100;
% cfg0_neighb = [];
% cfg0_neighb.method    = 'template';
% cfg0_neighb.template  = 'ctf275_neighb.mat';
% cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);
% nsub = size(megdat.freq(1).powspctrm,1);
% design = zeros(2,2*nsub);
% for i = 1:nsub
%   design(1,i) = i;
% end
% for i = 1:nsub
%   design(1,nsub+i) = i;
% end
% design(2,1:nsub)        = 1;
% design(2,nsub+1:2*nsub) = 2;
%
% cfg.design   = design;
% cfg.uvar     = 1;
% cfg.ivar     = 2;
%
% stat = {};
% for idiff = 3%:4;
%   for imod = 1%:3 % mod or latr
%     for idrug = 4 %1:4
%       for ifreq = 1:2 %2:-1:1
%         if ifreq==2
%           cfg.frequency = [35 100];
%         else
%           cfg.frequency = 'all';
%         end
%         for itrig = 1:2% :2
%           if itrig == 1
%             cfg.latency = [-0.25 0.4];
%           else
%             cfg.latency = [-0.4 0.2];
%           end
%           freq = ft_selectdata(cfg, megdat.freq(idrug, imod, itrig, ifreq, idiff));
%
%           freqzero = freq; %create zero freq to test against
%           freqzero.powspctrm = zeros(size(freq.powspctrm));
%           stat{idrug, imod, itrig, ifreq, idiff} = ft_freqstatistics(cfg, freq, freqzero);
%           stat{idrug, imod, itrig, ifreq, idiff}.powspctrm = squeeze(mean(freq.powspctrm));
%         end
%       end
%     end
%   end
% end
% megdat.stat = stat;
%
