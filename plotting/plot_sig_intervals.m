function plot_sig_intervals(stat, taxis, iteration, SIGCOL, BARWID)
% plotting
% iteration: for plotting sig lines underneath each other

THR=0.025;
nsteps = 20;

YLIM = get(gca, 'ylim');
ysteps = YLIM(1):diff(YLIM)/nsteps:YLIM(2); % # steps
BARPOS = ysteps(iteration+1); %

sigint={};
i = 1; sigint{i} = [];
pval = stat.prob';  % stat from ft_freqstatistics
for ismp = 1:length(pval)
    if pval(ismp) < THR
        sigint{i} = [sigint{i} ismp];
    end % concatenate significant samples in sigint{iint}
    if (ismp < size(pval,1)) && (pval(ismp+1) >= THR) && ~isempty(sigint{i}) % jump to next interval if next sample is not sig
        i= i + 1;
        sigint{i} = [];
    end
end
% replot significant intervals in different color
for iint = 1:length(sigint)
    if ~isempty(sigint{iint})
        begsmp = sigint{iint}(1); %*2
        endsmp = sigint{iint}(end); %*2
        plot(taxis(begsmp:endsmp),...   %begsmp-1:endsmp+1
            BARPOS*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
            'Color',SIGCOL,'LineWidth',BARWID) %SIGCOL(ircond,idcond,:)
    end
end


% % stats
% cfg=[];
% cfg.channel          = 'p';
% cfg.parameter        = 'trial';
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.025;
% cfg.numrandomization = 1000;
% cfg.neighbours       = []; %in case no channel data present
% 
% design = zeros(2,2*nsub);
% for i = 1:nsub
%     design(1,i) = i;
% end
% for i = 1:nsub
%     design(1,nsub+i) = i;
% end
% design(2,1:nsub)        = 1;
% design(2,nsub+1:2*nsub) = 2;
% 
% cfg.design   = design;
% cfg.uvar     = 1;
% cfg.ivar     = 2;
% 
% timelock.dimord = 'rpt_chan_time';
% timelock1=timelock;
% timelock2=timelock;
% timelockzeros=timelock;
% % ises=3
% 
% if loadstat
%     disp('Loading stat from disk!')
%     load([PREIN 'respavg/stat_' trigger])
%     disp('Done')
% else
% %     stat=cell(2,5,3); % report distr ievent
%     for ircond= 1:4 %[1,2,4]
%         for idcond=1:5 %1:5
%             %         timelock1.trial = squeeze(respavgsubj(:,ises,:,:,itype,1));
%             %         timelock2.trial = squeeze(respavgsubj(:,ises,:,:,itype,2));
%             timelock1.trial = squeeze(respavg(:,:,:,ircond,idcond,1));
%             timelock2.trial = squeeze(respavg(:,:,:,ircond,idcond,2));
%             timelockzeros.trial = zeros(size(timelock1.trial));
%             stat{ircond,idcond,1} = ft_timelockstatistics(cfg, timelock1, timelockzeros);
%             stat{ircond,idcond,2} = ft_timelockstatistics(cfg, timelock2, timelockzeros);
%             stat{ircond,idcond,3} = ft_timelockstatistics(cfg, timelock1, timelock2);
%         end
%     end
%     save([PREIN 'respavg/stat_' trigger], 'stat')
% end

