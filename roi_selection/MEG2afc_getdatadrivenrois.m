% MEG2afc_getdatadrivenrois
% Do space - get mean median RT across subjects
% - TODO put low + high freq into respavg
% - make freq, concatenating low and high freq, stim and resp locked: take half median RT for stim- resp division
% - run space-freq-time freqstats across subjects
%  plot clusterplots of sig clusters, showing average freq, for different time points

%  OR plot clusterplots of sig clusters, integrating over frequencies,
%  showing time points | integrating over time points, showing frequencies

% can run this for power, or for correlation


MEG2afc_load_respavg

%
driftrates = MEG2afc_load_driftrates161116(SUBJ, PREOUT);

%% freqstats space freq time
%
% % stack low and high freq vertically
%
% temp = cat(3,  respavg(:,:, 1:length(frind{1}), 1:length(tind{1}), :, :, :, 1), respavg(:,:, 1:length(frind{2}), 1:length(tind{1}), :, :, :, 2)  );
% % testing
% % figure; imagesc(taxis{1}, faxis_all, squeeze(mean(mean(temp(:,:,:,:,3,3,1)))), [-0.05 0.05])
% close all
% figure; imagesc(taxis{1}, [], squeeze(mean(mean(temp(:,:,:,:,3,3,1)))), [-0.05 0.05])
% ylabel_locations = get(gca, 'YTick');
% set(gca, 'YTicklabel',  {faxis_all(ylabel_locations)})
% set(gca, 'YDir','normal');

%% stack data
% stack low and high freq vertically
temp=[];
temp = cat(3,  respavg(:,:, 1:length(frind{1}), :, :, :, :, 1), respavg(:,:, 1:length(frind{2}), :, :, :, :, 2)  );
% % stack stim and resplocked horizontally
% temp = cat(4,  temp(:,:, :, 1:length(tind{1}), :, :, 1), temp(:,:, :, 1:length(tind{2}), :, :, 2)  );

% close all
% figure; imagesc([], [], squeeze(mean(mean(temp(:,:,:,:,3,3,1)))), [-0.05 0.05])
% ylabel_locations = get(gca, 'YTick');
% set(gca, 'YTicklabel',  {faxis_all(ylabel_locations)})
% 
% xlabel_locations = get(gca, 'XTick');
% set(gca, 'XTicklabel',  {taxis_all(xlabel_locations)})
% 
% set(gca, 'YDir','normal');

loadstat = 0;
outfile = 'allstat_clusteralpha005';
if loadstat;
    disp('Loading allstat . . .')
    load(fullfile(PREIN, outfile));
    return
end

corrtype = 'Spearman';

cfg = [];
cfg.channel          = {'MEG'};
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_correlationT';
cfg.type             = corrtype;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg.avgoverfreq      = 'no';
%     cfg.subj             = 1:28;
% prepare_neighbours determines what sensors may form clusters
cfg_neighb.method    = 'template'; %TODO do based on data (chans missing)
cfg_neighb.template  = 'CTF275_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb);
% ft_neighbourplot(cfg_neighb, freq)

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

allstat=[];
for idiff = 3 
    for idrug = 3:4 
        parfor itrig = 1:2
            freq = [];
            freq.label = chlabel;
            freq.freq = faxis_all;
            %         freq.time = taxis_all;
            freq.time = taxis{itrig};    % 1:length(taxis_all); % to fool
            freq.dimord = 'subj_chan_freq_time';
            
            %         freq.powspctrm = squeeze(respavg(:,:,:,tind{itrig}, idrug,3,idiff,itrig)); %DIMS: sub chan freq tind
            freq.powspctrm = squeeze(temp(:,:,:,:, idrug, idiff, itrig)); %DIMS: sub chan freq tind
            
            freqbehav = freq;
            for isub=1:nsub
                freqbehav.powspctrm(isub,:,:,:) = squeeze(driftrates(isub,idrug,idiff));
            end
            allstat{idrug, idiff, itrig} = ft_freqstatistics(cfg, freq, freqbehav);
            
        end
    end
end
save(fullfile(PREIN, outfile), 'allstat');


%% Collapse over freq or time, plot
THR = 0.05;
SAV = 1;
addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/cbrewer')
addpath('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/stats')

XLIM = [-0.25 1; -1 0.25];
% load('colormap170613.mat');
cmap_corr = cbrewer('div', 'RdBu',256);
cmap_corr = flipud(cmap_corr);

% clus = find(arrayfun(@(x) x.prob<=THR, allstat{end}.posclusters));
% cluslabels = allstat{end}.posclusterslabelmat;
% pvalues = arrayfun(@(x) x.prob, allstat{end}.posclusters(clus));

% rtemp = allstat{end}.rho;
% rtemp(cluslabels~=clus(1)) = 0;
% rtemp = squeeze(trapz(rtemp,2));%./numel(find(rtemp~=0));
% scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
itrig = 2;

trap_dim_leg = {'chan' 'freq' 'time'};
for iclus = 1:3%:3
    for trap_dim = 1 %1:3%:2
%         plotstat = allstat{end};
        plotstat = stat;
       
        temp = trapz(plotstat.rho, trap_dim);
        plotstat.rho = temp;
        plotstat.dimord = 'chan_freq_time';
        if trap_dim == 2 % collapse freqs
            plotstat.freq = 1;
        else
            plotstat.time = 1;
        end
        plotstat.posclusterslabelmat = any(plotstat.posclusterslabelmat == iclus, trap_dim); % select cluster 1, if any freq is part of the cluster, then include
        plotstat.negclusterslabelmat = zeros(size(plotstat.posclusterslabelmat));
        
        if trap_dim > 1
            % plotting with clusterplot
            % figure;    iplot=0; hold on
            % set(gcf, 'Position', [0 -200 375*3 210*4])
            %
            % load('colormap170613.mat');
            cmap_corr = cbrewer('div', 'RdBu',256);
            cmap_corr = flipud(cmap_corr);
            % subplot(223); colormap(cmap)
            
            close all
            cfg = [];
            cfg.alpha  = 0.025;
            cfg.parameter = 'rho';
            cfg.zlim   = [-15 15];
            cfg.layout = 'CTF275.lay';
            cfg.colorbar = 'no';
            cfg.subplotsize    = [5 5];
            cfg.highlightsizeseries = 1;
            %         cfg.highlightsymbolseries = '.';
            
            cfg.saveaspng = 'no';
            cfg.colormap = cmap_corr;
            cfg.maskparameter = 'posclusterslabelmat';
            cfg.style = 'straight';
            cd(PREOUT)
            
            ft_clusterplot(cfg, plotstat);
            
            colormap(cmap_corr);  %cmap = get(gcf, 'Colormap')
            set(gcf, 'Position', [0 -200 375*3 210*4])
            
            if SAV
                outpath = fullfile(PREOUT, 'topos');
                warning off; mkdir(outpath); warning on
                outfile = sprintf( '%s%stopoallstats_%slocked_%s_clus%d', outpath, filesep, trigger_leg{itrig}, trap_dim_leg{trap_dim}, iclus);
                disp(outfile)
                export_fig(outfile, '-pdf') %'-png', ,  '-depsc'  '-transparent'
            end;
            cd(outpath)
            
        else % TFR
%             close all
            figure;
            subplot(2,2,1)
            hold on
            %         load('colormap170613.mat');
            colormap(cmap_corr)
            %         colormap(cmap)
            
            ZLIM = [-90 90];
            dum = squeeze(plotstat.rho); %average over subj
            dum = showstatsTFR(dum, plotstat.posclusterslabelmat, 1:length(faxis_all), ZLIM, 1);
            
            imagesc(taxis{itrig}, faxis_all, dum,ZLIM);
%             YTICKS = [5:10:200];
%             YTICKS = [5:10:200];
            YTICKS = [10:10:150];
            set(gca,'Box','off','XTick',-1:0.5:2,...    [-0.5,0,0.5,1,1.5,2,2.5]
                'YDir','normal', 'YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                'TickDir','out', 'FontSize', 8);
            colorbar
            yaxis = [faxis_all(1), faxis_all(end)];
            plot([0,0],yaxis,'k',[0,0],yaxis,'k');
            xlim(XLIM(itrig,:))
            %             ylim([FREQLO 35])
            ylabel_locations = get(gca, 'YTick');
            freq_steps = round(length(faxis_all)/length(ylabel_locations));
            set(gca, 'YTicklabel',  {faxis_all(freq_steps:freq_steps:length(faxis_all))})

            
            ylabel('Frequency (Hz)');
            xlabel(sprintf('Time from %s (s)', trigger_leg{itrig}));
            title(sprintf('Cluster #%d, p = %g', iclus, plotstat.posclusters(iclus).prob))
            if SAV
                outpath = fullfile(PREOUT, 'topos');
                warning off; mkdir(outpath); warning on
                outfile = sprintf( '%s%sTFRallstats_%slocked_%s_clus%d', outpath, filesep, trigger_leg{itrig}, trap_dim_leg{trap_dim}, iclus);
                disp(outfile)
                export_fig(outfile, '-pdf') %'-png', ,  '-depsc'  '-transparent'
            end;
        end
    end
end
cd(outpath)




















% old


%     stats = {};
%     for ievent = 1:2
%         freq2stats = freq;
% %         freq2stats.freq = freq.freq(frind);
% %         freq2stats.dimord = 'subj_chan_freq_time';
%         freq2stats.powspctrm = squeeze((r(:,ievent,:,frind,:)));
%         freq2stats2 = freq2stats; %create zero freq to test against
%         freq2stats2.powspctrm = zeros(size(freq2stats.powspctrm));
%         stats{ievent} = ft_freqstatistics(cfg,freq2stats,freq2stats2);
%
%
%         cfg1 = [];
%         freq2stats = ft_freqdescriptives(cfg1,freq2stats);
%         freq2stats2= ft_freqdescriptives(cfg1,freq2stats2);
%
%         stats{ievent}.effect = freq2stats.powspctrm-freq2stats2.powspctrm;
%     end
%     save(statfile, 'stats');




% %%%old
% allstat=[];
% for idiff = 3 %:4 % 1:4 %[1,2,4] %1:4 %1:4
%     for idrug =  4 %1:4 %1:4 %1:4 %[1,2,4] %1:4
%         for itrig = 2 %1:2 %1:4 %1:4 %[1,2,4] %1:4
%             %             for iroi = 1:2 % beta 1 and 2
%             freq.time = taxis{itrig};
%             cfg.latency          = XLIM(itrig,:);
%
%             freq = [];
%             freq.dimord = 'subj_chan_freq_time';
%             freq.freq = faxis;
%             freq.time = taxis{itrig};
%             freq.label = chlabel;
%
%             freq.powspctrm = squeeze(respavg(:,:,:,tind{itrig}, idrug,3,idiff,itrig)); %DIMS: sub chan freq tind
%             freq.powspctrm = squeeze(temp(:,:,:,:, idrug, 3)); %DIMS: sub chan freq tind
%
%             freqbehav = freq; %create zero freq to test against
%             for isub=1:nsub
%                 freqbehav.powspctrm(isub,:,:,:) = squeeze(driftrates(isub,idrug,idiff));
%             end
%             allstat{idrug, idiff, itrig} = ft_freqstatistics(cfg, freq, freqbehav);
%
%             %             freqzero = freq;
%             %             freqzero.powspctrm = zeros(size(freq.powspctrm));
%             %
%             %             allstat{idrug, idiff, itrig} = ft_freqstatistics(cfg, freq, freqzero);
%             %             %             end
%         end
%     end
% end
% save(fullfile(PREIN, outfile), 'allstat');
