% MEGsplitCorrelation_per_cond_Hipp_method

% This script computes correlation between single trial measures and MEG.
% It does so separately for every subject, channel, frequency and time bin.
% Afterwards the correlations are averaged over relevant dimentions and
% tested for significance.
% It takes data from the script
% plotTFR_SurpriseReplayexp_Bestmann_surprise_split_per_cond as input.

SAV = 1;

overwrite = 1;
overwrite_stat = 1;
proj_outRT = 0;
regressOut = 1;
baseline = 1;


corrfactor      =  {'pupweights','surpriseGamma','entropyGamma','KLdivGamma','prev_trl','RT',...
    'surpriseSB', 'entorpySB', 'KLdivSB', 'Random', 'pupbasedat', 'Hazard', 'LP', 'TDsurprise', 'TDentropy'};

triggers = {'resp', 'stim'};
analysistypes = {'low', 'high'};


for iana = 1:2
    for itrg = 2:-1:1%1:2
        trigger = triggers{itrg};
        analysistype = analysistypes{iana};
        regime          = 'all';
        
        cond = {'all'};%{'FlatH', 'ShGau', 'LoGau'};
        event = {'Off', 'On', 'On-Off'};
        sign = {'Pos', 'Neg'};
        
        SUBJ = {
            'SB' 'SR' 'SN' 'MR' 'BP' 'JW' 'KH' 'CR' 'LR' 'LK' 'FH' 'JP' 'UD' 'TW' 'MW' 'CS' 'IP' 'MS' 'AT' 'DH'...
            'LB', 'VS', 'AE', 'FV', 'AR', 'MSa', 'MB', 'AK'};
        
        soiUnilateral_cmb_CABMSI_optimal_occ % get the sensors of interest (soi) list
        cd(sprintf('~/data/MEG/freq/%s/SB/160113/%s/', analysistype, trigger))
        examplefreq = dir('*200*totalpow_freq.mat');
        load(examplefreq(2).name);
        chlabel = freq.label;

        outdir = sprintf('~/data/MEG/freq/%s/MEGsplit/%s', analysistype, trigger);
        cd(outdir);
        %%
        for icorr = [2 3 1 11];%[1 2 3 13 14];

            outdir = sprintf('~/data/MEG/freq/%s/MEGsplit/%s', analysistype, trigger);
            cd(outdir);
            for isub = 1:length(SUBJ)
                subr = nan(4,2,length(freq.label), length(freq.freq), length(freq.time));
                subrespavg = subr;
                
                if regressOut %&& (icorr == 14 || icorr == 1)
                    subfile = sprintf('%s_corr_respavg_%s_per_cond_baseline%d_noRT%d_noSup.mat', SUBJ{isub}, corrfactor{icorr}, baseline,proj_outRT);
                else
                    subfile = sprintf('%s_corr_respavg_%s_per_cond_baseline%d_noRT%d.mat', SUBJ{isub}, corrfactor{icorr}, baseline,proj_outRT);
                end
                startfile = [outdir,filesep SUBJ{isub}, '_',analysistype,'_',trigger, '_',corrfactor{icorr},'_started.mat'];
                
                if (~exist(subfile)  &&  ~exist(startfile)) || (overwrite && ~exist(startfile))
                    a = [];
                    save(startfile,'a');
                    fprintf('%s\n',startfile)
                    
                    tic
                    for ievent = 1:2
                        for itype = 1:3
                            fprintf('Running %s, event %d\n', SUBJ{isub}, ievent)
                            file2load = dir(sprintf('%s_*%s_stb%d_type%d_event%d.mat',SUBJ{isub},cond{1}, baseline,itype, ievent));
                            if ~isempty(file2load)
                                load(file2load.name);
                                
                                for i = 1:length(correspdat.labels)
                                    if strcmp(correspdat.labels{i}, corrfactor{icorr})
                                        indx = i;
                                    end
                                end
                                
                                data = correspdat.data;
                                
                                
                                b = squeeze(correspdat.splitmeasures(:,indx));
                                ind = find(b~=0);
                                if  proj_outRT
                                    normRT = squeeze(correspdat.splitmeasures(:,6))/norm(squeeze(correspdat.splitmeasures(:,indx)));
                                    normRT = normRT(ind);
                                    X = [ones(length(b(ind)),1),normRT];
                                    [beta,~,res] = regress(b(ind),X);
%                                     res = res+beta(1);
                                    b(ind) = res;
                                end
                                
                                if regressOut && ~isempty(ind) %&& (icorr == 1 || icorr==14)
%                                     if icorr == 14
%                                         normVar = squeeze(correspdat.splitmeasures(:,15))/norm(squeeze(correspdat.splitmeasures(:,indx)));
%                                     elseif icorr == 1
%                                         normVar = squeeze(correspdat.splitmeasures(:,11))/norm(squeeze(correspdat.splitmeasures(:,indx)));
%                                     end
%                                     normVar = normVar(ind);
%                                     X = [ones(length(b(ind)),1),normVar];
%                                     [beta,~,res] = regress(b(ind),X);
%                                     res = res+beta(1);
%                                     b(ind) = res;
                                    if icorr == 14; 
                                        prj = 15; 
                                    elseif icorr == 15;
                                        prj = 14;
                                    elseif icorr == 1; 
                                        prj = 11; 
                                    elseif icorr == 11;
                                        prj = 1;
                                    elseif icorr == 2;
                                        prj = 3;
                                    elseif icorr == 3;
                                        prj = 2;
                                    end
                                    
                                    
                                    b(ind) = projectout(b(ind),correspdat.splitmeasures(ind,prj));
                                end
                                
                                %                                     if icorr == 3
                                %                                         b(2:end) = b(1:end-1);
                                %                                         b(1) = 0;
                                %                                     end
                                
                                tic
                                for i = 1:size(data,2)
                                    for j = 1:size(data,3)
                                        for k = 1:size(data,4)
                                            a = squeeze(data(:,i,j,k));
                                            ind = find(~isnan(a)&~isnan(b)&a~=0&b~=0);
                                            a = a(ind);
                                            bind = b(ind);
                                            if ~isempty(ind)
                                                if length(ind) >2
                                                    [subr(itype,ievent,i,j,k), ~] = corr(a, bind);
                                                end
                                            end
                                        end
                                    end
                                    
                                end
                                subrespavg(itype,ievent,:,:,:) = squeeze(nanmean(data));
                            end
                        end
                    end
                    save(subfile, 'subr', 'subrespavg');
                    time_passed = toc;
                    fprintf('Computation time: %1.2f minutes.\n', time_passed/60);
                end
            end
            
            if isCommandWindowOpen<2
                
                r = NaN(length(SUBJ), 4,2, length(freq.label), length(freq.freq), length(freq.time));
                p = r;
                respavg = r;
                cd(outdir)
                
                list = dir(sprintf('*_%s_%s_%s_started.mat', analysistype, trigger, corrfactor{icorr}));
                sprintf('Deleting %d _started-files\n', length(list))
                for i = 1:length(list)
                    delete(list(i).name);
                end
                
                for isub = 1:length(SUBJ)
                    if regressOut && (icorr == 14 || icorr == 1)
                        subfile = sprintf('%s_corr_respavg_%s_per_cond_baseline%d_noRT%d_noSup.mat', SUBJ{isub}, corrfactor{icorr}, baseline,proj_outRT);
                    else
                        subfile = sprintf('%s_corr_respavg_%s_per_cond_baseline%d_noRT%d.mat', SUBJ{isub}, corrfactor{icorr}, baseline,proj_outRT);
                    end
                    fprintf('Loading %s\n',subfile)
                    load(subfile)
                    r(isub,1:size(subr,1),:,:,:,:) = subr;
                    respavg(isub,1:size(subrespavg,1),:,:,:,:) = subrespavg;
                end
                respavg = squeeze(nanmean(respavg,2));
                r = squeeze(nanmean(r,2));
                
                %% Start of stats
                
                
                addpath('/mnt/homes/home024/meindertsmat/Documents/MATLAB/fieldtrip-20130305/');
                rmpath('/mnt/homes/home024/meindertsmat/Documents/MATLAB/fieldtrip-20150811/');
                ft_defaults;
                
                switch analysistype
                    case 'low'
                        foi = [3 35];
                    case 'high'
                        foi = [60 140];
                end
                frind = find(freq.freq>=foi(1) & freq.freq<=foi(2));
                
                statfile = sprintf('%s/Stats3D_%s_%s_%s_regress%d.mat',outdir, corrfactor{icorr}, analysistypes{iana}, triggers{itrg}, regressOut);
                if ~exist(statfile,'file') || overwrite_stat
                    cfg0 = [];
                    cfg0.channel          = {'MEG'};
                    cfg0.latency          = 'all';
                    cfg0.frequency        = foi;
                    cfg0.method           = 'montecarlo';
                    cfg0.statistic        = 'depsamplesT';
                    cfg0.correctm         = 'cluster';
                    cfg0.clusteralpha     = 0.05;
                    cfg0.clusterstatistic = 'maxsum';
                    cfg0.minnbchan        = 2;
                    cfg0.tail             = 0;
                    cfg0.clustertail      = 0;
                    cfg0.alpha            = 0.025;
                    cfg0.numrandomization = 1000;
                    cfg0.avgoverfreq      = 'no';
                    cfg0.subj             = 1:28;
                    cfg_pn               = [];
                    cfg_pn.method        = 'template';
                    cfg0.neighbours       = ft_prepare_neighbours(cfg_pn, freq);
                    
                    nsub = size(r,1);
                    design = zeros(2,2*nsub);
                    for i = 1:nsub
                        design(1,i) = i;
                    end
                    for i = 1:nsub
                        design(1,nsub+i) = i;
                    end
                    design(2,1:nsub)        = 1;
                    design(2,nsub+1:2*nsub) = 2;
                    
                    cfg0.design   = design;
                    cfg0.uvar     = 1;
                    cfg0.ivar     = 2;
                    
                    stats = {};
                    for ievent = 1:2
                        freq2stats = freq;
                        freq2stats.freq = freq.freq(frind);
                        freq2stats.dimord = 'subj_chan_freq_time';
                        freq2stats.powspctrm = squeeze((r(:,ievent,:,frind,:)));
                        freq2stats2 = freq2stats; %create zero freq to test against
                        freq2stats2.powspctrm = zeros(size(freq2stats.powspctrm));
                        stats{ievent} = ft_freqstatistics(cfg0,freq2stats,freq2stats2);
                        
                        
                        cfg1 = [];
                        freq2stats = ft_freqdescriptives(cfg1,freq2stats);
                        freq2stats2= ft_freqdescriptives(cfg1,freq2stats2);
                        
                        stats{ievent}.effect = freq2stats.powspctrm-freq2stats2.powspctrm;
                    end
                    save(statfile, 'stats');
                else
                    load(statfile);
                end
                
                
                %%
                THR = 0.05;
                trap = 1;

                for ievent = 1:2
                    for isign =1:2
                        if isign == 1
                            clus = find(arrayfun(@(x) x.prob<=THR, stats{ievent}.posclusters));
                            cluslabels = stats{ievent}.posclusterslabelmat;
                            pvalues = arrayfun(@(x) x.prob, stats{ievent}.posclusters(clus));
                        else
                            clus = find(arrayfun(@(x) x.prob<=THR, stats{ievent}.negclusters));
                            cluslabels = stats{ievent}.negclusterslabelmat;
                            pvalues = arrayfun(@(x) x.prob, stats{ievent}.negclusters(clus));
                        end
                        
                        if ~isempty(clus)
                            for iclus = clus
                                figure
                                
                                % for TFR
                                rtemp = squeeze(nanmean(r(:,ievent,:,frind,:)));
                                if trap
                                    
                                    rtemp(cluslabels~=clus(iclus)) = 0;
                                    rtemp = squeeze(trapz(rtemp,1));%./numel(find(rtemp~=0));
                                    scale = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                                else
                                    rtemp(cluslabels~=clus(iclus)) = NaN;
                                    rtemp = squeeze(nanmean(rtemp,1));%./numel(find(rtemp~=0));
                                    rtemp(isnan(rtemp))=0;
                                    scale = [-0.05 0.05];
                                end
                                cmap = cbrewer('div', 'RdBu',256);
                                cmap = flipud(cmap);
                                subplot(223); colormap(cmap)
                                %                                         rtemp = smooth2a(rtemp,1);
                                ft_plot_matrix(freq.time,freq.freq(frind), rtemp,'clim',scale);%, 'highlight', thrcluster, 'highlightstyle', 'saturation', 'clim', scale) %,taxis,faxis
                                hold on
                                set(gca,'YDir','normal');
                                plot([0,0],foi,'k',[0,0],foi,'k');
                                ylabel('Frequency (Hz)');
                                xlabel(sprintf('Time from %s (s)', trigger));
                                
                                
                                rtemp_time = squeeze(trapz(rtemp,1));
                                rtemp_freq = squeeze(trapz(rtemp,2));
                                
                                subplot(221);
                                area(freq.time, rtemp_time); hold on;
                                plot([freq.time(1), freq.time(end)], [0 0],'k');
                                xlim([freq.time(1), freq.time(end)]);
                                ylim([-(ceil(max(abs(rtemp_time)))).*1.2, ceil(max(abs(rtemp_time))).*1.2]);
                                xlabel(sprintf('Time from %s (s)', trigger));
                                ylabel('Integrated correlation');
                                title(sprintf('%s Corr %s %s-%d\n[%d %d] p=%1.4f', corrfactor{icorr}, event{ievent}, sign{isign}, iclus, scale(1), scale(2), pvalues(iclus)));
                               
                                 
                                subplot(224)
                                area(freq.freq(frind), rtemp_freq); hold on;
                                plot([freq.freq(frind(1)) freq.freq(frind(end))],[0 0],'k');
                                xlim([freq.freq(frind(1)) freq.freq(frind(end))]);
                                ylim([-(ceil(max(abs(rtemp_freq)))).*1.2, ceil(max(abs(rtemp_freq))).*1.2]);
                                camroll(270)
                                camorbit(0,180)
                                xlabel('Frequency (Hz)');
                                ylabel('Integrated correlation');
                                
                                % for topo
                                % CFG
                                subplot(222); colormap(cmap);
                                cfg = [];
                                cfg.layout = 'CTF275.lay';
                                cfg.comment = 'no';
                                cfg.marker = 'off';
                                cfg.shading = 'flat';
                                cfg.style = 'straight'; %both
                                cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
                                cfg.markersize = 1;
                                
                                rtemp = squeeze(nanmean(r(:,ievent,:,frind,:)));
                                if trap
                                    rtemp(cluslabels~=clus(iclus)) = 0;
                                    rtemp = squeeze(trapz(rtemp,2));
                                    rtemp = squeeze(trapz(rtemp,2));
                                    cfg.zlim = [-(ceil(max(abs(rtemp(:))))), ceil(max(abs(rtemp(:))))];
                                else
                                    rtemp(cluslabels~=clus(iclus)) = NaN;
                                    rtemp = squeeze(nanmean(rtemp,2));
                                    rtemp = squeeze(nanmean(rtemp,2));
                                    rtemp(isnan(rtemp))=0;
                                    cfg.zlim = [-0.25 0.25];
                                end
                                
                                freq2 = freq;
                                freq2.dimord = 'chan_freq_time';
                                freq2.powspctrm = repmat(rtemp,[1,2,2]);
                                freq2.time = [1 2];
                                freq2.freq = [1 2];
                                warning off;
                                
                                ft_topoplotTFR(cfg,freq2);
                                
                                if SAV
                                    outdir = sprintf('~/data/MEG/freq/%s/plots/MEGsplit/',analysistype);
                                    outfile = sprintf('%sHipp_corr_%s-locked_%s_%s_%s%d', outdir, trigger, corrfactor{icorr}, event{ievent}, sign{isign}, clus(iclus));
                                    display(outfile)
                                    print('-dpdf',outfile)
                                end
                            end
                        end
                    end
                end
                
            end
            
        end
    end
end
