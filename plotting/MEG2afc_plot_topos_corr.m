%% Plot topo's of the modulation-driftrate correlations in beta1 and beta2

MEG2afc_load_respavg

driftrates = MEG2afc_load_driftrates161116(SUBJ, PREOUT);

%% Do stats on TF windows (rois)
tois = [-0.5 0; -0.9 -0.25];
fois = [13 17;   19 25];
corrtype = 'Spearman';

SAV= 1;

loadstat = 0;
outfile = 'topocorrstat';
if loadstat;
    disp('Loading topocorrstat . . .')
    load(fullfile(PREIN, outfile));
    return
end

cfg = [];
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_intersubcorr';
cfg.type             = corrtype;
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.1;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
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

topocorrstat=[];
for idiff = 3 %:4 % 1:4 %[1,2,4] %1:4 %1:4
    for idrug = 1:4 %1:4 %1:4 %[1,2,4] %1:4
        for itrig = 2 %1:2 %1:4 %1:4 %[1,2,4] %1:4
            for iroi = 1:2 % beta 1 and 2
                freq.time = taxis{itrig};
                cfg.latency          = XLIM(itrig,:);
                
                freq = [];
                freq.dimord = 'subj_chan_freq_time';
                freq.freq = faxis;
                freq.time = taxis{itrig};
                freq.label = chlabel;
                freq.powspctrm = squeeze(respavg(:,:,:,tind{itrig}, idrug,3,idiff,itrig)); %DIMS: sub chan freq tind
                
                cfg2=[];
                cfg2.avgovertime = 'yes';    cfg2.avgoverfreq = 'yes';
                cfg2.latency     = tois(iroi,:); % resp-locked beta2
                cfg2.frequency   = fois(iroi,:);
                
                freq = ft_selectdata(cfg2, freq);
                freq.dimord = 'subj_chan';
                
                freqbehav = freq; %create zero freq to test against
                for isub=1:nsub
                    freqbehav.powspctrm(isub,:) = squeeze(driftrates(isub,idrug,idiff));
                end
                topocorrstat{idrug, idiff, itrig, iroi} = ft_freqstatistics(cfg, freq, freqbehav);
            end
        end
    end
end
save(fullfile(PREIN, outfile), 'topocorrstat');


%% plotting with clusterplot
% % figure;    iplot=0; hold on
% % set(gcf, 'Position', [0 -200 375*3 210*4])
% %
% close all
% load('colormap170613.mat');
% cfg = [];
% cfg.alpha  = 0.025;
% cfg.parameter = 'rho';
% cfg.zlim   = [-1 1];
% cfg.layout = 'CTF275.lay';
% cfg.colorbar = 'yes';
% cfg.subplotsize    = [1 1];
% % figure
% for idiff = 3 %:4 % 1:4 %[1,2,4] %1:4 %1:4
%     for idrug = [1,2,4] %1:4
%         for itrig = 2 %1:2 %1:4 %1:4 %[1,2,4] %1:4
% %             subplot(2,3,idrug); hold on
%             test = ft_clusterplot(cfg, topocorrstat{idrug, idiff, itrig});
%             % title('henkie', 'FontSize', 20)
%             colormap(cmap);  %cmap = get(gcf, 'Colormap')
%         end
%     end
% end

%% plotting
close all

SAV = 1;
freq.dimord = 'chan';
NKsensorselection

cfg = [];
% cfg.layout = '/home/mpib/kloosterman/MATLAB/toolbox/fieldtrip-20150803/template/CTF275.lay';       % neuromag306cmb neuromag306all neuromag306mag neuromag306planar
cfg.layout = 'CTF275.lay';       % neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat
% cfg.layout = 'CTF275_helmet.mat';       % neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat
cfg.comment = 'no'; %sprintf('zlim: [%g %g]', ZLIM(1),ZLIM(2)); %no xlim
cfg.marker = 'on';
cfg.shading = 'flat';
cfg.style = 'straight'; %both
% cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfg.markersymbol = '.';
cfg.markersize = 1;
cfg.colorbar = 'no';

imotor = 3;

figure;    iplot=0; hold on
set(gcf, 'Position', [0 -200 375*3 210*4])

load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
for idiff= 3 % [1,2,4]
    for iroi = 1:2 % beta 1 and 2
        for idrug = [1,2,4] % 1:4 % 1:2
            for itrig = 2 %1:2
                
                %             ZLIM    = [-1 1];
                ZLIM    = [-0.7 0.7];
                
                iplot = iplot+1; subplot(2,3,iplot); hold on
                
                freq.powspctrm = topocorrstat{idrug, idiff, itrig, iroi}.rho;
                
                %             cfg.xlim = tois;
                %             cfg.ylim = bandoi;
                cfg.zlim = ZLIM;
                
                cfg.highlight          = 'on';
                cfg.highlightsymbol    = '*';
                cfg.markersize         = 0.5;
                %                     cfg.markersize         = 5;
                cfg.highlightfontsize = 12;
                % cfg.markercolor        = [1 0 0];
                %                     cfg.highlightcolor        = [1 0 0];
                cfg.highlightcolor        = [0 0 0];
                cfg.highlightsize      = 6;
                %                     cfg.highlightchannel = occind;
                try
                    cfg.highlightchannel = find(topocorrstat{idrug, idiff, itrig, iroi}.posclusterslabelmat > 0);
                catch
                    cfg.highlightchannel = [];
                end
                
                test = ft_topoplotTFR(cfg,freq);
                
                title({sprintf('%s %s''s r= [%g to %g]',  pharm_conds{idrug}, corrtype,  ZLIM) ; % motor_conds{imotor},  diff_conds{idiff},
                    sprintf('%g-%g s wrt %s, %d-%d Hz', tois(iroi,:), trigger_leg{itrig}, fois(iroi,:)  )}, 'Fontsize', 14); % 'FontWeight','bold'
                
            end
        end
    end
end
if SAV
    outpath = fullfile(PREOUT, 'topos');
    warning off; mkdir(outpath); warning on
    
%     outfile = sprintf( '%s%stopostats_%slocked_%s_freq_%g-%gms_%d-%dHz', outpath, filesep, trigger_leg{itrig}, ...
%         motor_conds{imotor}, tois*1000, bandoi  );
    outfile = sprintf( '%s%stopostats_%slocked', outpath, filesep, trigger_leg{itrig});
    disp(outfile)
    
    %             orient portrait
    orient landscape
    %             print('-dpdf', outfile)
    print('-dpng', outfile)
                export_fig(outfile, '-pdf') %'-png', ,  '-depsc'  '-transparent'
    
end;
cd(outpath)
