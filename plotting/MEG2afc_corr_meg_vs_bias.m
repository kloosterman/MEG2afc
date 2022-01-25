% % addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
% % addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
% % ft_defaults

trigger = 'resp'

MEG2afc_load_respavg

%% load bias, split based on bias type and plot bars
close all
SAV=1;
% JW's 0.5 based proportion right responses
% biasdat =  csvread('drug_bias.csv', 1,1); % 19 subj, col 1 = drug

% % c_drug c_plac d_drug d_plac
% biasdat = csvread('drug_bias-critdprime.csv', 1,1); % 17 subj, col 1 = drug
% % biasdat = biasdat(:,1:2); % c_drug c_plac
% biasdat = biasdat(:,3:4); % d_drug d_plac

% dc_drug dc_drug_easy dc_drug_hard dc_plac dc_plac_easy dc_plac_hard
biasdat = csvread('drug_drift_criterion.csv', 1,1); % 17 subj, col 1 = drug
% biasdat = biasdat(:,1:2); % c_drug c_plac
% biasdat = biasdat(:,[1,4]); % d_drug d_plac
biasdat = biasdat(:,[2,5]); % d_drug d_plac


biasdat = fliplr(biasdat); % ALARM! now col 1 = plac

split = (biasdat(:,1) > 0) + 1; % split wrt bias direction: 12 vs 7 subj

%  %split based on plac data.
% split = (biasdat(:,1) > 0.5) + 1; % split wrt bias direction: 12 vs 7 subj
% split = (biasdat(:,1) > median(biasdat(:,1))) + 1; %10 have a left bias, 9 a right; JW split met median! % split = right_button_plac < np.median(right_button_plac)

bias_leg = { 'Left', 'Right' };
figure
set(gcf, 'Position', [0 -200 375*3 210*4])

for isplit = 1:2
    subplot(1,3,isplit)
    dum = biasdat(isplit == split, :);
    dumstd = std(dum) / sqrt(length(dum));
    
    barweb(mean(dum), dumstd, 0.75, [], [], [], 'Fraction right', [ 0 0 1; 1 0 0]) %, pharm_conds(1:idrug)
    hold on; axis square

    plot(get(gca, 'XLim'), [0.5 0.5], '--k')
%     ylim([0.45 0.55])
    ylim([-0.1 0.1])
    legend(pharm_conds{2:-1:1})
    
%     pval = randtest1d(dum(:,1), dum(:,2), 0, 10000);
                [~,pval] = ttest(dum(:,1), dum(:,2));

    title(sprintf('%s-hand bias, N = %d, p = %g', bias_leg{isplit}, length(dum), pval))
end
subplot(1,3,3)

scatter(biasdat(:,1), biasdat(:,2) - biasdat(:,1));
hold on
[r,p] = corr(biasdat(:,1), biasdat(:,2) - biasdat(:,1));
title(sprintf('r = %g p = %g', r, p));
lsline
axis square
% xlim([0.44 0.56])
%     xlim([-0.2 0.2])

% plot( [0.5 0.5], get(gca, 'YLim'), '--k')
xlabel('Button bias during placebo')
% text(0.46, -0.08, 'LEFT')
% text(0.52, -0.08, 'RIGHT')

ylabel('Bias change drug vs. placebo')

%to correlate drug - plac MEG vs bias
% calculate abs change 
bias_change_abs = abs( biasdat(:,2) - biasdat(:,1))

% bias_change_abs = abs( biasdat(:,2)-0.5 ) - abs(biasdat(:,1)-0.5 ); % weird

% calculate  change 
% bias_change_abs =  biasdat(:,2) - biasdat(:,1)

% % TODO NK2, 
% bias_change_abs = bias_change_abs([1,3:end],:)


if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
%     outfile = sprintf( '%s%sResponse_bias', outpath, filesep ); % motor_conds{imotor}, pupilconds{ipup}
    outfile = sprintf( '%s%sCriterion', outpath, filesep ); % motor_conds{imotor}, pupilconds{ipup}
    disp(outfile)
            export_fig(outfile, '-pdf', '-painters', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
    %         export_fig(outfile, '-eps') %'-png',  '-pdf',
%     print(outfile, '-dpsc', '-painters') %'-png',  '-pdf',
end;

%% do stats on correlation TFR's and bias
% load driftrates
% driftrates(:,4,:) = driftrates(:,1,:) - driftrates(:,2,:);
% driftrates(:,:,4) = driftrates(:,:,1) - driftrates(:,:,2);

loadstat = 0;
outfile = ['corrstat_' trigger];
if loadstat;
    disp('Loading corrstat . . .')
    load(fullfile(PREIN, 'respavg', outfile));
    return
end

freq=[];
freq.time = taxis;
freq.freq = faxis;
freq.dimord = 'chan_subj_freq_time';
freq.label = {'custompooling'};

cfg = [];
cfg.latency          = XLIM;
cfg.frequency        = YLIM;
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_intersubcorr'; %%% corr statfun!
cfg.correctm         = 'cluster';
% % % %     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 100;
cfg.neighbours       = []; %in case no channel data present

cfg.type = 'Pearson'

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


corrstat=[];
for isoi = 1:2 %1:3%:2 %:2 % occ and motor
    for idiff = 3 % 1:2 %[1,2,4] %1:4 %1:4
        for idrug = 4 %[1,2,4] %1:4
            
            freq.powspctrm = squeeze(respavg(:,sens.ind{isoi},:,:, idrug,3,idiff,3,3)); %DIMS:
            if isoi < 3
                freq.powspctrm = squeeze(mean(freq.powspctrm,2)); % avg across sensors
            end
            freq.powspctrm = shiftdim(freq.powspctrm, -1);
            
            freqbehav = freq; %create zero freq to test against
            for isub=1:nsub
                %                     freqbehav.powspctrm(isub,:,:) = squeeze(driftrates(isub,idrug,idiff));
                freqbehav.powspctrm(1,isub,:,:) = squeeze(bias_change_abs(isub));
            end
            corrstat{isoi, idrug, idiff} = ft_freqstatistics(cfg, freq, freqbehav);
        end
    end
end
save(fullfile(PREIN, 'respavg', outfile), 'corrstat');


%% Correlate TFR's with abs bias and plot
% close all
SAV=1;
% % loadlat = 1;
% % outfile = ['respavglat_' trigger];
% % if loadlat;
% %     disp('Loading respavglat_ . . .')
% %     load(fullfile(PREIN, 'respavg', outfile));
% %     
% % end
% 
close all
SAV=1;
set(0,'DefaultAxesFontsize',12)

showstats = 0;
imotor = 3;
istim = 3;
iresp = 3;

TYP = '2afc';
figh=figure;    iplot=0; hold on;
set(gcf, 'Position', [0 -200 375*4 210*5])
load('colormap170613.mat');
colormap(cmap);  %cmap = get(gcf, 'Colormap')
for isoi = [2, 7:9] 
% for isoi = [1,4:6] 
    for idiff= 3 % 1:3 %[1:2,4]
        for idrug = [1:2,4] % 1:2

            ZLIM_pow = [-0.025 0.025];
            if idrug == 4
                ZLIM_pow = [-0.02 0.02];
            end
                
            ZLIM_corr = [-1 1];
            iplot = iplot+1; subplot(4,4,iplot); hold on
            dum = squeeze(respavg(:,sens.ind{isoi},frind,:,idrug, imotor, idiff, istim, iresp)); %average over subj
            if isoi < 3
                dum = squeeze(mean(dum,2)); % average over sensors
            else
                dum = squeeze(dum); % just sq when 1 sensor
            end
            
            rcoeff=[];
            for ifreq = 1:size(dum,2)
                for itind = 1:size(dum,3)
                    rcoeff(ifreq,itind) = corr( bias_change_abs , squeeze(dum(:,ifreq,itind)), 'type', 'Pearson' );
                end
            end
            if showstats
                rcoeff = showstatsTFR(rcoeff, squeeze(corrstat{isoi, idrug, idiff}.prob), frind, ZLIM_pow, showstats);
            end
            
            imagesc(taxis, faxis, squeeze(mean(dum)), ZLIM_pow);

            title({sprintf('power; %s %s [%g]', sens.leg{isoi}, trigger, ZLIM_pow(2)) ;
                sprintf('%s', pharm_conds{idrug})}); % 'FontWeight','bold'

            xlim(XLIM)
            ylim([FREQLO 150])
            %             colorbar
            yaxis = [FREQLO,FREQHI];
            
            plot([0,0],yaxis,'k',[0,0],yaxis,'k');
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            
            if idrug == 4 %plot correlation TFR
                iplot = iplot+1; subplot(4,4,iplot); hold on
                imagesc(taxis, faxis, rcoeff, ZLIM_corr);
                
                title({sprintf('power-abs.bias corr; %s [%g]', sens.leg{isoi}, ZLIM_corr(2)) ;
                    sprintf('%s  %s',  diff_conds{idiff},  pharm_conds{idrug})}); % 'FontWeight','bold'
                hold on
                plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                
                xlim(XLIM)
                ylim([FREQLO 150])
                YTICKS = [0:50:200];
                set(gca,'Box','off','XTick',-2:0.2:2,...    [-0.5,0,0.5,1,1.5,2,2.5]
                    'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                    'TickDir','out', 'FontSize', 12);
                LAB=1;
                if LAB %&& (iplot == 9 || iplot == 19)
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    zlabel('Response (%)');
%                     h=colorbar;
                end
            end
        end
    end
end

if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%sTFR_corr_criterion_%slocked_%s_%sfreq_%s', outpath, filesep, trigger, TYP, analysistype, sens.leg{isoi} ); % motor_conds{imotor}, pupilconds{ipup}
    disp(outfile)
    export_fig(outfile, '-pdf', '-painters', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
    %         export_fig(outfile, '-eps') %'-png',  '-pdf',
    %         print(outfile, '-dpsc', '-painters') %'-png',  '-pdf',
end;

%% Correlate drift rate with significant TF clusters
% 
close all
SAV=1;
freq=[];
freq.time = taxis;
freq.freq = faxis;
freq.dimord = 'subj_chan_freq_time';
freq.label = chlabel;

% tois = [-0.4 0.2; -0.4 -0.2]; %mod then lat
% fois = [0 40; 30 50];
tois = [-0.2 0]; %mod then lat
fois = [20 30];
cfg=[];
cfg.trials      = 'all';         
cfg.avgovertime = 'yes';      
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'yes';

corrstat=[];
figure; iplot=0;
set(gcf, 'Position', [0 -100 375*3 210*4])
for isoi = 4 
    freq.label = chlabel(sens.ind{isoi});
%     for imeas = 1%:2 % modulation or lateralization
        for idiff = 3 %1:3 %[1,2,4] %1:4 %1:4
            for idrug = 4 % [1,2,4] %1:4
                
                freq.powspctrm = squeeze(respavg(:,sens.ind{isoi},:,:, idrug,3,idiff,3,3)); %DIMS:
                if isoi == 4 % make dummy chan dimension
%                     freq.powspctrm = squeeze(respavg(:,sens.ind{isoi},:,:, idrug,3,idiff,3,4)); %DIMS:
                    pw_size = size(freq.powspctrm);
                    temp = nan([pw_size(1) 1 pw_size(2:3)] );
                    temp(:,1,:,:) = freq.powspctrm;
                    freq.powspctrm = temp;
                end

                cfg.latency     = tois(1,:);
                cfg.frequency      = fois(1,:);
                
                freqavg = ft_selectdata(cfg, freq);
%                 [r, p] = corr(freqavg.powspctrm, squeeze(driftrates(:,idrug, idiff)));
%                 [p, r] = permuteCorr(freqavg.powspctrm, squeeze(driftrates(:,idrug, idiff)));
                [p, r] = permuteCorr(freqavg.powspctrm, bias_change_abs, 10000, 0);
                
                iplot=iplot+1; subplot(2,3,iplot); hold on; axis square; box on
                for isub=1:nsub
                    plot(freqavg.powspctrm(isub), bias_change_abs(isub), 'o', 'MarkerSize', 12, ...
                        'Color', 'white', 'markerfacecolor', [0 0 0], 'Linewidth', 1);
                end
%                 scatter(freqavg.powspctrm, squeeze(driftrates(:,idrug, idiff)))
                if p < 0.05
                    %                     lsline
                    beta=polyfit(freqavg.powspctrm, squeeze(bias_change_abs),1);
                    h_ls = refline(beta);
                    set(h_ls, 'Linewidth', 1.5, 'Color', [0 0 0])
                end
                title({sprintf('%s, r = %g, p = %g', diff_conds{idiff}, r,p), ...
                    sprintf('%d-%d Hz, %g-%g s', cfg.frequency, cfg.latency)})
                xlabel(sprintf('%s, %s', sens.leg{isoi}, pharm_conds{idrug}))
                ylabel(sprintf('Abs. Bias, %s', pharm_conds{idrug}))

               
                set(gca, 'Tickdir', 'out')
%                 freqbehav = freq; %create zero freq to test against
%                 for isub=1:nsub
%                     freqbehav.powspctrm(isub,:,:) = squeeze(driftrates(isub,idrug,idiff));
%                 end
%                 corrstat{isoi, idrug, idiff} = ft_freqstatistics(cfg, freq, freqbehav);
            end
        end
%     end
end
if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%sscatter_corr_driftrate_%slocked_%s_%sfreq', outpath, filesep, trigger, TYP, analysistype ); % motor_conds{imotor}, pupilconds{ipup}
    disp(outfile)
    export_fig(outfile, '-pdf', '-painters', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
%             export_fig(outfile, '-eps') %'-png',  '-pdf',
%     print(outfile, '-dpsc', '-painters') %'-png',  '-pdf',
end;



