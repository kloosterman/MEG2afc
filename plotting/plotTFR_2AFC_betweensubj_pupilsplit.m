%% Split subjects into low and high pupil group, plot TFR's for each condition as normal
% First run plotTFR_2AFC_poolings FIRST CELL to get respavg, varnames etc.
% Compute respavgpool here with 2 levels (pupil lo and hi)

% load jw's pupil data
pupil = load('/mnt/homes/home022/nkloost1/pupildata/pupil_data.mat');

% % split based on pupil placebo decision
% [~, rank_ind] = sort(pupil.pupil_decision_placebo([1, 3:end])); % drop subj 2, incomplete
% pupil_split_lohi = {rank_ind(1:9) rank_ind(10:end)}; % 9 sub in lo, 9 in hi

% split based on pupil drug-placebo decision
pupil_diff =  pupil.pupil_decision_drug - pupil.pupil_decision_placebo;
[~, rank_ind] = sort(pupil_diff([1, 3:end]), 'descend'); % drop subj 2, incomplete
pupil_split_lohi = {rank_ind(1:9) rank_ind(10:end)}; % 9 sub in lo, 9 in hi

% respavgsoi = nan([ 18   2    73    length(tind)     4     4     4     4     4]); % collapse sens dim
% NKsensorselection
% respavgsoi(:,1,:,:,1:2, 1:2,1:3,1:3,1:3) = squeeze(mean(respavg(:,occind,:,:,:, :,:,:,: ), 2 )); %average over sens
% respavgsoi(:,2,:,:,1:2, 1:2,1:3,1:3,1:3) = squeeze(mean(respavg(:,motorind,:,:,:, :,:,:,: ), 2 )); %average over sens
% respavgsoi(:,:,:,:,:,3,:,:,:) = mean(respavgsoi(:,:,:,:,:,1:2,:,:,:), 6); %avg over motor regime
% respavgsoi(:,:,:,:,3,:,:,:,:) = mean(respavgsoi(:,:,:,:,1:2,:,:,:,:), 5); %avg over pharma
% respavgsoi(:,:,:,:,4,:,:,:,:) = respavgsoi(:,:,:,:,1,:,:,:,:) - respavgsoi(:,:,:,:,2,:,:,:,:) ; %pharma diff
% respavgsoi(:,:,:,:,:,4,:,:,:) = respavgsoi(:,:,:,:,:,1,:,:,:) - respavgsoi(:,:,:,:,:,2,:,:,:) ; %motor ipsi - contra
% respavgsoi(:,:,:,:,:,:,4,:,:) = respavgsoi(:,:,:,:,:,:,1,:,:) - respavgsoi(:,:,:,:,:,:,2,:,:) ; %easy-hard
% respavgsoi(:,:,:,:,:,:,:,4,:) = respavgsoi(:,:,:,:,:,:,:,1,:) - respavgsoi(:,:,:,:,:,:,:,2,:) ; %stim left-right
% respavgsoi(:,:,:,:,:,:,:,:,4) = respavgsoi(:,:,:,:,:,:,:,:,1) - respavgsoi(:,:,:,:,:,:,:,:,2) ; %resp left-right
% 
% respavgpool = nan([ 2   2    73    61     4  4   4   4   4]); % replace subj dim w pupil dim
% disp('matrices preallocated')

%% TODO stats
FREQLO = 5;
FREQHI = 149;
if strcmp(trigger, 'stim')
    TIMLO          = -0.2;
    TIMHI          =  0.4;
else
    TIMLO          = -0.4 ;
    TIMHI          =  0.2;
end
tind = find(freq.time <= TIMHI & freq.time >= TIMLO );
XLIM = [TIMLO TIMHI];
YLIM = [FREQLO FREQHI];

nsub=9; % stats done per pupil group

cfg = [];
cfg.latency          = XLIM;
cfg.frequency        = YLIM;
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
%     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000;

cfg.neighbours       = []; %in case no channel data present

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

freq2stats = freq;
freq2stats.dimord = 'subj_freq_time';
freq2stats.label = {'custompooling'};

poolstat=[];
for isoi = 1:2 % occ and motor
    for idiff=3
        for idrug= 1:4
            for ipup = 1:2
                
                freq2stats.powspctrm = squeeze(respavg(pupil_split_lohi{ipup},isoi,:,:, idrug,3,idiff,3,3,3)); %DIMS:
                freq2statszero = freq2stats; %create zero freq to test against
                freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
                poolstat{ipup, isoi, idrug, idiff} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
            end
        end
%         freq2stats.powspctrm = squeeze(respavgsoi(1,isoi,:,:, 4,3,idiff,3,3)); %DIMS:
%         freq2stats2 = freq2stats; %create zero freq to test against
%         freq2stats2.powspctrm = zeros(size(freq2stats.powspctrm));
%         poolstat{ipup, isoi, idrug, idiff} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
        

    end
end
% do independent t stat 

cfg.statistic        = 'indepsamplesT';
% cfg = rmfield(cfg, 'uvar');
freq2stats2 = freq2stats;
isoi = 1;
idiff= 3;

for idrug= 1:4
    freq2stats.powspctrm = squeeze(respavg(pupil_split_lohi{1},isoi,:,:, idrug,3,idiff,3,3)); %DIMS:
    freq2stats2.powspctrm = squeeze(respavg(pupil_split_lohi{2},isoi,:,:, idrug,3,idiff,3,3)); %DIMS:
    poolstat{3, isoi, idrug, idiff} = ft_freqstatistics(cfg, freq2stats, freq2stats2);
end

%% plotting

close all
SAV=1;
set(0,'DefaultAxesFontsize',8)

showstats = 1;
imotor =3;
idiff=3;
istim=3;
iresp=3;
frind = find(faxis <= YLIM(2) & faxis >= YLIM(1));

TYP = '2afc';

for isoi = 1%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
%         for idiff=1:4

    figure;    iplot=0; hold on
    set(gcf, 'Position', [0 -200 375*3 210*4])
    load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
    for ipup = 1:3 % hi, lo
        for idrug = 1:4 % 1:2
            ZLIM = [-0.1 0.1];
            if idrug == 4;          
                ZLIM = [-0.05 0.05];      
%                 showstats = 1;
            end
            iplot = iplot+1; subplot(3,4,iplot); hold on
            
            if ipup<3
                dum = squeeze(mean(respavg(pupil_split_lohi{ipup}, isoi,frind,tind,idrug, imotor, idiff, istim, iresp,3))); %average over subj
            else
                dum = squeeze(mean(respavg(pupil_split_lohi{2}, isoi,frind,tind,idrug, imotor, idiff, istim, iresp,3))) ...
                    - squeeze(mean(respavg(pupil_split_lohi{1}, isoi,frind,tind,idrug, imotor, idiff, istim, iresp,3))); %average over subj
            end         
            
            dum = showstatsTFR(dum, poolstat{ipup, isoi, idrug, idiff}.prob, frind, ZLIM, showstats);
            
            imagesc(taxis(tind), faxis(frind), dum,ZLIM);
            
%             if ipup == 3, ipup = 4; end
            title({sprintf('%s %s %slocked [%g]', motor_conds{imotor}, diff_conds{idiff}, trigger, ZLIM(2)) ;
                sprintf('%s %s %s', sens.leg{isoi}, pharm_conds{idrug}, pupilconds{ipup} )}); % 'FontWeight','bold'
            hold on
            yaxis = [CUTLO,CUTHI];
            plot([0,0],yaxis,'k',[0,0],yaxis,'k');
            
            xlim(XLIM)
            ylim([CUTLO 150])
            YTICKS = [0:50:200];
            set(gca,'Box','off','XTick',-1:0.2:1,...    [-0.5,0,0.5,1,1.5,2,2.5]
                'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                'TickDir','out', 'FontSize', 8);
            LAB=1;
            if LAB %&& (iplot == 9 || iplot == 19)
                xlabel('Time (s)');
                ylabel('Frequency (Hz)');
                zlabel('Response (%)');
                %                                         h=colorbar;
            end
        end %ises
    end
    
    if SAV
        outpath = fullfile(PREOUT, 'poolings');
        warning off; mkdir(outpath); warning on
        outfile = sprintf( '%s%sTFRpooling_%slocked_%s_%s_%sfreq_%s', outpath, filesep, trigger, TYP, sens.leg{isoi}, analysistype, motor_conds{imotor}  );
        disp(outfile)
        export_fig(outfile, '-pdf', '-transparent', '-depsc') %'-png',  '-pdf',
    end;
end

%% plot spectra -0.2 to 0.2 pupilhi and low, 
switch trigger
    case 'stim'
        TIMLO = 0.2;
        TIMHI = 0.4;
    case 'resp'
        TIMLO = -0.2;
        TIMHI = 0.2;
end
tindoi = find((freq.time>= TIMLO) & (freq.time<=TIMHI));

isoi=1;
idrug = 4;
imotor = 3;
idiff = 3;
spectdat = squeeze(mean(respavg(:,isoi,:,tindoi, idrug, imotor, idiff, 3, 3), 4));
% Compute stats for spectra PERMUTATION
stat=[];
nrand = 1000;
nsub=9;
cfg = [];
cfg.frequency        = 'all';
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
%     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = nrand;
cfg.neighbours       = []; %in case no channel data present

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

freq2stats = freq;
freq2stats.dimord = 'subj_freq';
freq2stats.label = {'custompooling'};
freq2stats= rmfield(freq2stats, 'time');

for ipup=1:2 %
    freq2stats.powspctrm = squeeze(spectdat(pupil_split_lohi{ipup},:)); %DIMS: subj sess poi nfreq toi type event ltr
    freq2statszero = freq2stats; %create zero freq to test against
    freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
    stat{ipup} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
    %                             stat{1,ises,ipoi,itoi,itype,ievent,iltr} = 1;
end

% test against each other:
cfg.statistic        = 'indepsamplesT';
freq2stats2 = freq2stats;
cfg = rmfield(cfg, 'uvar');

freq2stats.powspctrm = squeeze(spectdat(pupil_split_lohi{1},:)); %DIMS: subj sess poi nfreq toi type event ltr
freq2stats2.powspctrm = squeeze(spectdat(pupil_split_lohi{2},:)); %DIMS: subj sess poi nfreq toi type event ltr
stat{3} = ft_freqstatistics(cfg, freq2stats, freq2stats2);

% freq2stats1 = freq2stats;
% freq2stats2 = freq2stats;
% freq2stats1.powspctrm = squeeze(respavglat(:,ises,ipoi,:, itoi,itype,1,iltr)); %DIMS: subj sess poi nfreq toi type event ltr
% freq2stats2.powspctrm = squeeze(respavglat(:,ises,ipoi,:, itoi,itype,2,iltr)); %DIMS: subj sess poi nfreq toi type event ltr

                        
%%
%%%%%
dat=[]; sem=[];
dat(1,:) = squeeze(mean(spectdat(pupil_split_lohi{2},:,:), 1));
dat(2,:) = squeeze(mean(spectdat(pupil_split_lohi{1},:,:), 1));
sem(1,:) = squeeze(std(spectdat(pupil_split_lohi{2},:,:), 0, 1)) / sqrt(9);
sem(2,:) = squeeze(std(spectdat(pupil_split_lohi{1},:,:), 0, 1)) / sqrt(9);

showstats=1;
linecolors = {'-r', '-b'};
THR= 0.05;
close all
figure; hold on
for ipup=2:-1:1
    h(ipup) = shadedErrorBar(faxis,dat(ipup,:),sem(ipup,:), linecolors{ipup}, 1)
end
legend([h.mainLine], 'High pupil', ' Low pupil')
legend boxoff
plot([0 faxis(end)], [0 0], '--k' )
set(gca, 'TickDir', 'Out', 'XTick', 0:25:150)
ylim([-0.05 0.2])
xlabel('Frequency (Hz)')
ylabel('Modulation (%)')
title({sprintf('%s %s %slocked', motor_conds{imotor}, diff_conds{idiff}, trigger ) ;
    sprintf('%s %s %g - %g s', sens.leg{isoi}, pharm_conds{idrug}, TIMLO, TIMHI )}); % 'FontWeight','bold'


if showstats
    for ipup=1:3
        pval = stat{ipup}.prob';  % stat from ft_freqstatistics
        
        sigint={};
        BARPOS = [-0.03 -0.02 -0.01];
        BARWID = 8;
        SIGCOL = [0 0 1; 1 0 0; 0 0 0];
        i = 1; sigint{i} = [];
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
                plot(faxis(begsmp:endsmp),...   %begsmp-1:endsmp+1
                    BARPOS(ipup)*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
                    'Color',SIGCOL(ipup,:), 'LineWidth', BARWID)
            end
        end
    end
end