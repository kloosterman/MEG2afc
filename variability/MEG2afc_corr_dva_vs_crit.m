% MEG2afc_corr_dva_vs_crit
cd('/Users/kloosterman/gridmaster2012/MATLAB/MEG_HH_analysis/plotting')

% MEG2afc_setup_paths
trigger = 'stim'

MEG2afc_load_dvadata

%% correlate beta effect, frontal gamma and posterior gamma lateralization effects with rep prob:
% TODO compute lateralization wrt choice (for gamma posterior)

load Repetition_probabilities.mat  % 1=NKID, 2=atomoxetine, 3=placebo, 4=atomoxetine - placebo.
Repetitionprobabilities = Repetitionprobabilities([1, 3:end],:);
% Repetitionprobabilities = Repetitionprobabilities([1, 3:6, 8:end],:);

if length(SUBJ) == 8
    Repetitionprobabilities = Repetitionprobabilities(subj_drugfirst,:);
elseif length(SUBJ) == 10
    Repetitionprobabilities = Repetitionprobabilities(subj_placfirst,:);
else
    warning('Assuming using all subjects!')
%     error('drug of plac first goed zetten!')
end
    

figure
barweb( mean(Repetitionprobabilities(:,2:3)),  std(Repetitionprobabilities(:,2:3)) /sqrt(length(Repetitionprobabilities)),  0.5 )
axis square
ylim([0.4 0.6]); ylabel('Tendency to repeat previous choice')
legend({'atomoxetine' 'placebo' }); legend boxoff
title(sprintf('p = %g, N = %d', randtest1d(Repetitionprobabilities(:,4), zeros(size(Repetitionprobabilities(:,4))), 0, 1000), length(SUBJ) ))
outfile = fullfile(PREOUT, sprintf('Choice_history_behavior' ));
disp(outfile)
export_fig(outfile, '-pdf')

%

cfgsel=[];
cfgsel.avgoverfreq = 'no';
cfgsel.avgovertime = 'yes';
cfgsel.avgoverchan = 'yes';

% SELECT CONDITIONS
imotor = 3; idiff = 2;  istim = 3; iresp = 3;  % correlation should be strongest for hard

% % % frontal gamma
% idrug = 4;   imotor = 3; idiff = 3;  istim = 3; iresp = 3;
% % cfgsel.frequency        =	  [45 55];  %
% cfgsel.latency          =    [-0.4 0];
% cfgsel.channel = motorind;

% % gamma latr
% idrug = 4;   imotor = 3; idiff = 3;  istim = 3; iresp = 4;
% cfgsel.frequency        =	  [40 60];  % check
% cfgsel.latency          =    [-0.4 0.1];

freqsel = [];
freqsel.dimord = 'subj_chan_freq_time';
freqsel.label = chlabel;
freqsel.freq = faxis;
freqsel.time = taxis;

freq2stat=[];
freq2stat.freq = faxis;
freq2stat.dimord = 'chan_subj_freq';
freq2stat.label = {'custompooling'};
freq2stat.time = freqsel.time;

nsub=length(SUBJ);

% loadstat = 0;
% outfile = ['poolstat_' trigger];
% if loadstat;    load(fullfile(PREIN, 'respavg', outfile));    return
% else corrstat={}; end
corrstat={};
powstat={};

cfg = [];
% cfg.latency          = XLIM;
cfg.frequency        = [faxis(1) faxis(end)];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_intersubcorr'; % set below
cfg.correctm         = 'cluster';
% % % %     cfg.correctm         = 'no';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1;
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

for idrug = [2, 4];
    if idrug == 2; RPcol = 3; elseif idrug == 4; RPcol = 4; end  % 1=NKID, 2=atomoxetine, 3=placebo, 4=atomoxetine - placebo.
    
    for isoi=1:3
        cfgsel.channel = sens.ind{isoi};

        if isoi==1; cfgsel.latency = [-0.2 0.2]; elseif isoi==2; cfgsel.latency = [-0.3 0.0]; elseif isoi==3; cfgsel.latency = [-0.2 0.2]; end
        
        if isoi==3; iresp = 4;        else iresp = 3;               end        
        freqsel.powspctrm = squeeze(respavg(:,:,:,:, idrug, imotor, idiff, istim, iresp));
        
        freqselout = ft_selectdata(cfgsel,freqsel); 
        freqselout.dimord = 'subj_chan_freq';
        freq2stat = freqselout;
        
        freqbehav = freq2stat; %create behav freq to correlate with
        for isub=1:nsub
            freqbehav.powspctrm(isub,1,:) = squeeze(Repetitionprobabilities(isub,RPcol));
        end
        cfg.statistic        = 'ft_statfun_intersubcorr'; % ft_statfun_intersubcorr
        corrstat{idrug,isoi} = ft_freqstatistics(cfg, freq2stat, freqbehav);
        % do stats for power
        cfg.statistic        = 'depsamplesT';
        freqzero = freq2stat; %create zero freq to test against
        freqzero.powspctrm = zeros(size(freq2stat.powspctrm));
        powstat{idrug,isoi} = ft_freqstatistics(cfg, freq2stat, freqzero);
        powstat{idrug,isoi}.powspctrm = squeeze(freq2stat.powspctrm);
        % put data in scatterplot var
        cfgsel2 = [];
        cfgsel2.frequency = [40 80];
        cfgsel2.avgoverfreq = 'yes';
        freqscatter = ft_selectdata(cfgsel2, freq2stat)
        scatterdat{idrug,isoi} = [freqscatter.powspctrm  Repetitionprobabilities(:,RPcol) ];
    end
end

%% plotting
close all
figure; iplot=0;
set(gcf, 'Position', [0 -200 375*3 210*4])

for isoi=1:3
    
    for idrug = [2,4]
        iplot=iplot+1; subplot(3,2,iplot); grid on
        
        %         plot(faxis, corrstat{idrug, isoi}.rho);
        [AX,H1,H2] = plotyy(faxis, mean(powstat{idrug, isoi}.powspctrm), faxis, corrstat{idrug, isoi}.rho);
        hold(AX(1)); hold(AX(2))
        plot(AX(1), [0 faxis(end)], [0 0], 'Color', get(H1, 'Color'), 'Linestyle', '--' )
        plot(AX(2), [0 faxis(end)], [0 0], 'Color', get(H2, 'Color'), 'Linestyle', '--' )
        ylim(AX(2), [-1 1])
        set(AX, 'TickDir', 'out', 'XTick', 0:10:150 )
        set(AX(2), 'YTick', -1:0.5:1)
        xlabel('Frequency (Hz)'); ylabel(AX(1), 'Power modulation (%)' ); ylabel(AX(2), 'Choice history correlation (R)' );
        
        % [x,y] = find(corrstat.prob < 0.05);
        % scatter(faxis, corrstat.prob < 0.05)
        % spy(corrstat.prob < 0.05)
        
        THR=0.05; BARWID = 5;  
        for im = 1:2
            if im==1;        pval = powstat{idrug, isoi}.prob';   SIGCOL = get(H1, 'Color'); BARPOS = -0.9;
            else   pval = corrstat{idrug, isoi}.prob';   SIGCOL = get(H2, 'Color');  BARPOS = -0.8;
            end
            
            sigint={};
            i = 1; sigint{i} = [];
            for ismp = 1:length(pval)
                if pval(ismp) < THR
                    sigint{i} = [sigint{i} ismp];
                end % concatenate significant samples in sigint{iint}
                if (ismp < size(pval,1)) && (pval(ismp+1) >= THR) && ~isempty(sigint{i}) % jump to next interval if next sample is not sig
                    i = i + 1;
                    sigint{i} = [];
                end
            end
            % replot significant intervals in different color
            for iint = 1:length(sigint)
                if ~isempty(sigint{iint})
                    begsmp = sigint{iint}(1); %*2
                    endsmp = sigint{iint}(end); %*2
                    plot(AX(2), faxis(begsmp:endsmp),...   %begsmp-1:endsmp+1
                        BARPOS*(ones(1,length(begsmp:endsmp))),...  %begsmp-1:endsmp+1
                        'Color',SIGCOL,'LineWidth',BARWID)
                end
            end
            %         if isoi==1; cfgsel.latency = [-0.2 0.2]; elseif isoi==2; cfgsel.latency = [-0.3 0.0]; end
        end
        if isoi==1; cfgsel.latency = [-0.2 0.2]; elseif isoi==2; cfgsel.latency = [-0.3 0.0]; elseif isoi==3; cfgsel.latency = [-0.2 0.2]; end

        title(sprintf('%s %s %g - %g s, N = %d',  sens.leg{isoi}, pharm_conds{idrug}, cfgsel.latency, nsub ))
    end
end

outfile = fullfile(PREOUT, sprintf('Corr_history_vs_meg' ));
disp(outfile)
export_fig(outfile, '-pdf')


%% plot scatter plots

figure; iplot=0;
set(gcf, 'Position', [0 -200 375*3 210*4])

for isoi=1:3
    for idrug = [2,4]
        iplot=iplot+1; subplot(3,2,iplot); 
        scatter(scatterdat{idrug, isoi}(:,1), scatterdat{idrug, isoi}(:,2));
        axis square; box on
        [r,p] =         corr(scatterdat{idrug, isoi}(:,1), scatterdat{idrug, isoi}(:,2));
        title(sprintf('r = %g, p = %g', r, p))
        if pval<0.05; lsline; end
    end
end
        
