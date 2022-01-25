% load in basespec for drug/plac and contra ipsi regime

close all

if ismac
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman'; % on the cluster
else
    basepath = '/home/mpib/kloosterman'; % on the cluster
end
PREOUT = fullfile(basepath, 'plots'); mkdir(PREOUT)

freqtype = 'low';

trigger = 'stim';
PREIN = fullfile('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq', freqtype, 'stim/respavg');
load(fullfile(PREIN,  sprintf('basespec_%s_ses%d.mat', 'NK1', 1)))

if strcmp(freqtype, 'low')
    load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/low/stim/NK1/drug_contra/NK1_ses3_20141129_drug_contra_run1_stim_totalpow_freq.mat')
else
    load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/full/stim/NK1/drug_contra/NK1_ses3_20141129_drug_contra_run1_stim_totalpow_freq.mat')
end
NKsensorselection

SUBJ  = { 'NK1'  'NK2'  'NK3'   'NK4'   'NK5'   'NK7'   'NK8'   'NK9'  'NK11'   'NK12'   'NK13'   'NK14' ...
    'NK15'     'NK16'     'NK17'     'NK18'     'NK19'     'NK20'     'NK21' }; %  'NK6' out,  'NK10' don't exist, 'NK2' incomplete

basespec_all = nan(length(SUBJ), size(basespec,1), size(basespec,2), 4);

for isub = 1:length(SUBJ)
    for ises = 1:4
        fprintf('%s ses %d . . .\n',  SUBJ{isub}, ises)
        filesave = fullfile(fileparts(PREIN), 'respavg', sprintf('basespec_%s_ses%d.mat', SUBJ{isub}, ises));
        if ~exist(filesave)
            warning(sprintf('%s not found', filesave ))
            continue
        end            
        load(filesave)
        basespec_all(isub,:,:,ises) = basespec;
    end
end

basespec_avg(:,:,:,1) = mean(basespec_all(:,:,:,1:2),4); % drug
basespec_avg(:,:,:,2) = mean(basespec_all(:,:,:,3:4),4); % plac
basespec_avg(:,:,:,4) = basespec_avg(:,:,:,1) - basespec_avg(:,:,:,2);

% basespec_subjavg = squeeze(mean(basespec_avg)); %avg over subj
% 
% basespec_avg = squeeze(mean(basespec_avg(:,motorind,:),2)); % TODO, select motor, occ

%% %plot basespecs for the 4 conds
SAV=1;
LEG = {'Drug', 'Placebo', '', 'Drug-Placebo'};
linecolors= {{'r','markerfacecolor','r'} {'b','markerfacecolor','b'} {'k','markerfacecolor','k'} {'k','markerfacecolor','k'} };
bandoi = [2 35; 40 150]
bandnames = {'lowfreq'  'highfreq'};

close all
figure;
set(gcf, 'Position', [0 -200 375*3 210*4])
iplot=0;
for iband=1%:2
    for isoi=1:3
        dat = squeeze(nanmean(basespec_avg(:,sens.ind{isoi},:,:),2));
        avg = squeeze(nanmean(dat));
        sem = squeeze(nanstd(dat)/sqrt(length(SUBJ)));
        
        %     plot(freq.freq, basespec_avg, 'Linewidth', 2)
        iplot=iplot+1;
        subplot(2,3,iplot); hold on
        legh = [];
        for idrug = [1:2]
            freqind = find(freq.freq >= bandoi(iband,1) & freq.freq <= bandoi(iband,2));
            h = shadedErrorBar(freq.freq(freqind), avg(freqind,idrug), sem(freqind,idrug), linecolors{idrug}, 1);
            legh(idrug) = h.mainLine;
            set(gca, 'Tickdir', 'out')
        end
        
        xlim(bandoi(iband,:))
        xlim([bandoi(iband,1)-5  bandoi(iband,2)])
        if idrug == 4; plot(get(gca, 'XLim'), [0 0], '--k'); end;
        box on
        
        legend(legh, LEG); legend BOXOFF
        title(sprintf('%s %s', sens.leg{isoi}, bandnames{iband}))
        xlabel('Frequency (Hz)')
        ylabel('Power (T)')
    end
    
    % TODO stats on diff
    
end

if SAV
    outfile = sprintf( 'baselinespectra_%s_%sfreq', LEG{idrug}, freqtype );
    out = fullfile(PREOUT, 'basespec', outfile);
    warning off; mkdir(fileparts(out)); warning on
    disp(out)
    export_fig(out, '-pdf') %'-png',  '-pdf',
end

%% Correlate baseline diff with bias diff
% See MEG2afc_corr_meg_vs_bias

close all
biasdat =  csvread('drug_bias.csv', 1,1); % 19 subj, col 1 = drug
biasdat = fliplr(biasdat); % ALARM! now col 1 = plac
 %split based on plac data.
split = (biasdat(:,1) > 0.5) + 1; % split wrt bias direction: 12 vs 7 subj
% split = (biasdat(:,1) > median(biasdat(:,1))) + 1; %10 have a left bias, 9 a right; JW split met median! % split = right_button_plac < np.median(right_button_plac)

%to correlate drug - plac MEG vs bias
% calculate abs change 
bias_change_abs = abs( biasdat(:,2) - biasdat(:,1));

%
SAV=1;
LEG = {'Drug', 'Placebo', '', 'Drug-Placebo'};
linecolors= {{'r','markerfacecolor','r'} {'b','markerfacecolor','b'} {'k','markerfacecolor','k'} {'k','markerfacecolor','k'} };
bandoi = [5 35; 40 150]
bandnames = {'lowfreq'  'highfreq'};

close all
figure;
set(gcf, 'Position', [0 -200 375*3 210*4])
iplot=0;
for iband=1:2%:2
    freqind = find(freq.freq >= bandoi(iband,1) & freq.freq <= bandoi(iband,2));
    for isoi=1:3
        dat = squeeze(nanmean(basespec_avg(:,sens.ind{isoi},:,:),2));
        avg = squeeze(nanmean(dat));
        sem = squeeze(nanstd(dat)/sqrt(length(SUBJ)));
        
        iplot=iplot+1;        subplot(2,3,iplot); hold on
        for idrug = 4 %[1:2]

            h = shadedErrorBar(freq.freq(freqind), avg(freqind,idrug), sem(freqind,idrug), linecolors{idrug}, 1);
            
            rcoeff = [];
            for ifreq = freqind
                dum = squeeze(dat(:,ifreq,idrug));               
                rcoeff(ifreq) = corr( bias_change_abs(~isnan(dum)), dum(~isnan(dum)), 'type', 'Pearson' );
            end
            [AX,H1,H2] = plotyy(freq.freq(freqind), nan(size(freqind)), freq.freq(freqind), rcoeff(freqind));
            hold(AX(2))
            plot(AX(2), bandoi(iband,:), [0 0], 'Color', 'r', 'Linestyle', '--');

%             [AX,H1,H2] = plotyy(freq.freq(freqind), avg(freqind,idrug), freq.freq(freqind), rcoeff(freqind));
            set(AX,{'ycolor'},{'k';'r'})  % Left color red, right color blue...
            set(AX(2), 'YLIM', [-1 1])
            set(AX(2), 'YTick', -1:0.5:1)
%             legh(idrug) = h.mainLine;



        end
        
%         xlim(bandoi(iband,:))
        xlim(AX(1),[bandoi(iband,1)-5  bandoi(iband,2)])
        xlim(AX(2),[bandoi(iband,1)-5  bandoi(iband,2)])
        box off
        
%         legend(legh, LEG); legend BOXOFF
        title(sprintf('%s %s', sens.leg{isoi}, bandnames{iband}))
        xlabel('Frequency (Hz)')
        ylabel(AX(1), 'Power (T)')
        ylabel(AX(2), 'Correlation with abs. bias')
    end
    
    
end

if SAV
    outfile = sprintf( 'baselinespectra_absbiascorr_drug-plac_%s', LEG{idrug} );
    out = fullfile(PREOUT, 'basespec', outfile);
    warning off; mkdir(fileparts(out)); warning on
    disp(out)
    export_fig(out, '-pdf') %'-png',  '-pdf',
end


