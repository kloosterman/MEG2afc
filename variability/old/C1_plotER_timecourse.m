function [  ] = C1_plotER_timecourse()
%  plots evoked data. The split files have already been planar transformed
%  and combined, so can do stats + plotting directly on the planar
%  gradiometers.

close all; clc; clear;
dbstop if error;

addpath(genpath('~/code/MEG/NewPipeline'));
addpath(genpath('~/code/Tools'));
addpath('~/Documents/fieldtrip/');
ft_defaults; warning off;

% get all those subjects who have a folder
cd ~/Data/MEG-PL;
f = dir('P*');
f = {f(:).name};
for i = 1:length(f), subjects(i) = str2num(f{i}(2:3)); end

subjects = {'GA_pooled', 'GA_placebo', 'GA_atomoxetine', 'GA_donepezil'}; % grand average over drug groups
%subjects = {'GA_pooled'};

% define if we want to plot the mean evoked fields or variance over trials
avgorvar = 'avg';

subjects = 18;

% ==================================================================
% make layout
% ==================================================================

cfg = [];
cfg.layout              = 'CTF275';
lay                     = ft_prepare_layout(cfg);

% colormap
load('~/code/Tools/plotting/colormapsmetgrijs/colormapsmooth4.mat'); % from Niels

for sj = unique(subjects),
    
    close all;
    clearvars -except sj subjects avgorvar lay cmap;
    try
        disp(['Analysing subject ' num2str(sj)]);
    catch
        disp(['Analysing ' sj{1}]);
    end
    
    if isnumeric(sj),
        subjectdata = subjectspecifics(sj);
        cd(subjectdata.lockdir);
    else
        cd('~/Data/MEG-PL/GrandAverage');
    end
    
    % ==================================================================
    % loop over preselected sensor groups
    % ==================================================================
    
    % load file that was created in B3 (on grand average)
    load('~/Data/MEG-PL/sensorselection.mat');
    
    % add the pupil!
    chans(length(chans)+1).group = 'pupil';
    [data, ~, ~, ~] = getData(sj, 'mean', 'avg'); % make sure to select the right chan
    chans(length(chans)).sens = find(~cellfun(@isempty, strfind(data.label, 'EyePupil')));
    
    for thesesens = 1:length(chans),

        clear data;
        
        % baseline correction has already been done in A5
        disp('loading files');
        fields = {'mean', 'stimstrong', 'stimweak', 'respstrong', 'respweak', ...
            'stimweak_respweak', 'stimweak_respstrong', 'stimstrong_respstrong', 'stimstrong_respweak', ...
            'correct', 'incorrect'};
        
        % ==================================================================
        % loop over conditions of interest
        % ==================================================================
        
        for f = 1:length(fields),
            
            % load in data
            [refdata, stimdata, respdata, fbdata] = getData(sj, fields{f}, avgorvar);
            
            % check for the pupil
            if thesesens == length(chans),
                assert(strcmp(refdata.label(chans(thesesens).sens), 'EyePupil'), 'did not find the pupil chan');
            end
            
            if isnumeric(sj),
                
                % average the power spectra over chan, and concatenate them for plotting
                data.(fields{f}).avg = cat(2, ...
                    squeeze(nanmean(refdata.avg(chans(thesesens).sens, :), 1)), ...
                    squeeze(nanmean(stimdata.avg(chans(thesesens).sens, :), 1)), ...
                    squeeze(nanmean(respdata.avg(chans(thesesens).sens, :), 1)), ...
                    squeeze(nanmean(fbdata.avg(chans(thesesens).sens, :), 1)));
                data.(fields{f}).var = cat(2, ...
                    squeeze(nanmean(refdata.var(chans(thesesens).sens, :), 1)), ...
                    squeeze(nanmean(stimdata.var(chans(thesesens).sens, :), 1)), ...
                    squeeze(nanmean(respdata.var(chans(thesesens).sens, :), 1)), ...
                    squeeze(nanmean(fbdata.var(chans(thesesens).sens, :), 1)));
                
                data.(fields{f}).timename   = [refdata.time stimdata.time respdata.time fbdata.time];
                % fool fieldtrip into thinking that the time axis increases
                data.(fields{f}).time       = 1:length(data.(fields{f}).timename);
                data.(fields{f}).dimord     = 'chan_time';
                data.(fields{f}).label      = {'selectedsens'};
                subj = 1;
                
            else % grand average
                
                % average over channels, concatenate over time
                data.(fields{f}).avg = cat(3, ...
                    (nanmean(refdata.individual(:, chans(thesesens).sens, :), 2)), ...
                    (nanmean(stimdata.individual(:, chans(thesesens).sens, :), 2)), ...
                    (nanmean(respdata.individual(:, chans(thesesens).sens, :), 2)), ...
                    (nanmean(fbdata.individual(:, chans(thesesens).sens, :), 2)));
                
                % also save some other important info
                data.(fields{f}).name       = fields{f};
                data.(fields{f}).timename   = [refdata.time stimdata.time respdata.time fbdata.time];
                % fool fieldtrip into thinking that the time axis increases
                data.(fields{f}).time       = 1:length(data.(fields{f}).timename);
                data.(fields{f}).dimord     = 'subj_chan_time';
                data.(fields{f}).label      = {'selectedsens'};
                
                switch fields{f}
                    % for the mean activity, do stats against 0
                    case 'mean'
                        % get the baseline to permute against
                        nulldata                    = data.(fields{f});
                        blsmp                       = find(refdata.time > -.3 & refdata.time < -.1);
                        bl                          = (mean(nulldata.avg(:,1,blsmp), 3));
                        nulldata.avg                = repmat(bl, [1 1 length(nulldata.time)]);
                        
                        subj = size(nulldata.avg, 1);
                        cfgstats                    = setupStats(subj);
                        
                        % compute the statistical against the null distribution
                        data.(fields{f}).stat   = ft_timelockstatistics(cfgstats, nulldata, data.(fields{f}));
                        
                        data.(fields{f}).var    = squeeze(std(data.(fields{f}).avg)) / sqrt(subj);
                        % average the avg across sj for plotting
                        data.(fields{f}).avg    = squeeze(mean(data.(fields{f}).avg));
                        
                    otherwise
                        % stats will be computed later
                end
                
            end % GA
        end % loaded in all subfields
        
        % ==================================================================
        % DEFINE CONTRASTS OF INTEREST
        % ==================================================================
        
        switch chans(thesesens).group
            case 'occipital'
                
                fields = {'mean', 'stim_strong_minus_stim_weak', ...
                    'stim_strong_minus_stim_weak_correct', ...
                    'stim_strong_minus_stim_weak_incorrect', ...
                    'choice_probability'};
                
                % create stronger - weaker stimulus difference
                data = createContrast_doStats(data, ...
                    'stim_strong_minus_stim_weak', 'stimstrong', 'stimweak', subj);
                data = createContrast_doStats(data, ...
                    'stim_strong_minus_stim_weak_correct', 'stimstrong_respstrong', 'stimweak_respweak', subj);
                data = createContrast_doStats(data, ...
                    'stim_strong_minus_stim_weak_incorrect', 'stimstrong_respweak', 'stimweak_respstrong', subj);
                data = createContrast_doStats(data, ...
                    'choice_probability', 'x', 'x', subj);
                
            case 'parietal'
                
                fields = {'mean', 'stim_strong_minus_stim_weak', ...
                    'correct_incorrect', ...
                    'choice_probability'};
                data = createContrast_doStats(data, ...
                    'stim_strong_minus_stim_weak', 'stimstrong', 'stimweak', subj);
                data = createContrast_doStats(data, ...
                    'correct_incorrect', 'correct', 'incorrect', subj);
                data = createContrast_doStats(data, ...
                    'choice_probability', 'x', 'x', subj);
                
            case 'central'
                fields = {'mean', 'correct_incorrect', ...
                    'choice_probability'};
                data = createContrast_doStats(data, ...
                    'correct_incorrect', 'correct', 'incorrect', subj);
                data = createContrast_doStats(data, ...
                    'choice_probability', 'x', 'x', subj);
                
            case 'motor'
                fields = {'mean', 'lateralisation'};
                
                if isnumeric(sj),
                    % get the mapping between strong/weak and hand from subjectspecifics
                    switch mod(sj,2)
                        case 0
                            chans_weak_contra      = rightchans;
                            chans_strong_contra    = leftchans;
                        case 1
                            chans_weak_contra      = leftchans;
                            chans_strong_contra    = rightchans;
                    end
                    
                    % get the full data again
                    [refdata, stimdata, respdata, fbdata] = getData(sj, 'respstrong', avgorvar);
                    
                    % average the power spectra over rpt and chan, and concatenate them for plotting
                    ContraVsIpsi_strong = [squeeze(nanmean(refdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanmean(stimdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanmean(respdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanmean(fbdata.avg(chans_strong_contra, :,:)))] ...
                        - [squeeze(nanmean(refdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanmean(stimdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanmean(respdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanmean(fbdata.avg(chans_weak_contra, :,:)))];
                    
                    ContraVsIpsi_strong_std = [squeeze(nanstd(refdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanstd(stimdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanstd(respdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanstd(fbdata.avg(chans_strong_contra, :,:)))] ...
                        - [squeeze(nanstd(refdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanstd(stimdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanstd(respdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanstd(fbdata.avg(chans_weak_contra, :,:)))];
                    
                    % get the full data again
                    [refdata, stimdata, respdata, fbdata] = getData(sj, 'respweak', avgorvar);
                    
                    % average the power spectra over rpt and chan, and concatenate them for plotting
                    ContraVsIpsi_weak = [squeeze(nanmean(refdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanmean(stimdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanmean(respdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanmean(fbdata.avg(chans_weak_contra, :,:)))] ...
                        - [squeeze(nanmean(refdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanmean(stimdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanmean(respdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanmean(fbdata.avg(chans_strong_contra, :,:)))];
                    
                    ContraVsIpsi_weak_std = [squeeze(nanstd(refdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanstd(stimdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanstd(respdata.avg(chans_weak_contra, :,:))), ...
                        squeeze(nanstd(fbdata.avg(chans_weak_contra, :,:)))] ...
                        - [squeeze(nanstd(refdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanstd(stimdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanstd(respdata.avg(chans_strong_contra, :,:))), ...
                        squeeze(nanstd(fbdata.avg(chans_strong_contra, :,:)))];
                    
                    % motor lateralisation TFR, averaged over both choices
                    data.lateralisation.avg = (ContraVsIpsi_strong + ContraVsIpsi_weak) ./2;
                    data.lateralisation.var = (ContraVsIpsi_strong_std + ContraVsIpsi_weak_std) ./2;
                    data.lateralisation.legtxt = {'contra-ipsi'};
                    
                else % compute lateralisation on GA
                    
                    [refdata, stimdata, respdata, fbdata] = getData(sj, 'lefthand', avgorvar);
                    
                    % average the power spectra chans, and concatenate them for plotting
                    contralat_lefthand = cat(2, ...
                        squeeze(nanmean(refdata.individual(:, rightchans,:), 2)), ...
                        squeeze(nanmean(stimdata.individual(:, rightchans, :), 2)), ...
                        squeeze(nanmean(respdata.individual(:, rightchans, :), 2)), ...
                        squeeze(nanmean(fbdata.individual(:, rightchans, :), 2)));
                    contralat_lefthand_var = cat(2, ...
                        squeeze(nanstd(refdata.individual(:, rightchans,:),[], 2)), ...
                        squeeze(nanstd(stimdata.individual(:, rightchans, :),[], 2)), ...
                        squeeze(nanstd(respdata.individual(:, rightchans, :),[], 2)), ...
                        squeeze(nanstd(fbdata.individual(:, rightchans, :),[], 2))) ...
                        ./ sqrt(subj);
                    ipsilat_lefthand = cat(2, ...
                        squeeze(nanmean(refdata.individual(:, leftchans, :), 2)), ...
                        squeeze(nanmean(stimdata.individual(:, leftchans, :), 2)), ...
                        squeeze(nanmean(respdata.individual(:, leftchans, :), 2)), ...
                        squeeze(nanmean(fbdata.individual(:, leftchans, :), 2)));
                    ipsilat_lefthand_var = cat(2, ...
                        squeeze(nanstd(refdata.individual(:, leftchans, :),[], 2)), ...
                        squeeze(nanstd(stimdata.individual(:, leftchans, :), [],2)), ...
                        squeeze(nanstd(respdata.individual(:, leftchans, :), [],2)), ...
                        squeeze(nanstd(fbdata.individual(:, leftchans, :),[], 2))) ...
                        ./ sqrt(subj);
                    
                    [refdata, stimdata, respdata, fbdata] = getData(sj, 'righthand', avgorvar);
                    
                    % average the power spectra over rpt and chan, and concatenate them for plotting
                    contralat_righthand = cat(2, ...
                        squeeze(nanmean(refdata.individual(:, leftchans, :), 2)), ...
                        squeeze(nanmean(stimdata.individual(:, leftchans, :), 2)), ...
                        squeeze(nanmean(respdata.individual(:, leftchans, :), 2)), ...
                        squeeze(nanmean(fbdata.individual(:, leftchans, :), 2)));
                    contralat_righthand_var = cat(2, ...
                        squeeze(nanstd(refdata.individual(:, leftchans, :),[], 2)), ...
                        squeeze(nanstd(stimdata.individual(:, leftchans, :),[], 2)), ...
                        squeeze(nanstd(respdata.individual(:, leftchans, :),[], 2)), ...
                        squeeze(nanstd(fbdata.individual(:, leftchans, :), [],2))) ...
                        ./ sqrt(subj);
                    ipsilat_righthand = cat(2, ...
                        squeeze(nanmean(refdata.individual(:, rightchans, :), 2)), ...
                        squeeze(nanmean(stimdata.individual(:, rightchans, :), 2)), ...
                        squeeze(nanmean(respdata.individual(:, rightchans,:), 2)), ...
                        squeeze(nanmean(fbdata.individual(:, rightchans, :), 2)));
                    ipsilat_righthand_var = cat(2, ...
                        squeeze(nanstd(refdata.individual(:, rightchans, :), [],2)), ...
                        squeeze(nanstd(stimdata.individual(:, rightchans, :),[], 2)), ...
                        squeeze(nanstd(respdata.individual(:, rightchans,:), [],2)), ...
                        squeeze(nanstd(fbdata.individual(:, rightchans, :),[], 2))) ...
                        ./ sqrt(subj);
                    
                    contralat   = (contralat_lefthand + contralat_righthand) ./2;
                    ipsilat     = (ipsilat_lefthand + ipsilat_righthand) ./2;
                    contralat_var = (contralat_lefthand_var + contralat_righthand_var) ./2;
                    ipsilat_var  = (ipsilat_lefthand_var + ipsilat_righthand_var) ./2;
                    
                    % create artificial TFRs and do stats
                    data.contra = nulldata;
                    data.contra.avg = permute(contralat, [1 3 2]);
                    data.contra.var = contralat_var;
                    data.ipsi   = nulldata;
                    data.ipsi.avg   = permute(ipsilat, [1 3 2]);
                    data.ipsi.var = ipsilat_var;
                    
                    data = createContrast_doStats(data, ...
                        'lateralisation', 'contra', 'ipsi', subj);
                    
                end
                
            case 'lateralfrontal',
                
                fields = {'mean', 'correct_incorrect', ...
                    'choice_probability'};
                data = createContrast_doStats(data, ...
                    'correct_incorrect', 'correct', 'incorrect', subj);
                data = createContrast_doStats(data, ...
                    'choice_probability', 'x', 'x', subj);
                
            case 'pupil',
                
                fields = {'mean', 'correct_incorrect', ...
                    'choice_probability'};
                data = createContrast_doStats(data, ...
                    'correct_incorrect', 'correct', 'incorrect', subj);
                data = createContrast_doStats(data, ...
                    'choice_probability', 'x', 'x', subj);
        end
        
        % ==================================================================
        % PLOT HEAD LAYOUT WITH SENSORS
        % ==================================================================
        
        figure;
        set(0, 'DefaultAxesFontSize', 5, 'DefaultAxesFontName', 'Helvetica');
        colormap(linspecer);
        
        cfgtopo                     = [];
        cfgtopo.style               = 'blank'; % only head shape
        cfgtopo.marker              = '.';
        cfgtopo.markersize          = 1;
        cfgtopo.layout              = lay;
        cfgtopo.comment             = 'no';
        cfgtopo.highlight           = 'on';
        cfgtopo.highlightsymbol     = '.';
        cfgtopo.highlightsize       = 6;
        cfgtopo.highlightchannel    = chans(thesesens).sens;
        
        subplot(5,4,1);
        topodata = refdata;
        if ~isfield(topodata, 'avg'),
            topodata.avg = squeeze(mean(topodata.individual));
        end
        topodata.dimord = 'chan_time';
        if thesesens < length(chans), % skip for pupil
            ft_topoplotER(cfgtopo, topodata);
        else
            axis off; % plot a fake eye?
        end
        title(gca, {sprintf('%s', chans(thesesens).group); 'sensors'});
        
        % ==================================================================
        % PLOT TIMECOURSE FOR THESE SENSORS
        % ==================================================================
        
        for f = 1:length(fields),
            
            % make the plot span several subplots for better visualization
            sh = subplot(5,4,[(f-1)*4+2 (f-1)*4+3 (f-1)*4+4]);
            
            % use my own shaded errorbar (boundedline, so there are no extra lines
            % on the outsides of the errorbars), not supported by FT
            thisdat  = data.(fields{f});
            
            if size(thisdat.var,2) ~= 1 && size(thisdat.avg,1)~=1,
                thisdat.var = permute(thisdat.var, [1 3 2]);
            elseif all(isnan(thisdat.var(:))), % pupil doesnt have any across-chan variance
                thisdat.var = zeros(size(thisdat.var));
            end
            ph       = boundedline(data.mean.time, thisdat.avg, thisdat.var);
            
            if thesesens < length(chans),
                switch avgorvar
                    case 'avg'
                        ylabel('fT');
                        %set(gca, 'YLim', [0 1.5*10^-13]);
                    case 'var'
                        ylabel('variance');
                        %set(gca, 'YLim', [0 1.5*10^-25]);
                end
            else
                ylabel('pupil size');
            end
            
            axis tight;
            
            hold on;
            try
                % add a line to indicate significant timewindow
                mask = thisdat.stat.mask;
                
                ylims = get(gca, 'Ylim'); mask = ((ylims(2)*0.05)+ylims(1))*mask; mask(mask==0) = nan;
                plot(1:length(mask), mask, '.r', 'MarkerSize', 3);
            end
            
            % layout
            set(gca, 'TickDir', 'out', 'YDir', 'normal', 'box', 'off');
            % subfunction for white lines
            plotLines(refdata, [-0.2 0 .5], stimdata, [-.2 0 0.5], respdata, [-.2 0], fbdata, [-0.2 0 0.5]);
            
            % legend
            switch fields{f}
                case 'mean'
                    % no legend
                otherwise
                    l = legend(thisdat.legtxt);
                    lpos = get(l, 'position');
                    lpos(1) = lpos(1) + 0.05;
                    set(l, 'PlotBoxAspectRatioMode', 'manual', ...
                        'PlotBoxAspectRatio', [1 1 1], 'box', 'off', ...
                        'FontSize', 4, 'position', lpos, 'interpreter', 'none');
            end
            
            if f == length(fields), xlabel('time (s)'); end
            if f == 1,
                if isnumeric(sj),
                    switch avgorvar
                        case 'avg'
                            title(sprintf('P%02d - AVERAGE', sj), 'fontweight', 'bold');
                        case 'var'
                            title(sprintf('P%02d - VARIANCE', sj), 'fontweight', 'bold');
                    end
                else
                    switch avgorvar
                        case 'avg'
                            title(sprintf('%s, N = %d - AVERAGE', ...
                                sj{:}, length(refdata.cfg.inputfile)), 'interpreter', 'none', 'fontweight', 'bold');
                        case 'var'
                            title(sprintf('%s, N = %d - VARIANCE', ...
                                sj{:}, length(refdata.cfg.inputfile)), 'interpreter', 'none', 'fontweight', 'bold');
                    end
                end
            else
                title(fields{f}, 'interpreter', 'none', 'fontweight', 'bold');
            end
            
        end
        
        % add an extra subplot without axes to make sure all plots have the
        % same size
        sh = subplot(5,4,20);
        set(gca, 'XColor', 'w', 'YColor', 'w');
        
        % wrap up
        set(gcf, 'PaperPositionMode', 'auto');
        if isnumeric(sj),
            print(gcf, '-depsc', '-painters', sprintf('~/Data/MEG-PL/Figures/P%02d_%s_%ssensors.eps', sj, avgorvar, chans(thesesens).group));
        else
            print(gcf, '-depsc', '-painters', sprintf('~/Data/MEG-PL/Figures/%s_%s_%ssensors.eps', sj{:}, avgorvar, chans(thesesens).group));
        end
        
    end % sensor groups
end % sj or GA pharma groups

end % function

% ==================================================================
% layout, plots lines to indicate event onset
% ==================================================================
function [] = plotLines(refdata, reftp, stimdata, stimtp, respdata, resptp, fbdata, fbtp)

xticks = []; xlabels = {};
for t = 1:length(reftp),
    xticks = [xticks dsearchn(refdata.time', reftp(t))];
    if reftp(t) == 0,
        xlabels = [xlabels 'ref'];
    else
        xlabels = [xlabels reftp(t)];
    end
end

for t = 1:length(stimtp),
    xticks = [xticks length(refdata.time) + dsearchn(stimdata.time', stimtp(t))];
    if stimtp(t) == 0,
        xlabels = [xlabels 'stim'];
    else
        xlabels = [xlabels stimtp(t)];
    end
end

for t = 1:length(resptp),
    xticks = [xticks length(refdata.time) + length(stimdata.time) + ...
        dsearchn(respdata.time', resptp(t))];
    if resptp(t) == 0,
        xlabels = [xlabels 'resp'];
    else
        xlabels = [xlabels resptp(t)];
    end
end

for t = 1:length(fbtp),
    xticks = [xticks length(refdata.time) + length(stimdata.time) + ...
        length(respdata.time) + dsearchn(fbdata.time', fbtp(t))];
    if fbtp(t) == 0,
        xlabels = [xlabels 'fb'];
    else
        xlabels = [xlabels fbtp(t)];
    end
end

set(gca, 'XTick', xticks, 'XTickLabel', xlabels);

% add white lines to indicate transitions between intervals
x = length(refdata.time)+.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);
x = length(refdata.time) + length(stimdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);
x = length(refdata.time) + length(stimdata.time) + length(respdata.time) +.5;
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'w', 'LineStyle', '-', 'LineWidth', 2);

% add dotted  black lines to indicate event onsets
x = dsearchn(refdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '--');
x = length(refdata.time) + dsearchn(stimdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '--');
x = length(refdata.time) + length(stimdata.time) + dsearchn(respdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '--');
x = length(refdata.time) + + length(stimdata.time) + ...
    length(respdata.time) + dsearchn(fbdata.time', 0);
l = line([x x], get(gca, 'YLim')); set(l, 'Color', 'k', 'LineStyle', '--');

end

% ==================================================================
% subfunction to load in data
% ==================================================================
function [refdata, stimdata, respdata, fbdata] = getData(sj, fieldname, avgorvar)

if isnumeric(sj),
    fprintf('loading in data for sj %s and fieldname %s \n', sj, fieldname);
    
    subjectdata = subjectspecifics(sj);
    cd(subjectdata.lockdir);
    
    switch avgorvar
        case 'avg'
            
            % load in the right files
            refdata     = load(sprintf('P%02d_ref_%s_lockbl.mat', sj, fieldname));
            refdata     = refdata.timelock;
            stimdata    = load(sprintf('P%02d_stim_%s_lockbl.mat', sj, fieldname));
            stimdata    = stimdata.timelock;
            respdata    = load(sprintf('P%02d_resp_%s_lockbl.mat', sj, fieldname));
            respdata    = respdata.timelock;
            fbdata      = load(sprintf('P%02d_fb_%s_lockbl.mat', sj, fieldname));
            fbdata      = fbdata.timelock;
        case 'var'
            
            % load in the right files
            refdata     = load(sprintf('P%02d_ref_%s_var.mat', sj, fieldname));
            refdata     = refdata.timelock;
            stimdata    = load(sprintf('P%02d_stim_%s_var.mat', sj, fieldname));
            stimdata    = stimdata.timelock;
            respdata    = load(sprintf('P%02d_resp_%s_var.mat', sj, fieldname));
            respdata    = respdata.timelock;
            fbdata      = load(sprintf('P%02d_fb_%s_var.mat', sj, fieldname));
            fbdata      = fbdata.timelock;
    end
    
    
else
    fprintf('loading in data for sj %s and fieldname %s \n', sj{:}, fieldname);
    
    cd('~/Data/MEG-PL/GrandAverage');
    switch avgorvar
        case 'avg'
            % grand average means
            refdata     = load(sprintf('%s_ref_%s_lockbl.mat', sj{:}, fieldname));
            refdata     = refdata.grandavg;
            stimdata    = load(sprintf('%s_stim_%s_lockbl.mat', sj{:}, fieldname));
            stimdata    = stimdata.grandavg;
            respdata    = load(sprintf('%s_resp_%s_lockbl.mat', sj{:}, fieldname));
            respdata    = respdata.grandavg;
            fbdata      = load(sprintf('%s_fb_%s_lockbl.mat', sj{:}, fieldname));
            fbdata      = fbdata.grandavg;
            
        case 'var'
            % grand average means
            refdata     = load(sprintf('%s_ref_%s_var.mat', sj{:}, fieldname));
            refdata     = refdata.grandavg;
            stimdata    = load(sprintf('%s_stim_%s_var.mat', sj{:}, fieldname));
            stimdata    = stimdata.grandavg;
            respdata    = load(sprintf('%s_resp_%s_var.mat', sj{:}, fieldname));
            respdata    = respdata.grandavg;
            fbdata      = load(sprintf('%s_fb_%s_var.mat', sj{:}, fieldname));
            fbdata      = fbdata.grandavg;
    end
    
end
end

% subfunction that creates a contrast between two conditions and runs stats
function data = createContrast_doStats(data, name, field1, field2, subj)

switch name
    case 'choice_probability'
        
        data.(name).var = ...
            (squeeze(std(data.stimstrong_respstrong.avg - data.stimstrong_respweak.avg)) + ...
            squeeze(std(data.stimweak_respstrong.avg - data.stimweak_respweak.avg))) ...
            ./ 2 * sqrt(subj);
        
        data.(name).avg = ...
            (squeeze(mean(data.stimstrong_respstrong.avg - data.stimstrong_respweak.avg)) + ...
            squeeze(mean(data.stimweak_respstrong.avg - data.stimweak_respweak.avg))) ...
            ./ 2;
           
        % make a text for the legend
        data.(name).legtxt      = 'choice_probability';
        
        %         % get params for stats
        %         cfgstats = setupStats(subj);
        %
        %         data.(name).stat      = ...
        %             ft_timelockstatistics(cfgstats, ...
        %             data.(field1), data.(field2));
        
    otherwise
        if subj > 1,
            data.(name).var    = ...
                [squeeze(std(data.(field1).avg))' / sqrt(subj);...
                squeeze(std(data.(field2).avg))' / sqrt(subj)];
            
            data.(name).avg = ...
                [squeeze(mean(data.(field1).avg))' ;...
                squeeze(mean(data.(field2).avg))' ];
            
            % get params for stats
            cfgstats = setupStats(subj);
            
            data.(name).stat      = ...
                ft_timelockstatistics(cfgstats, ...
                data.(field1), data.(field2));
            
        else
            % single subject, dont do stats for now
            data.(name).var    = ...
                [squeeze(mean(data.(field1).var, 1)) ;...
                squeeze(mean(data.(field2).var, 1)) ];
            
            data.(name).avg = ...
                [squeeze(mean(data.(field1).avg, 1)) ;...
                squeeze(mean(data.(field2).avg, 1)) ];
        end
        
        % make a text for the legend
        data.(name).legtxt      = {field1, field2};
end
end

function cfgstats = setupStats(subj)

% ==================================================================
% prepare for group-level stats
% ==================================================================

cfgstats                  = [];
cfgstats.method           = 'montecarlo'; % permutation test
cfgstats.statistic        = 'ft_statfun_depsamplesT';

% do cluster correction
cfgstats.correctm         = 'cluster';
cfgstats.clusteralpha     = 0.05;
cfgstats.clusterstatistic = 'maxsize'; % weighted cluster mass needs cfg.wcm_weight...
%cfgstats.minnbchan        = 1; % average over chans
cfgstats.tail             = 0;
cfgstats.clustertail      = 0;
cfgstats.alpha            = 0.025;
cfgstats.numrandomization = 1000;

% use only our preselected sensors for the time being
cfgstats.channel          = 'selectedsens';

% specifies with which sensors other sensors can form clusters
cfgstats.neighbours       = []; % only cluster over data and time

design = zeros(2,2*subj);
for i = 1:subj,  design(1,i) = i;       end
for i = 1:subj,  design(1,subj+i) = i;  end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfgstats.design   = design;
cfgstats.uvar     = 1;
cfgstats.ivar     = 2;

end