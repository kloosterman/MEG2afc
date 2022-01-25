%% plot TFR's for selected poolings
% function plotTFR_2AFC(~)
% plots SOI-TFRs
% average of selected sensors for multiple subjects
% sorted according to stimulus condition
% MEG hh conditions:

% drug placebo
% easy hard,
% pupil hi lo

% H, FA, M, CR: Focus on STIM LEFT

loadrespavg = 1; %load in subj respavg
trigger = 'stim';

ALLSUBJ  = { 'NK1'  'NK3'   'NK4'   'NK5'   'NK7'   'NK8'   'NK9'  'NK11'   'NK12'   'NK13'   'NK14' ...
    'NK15'     'NK16'     'NK17'     'NK18'     'NK19'     'NK20'     'NK21' }; %  'NK6' out,  'NK10' don't exist, 'NK2' incomplete

PREIN = fullfile('/mnt/homes/home022/nkloost1/projectdata/2afc/preproc/variance');

PREOUT = '/mnt/homes/home022/nkloost1/plots/';
mkdir(PREOUT)

% addpath(genpath('/Users/niels/Dropbox/PROJECTS/Hamburg_MEG/MEG_Hamburg_shared/scripts/MEG_HH_analysis/'))
addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
ft_defaults

% stim left right
% resp left or right
% NOTE: comb dim levels not put in b/c space lims
pharm_conds = {'atomox' 'placebo' 'drugscomb' 'atomox - placebo'};
motor_conds = {'pressipsi' 'presscontra' 'motorcomb' 'pressipsi - presscontra'};
diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};
pupilconds = {'pupillow' 'pupilhi' 'pupilcomb' 'pupilhi-lo'};
stim_conds = {'stimleft' 'stimright' 'allstim' 'stimleft-right'};
choice_conds = {'choiceleft' 'choiceright' 'allchoice' 'choiceleft-right'};

soi_optimal_occ % get the sensors of interest (soi) list
load('sensorselection.mat')
SOINsel = { 'occ-l'    'occ-r'  'occ-spec' 'frontal'}; %  'motor-l' 'motor-r' %         SOIN = {'occ-l', 'occ-r', 'occ-spec', 'frontal'};

resamplefs = 100;
BASELO = -0.25;
BASEHI = 0;
TIMLO = -1.5; TIMHI = 0.3; % resp

taxis = TIMLO:1/resamplefs:TIMHI;
basetind = find(taxis >= BASELO & taxis <= BASEHI)

% rt = cell(length(SUBJ),2,2,2,2,2,3); % subj pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3) 
baseline = 'trial';
if loadrespavg % todo load in respavg
    for isub = 1:length(ALLSUBJ)
%         fprintf('Loading subjrespavg %s %s freq trials %s . . . \n',SUBJ{isub}, TYP, AVG)
        fprintf('loading subrespavg %s\n', ALLSUBJ{isub})
        subjrespavgin = sprintf('%s/respavg/%s_%s.mat', PREIN, ALLSUBJ{isub}, baseline);
        load(subjrespavgin);
        if ~exist('respavg', 'var')
            respavg = nan(length(ALLSUBJ),length(chlabel),length(taxis),2,2,3,3,3,1);
        end
        %         respavg(isub,:,:,:, :,:,:,:,:) = subjrespavg;
                respavg(isub,:,:,:, :,:,:,:,:) = subjrespvar;
        
%         clear subjrespavg
    end
end


% flip motor dimension!
% for imotor=2 (contra) swap levels 1 and 2 of motor dim:
% b/c 1 means press left, but subject meant right, and vice versa! 
% I.e. dim 9 is CHOICE, not the button that was pressed
respavg(:,:,:, :, 2, :, :, 1:2, :)  = flipdim(respavg(:,:,:, :, 2, :, :, 1:2, :), 8);
%     respavgtemp = respavg;
%     respavg(:,:,:,:, :, 2, :, :, 2, :) =  respavgtemp(:,:,:,:, :, 2, :, :, 1, :);
%     respavg(:,:,:,:, :, 2, :, :, 1, :) =  respavgtemp(:,:,:,:, :, 2, :, :, 2, :);
%     clear respavgtemp


[nsub, nchan,ntim,npharm,nmotor,ndiff,nstim, nchoice,npup,]=size(respavg)
% respavgdimord = 'nsub_nchan_nfrq_ntim_npharm_nmotor_ndiff_npup_nsdt'
respavgdimord = 'nsub_nchan_ntim_npharm_nmotor_ndiff_nstim_nchoice_npup'
size(respavg)

% filesave = [respavgout filesep 'respavg_' trigger '.mat'];
% fprintf('Saving %s . . . \n', filesave)
% save(filesave, 'respavg',   '-v7.3' )


% soiUnilateral_cmb_CABMSI_optimal_occ % get the sensors of interest (soi) list
% % %
% % examplefreqpath = fullfile(PREIN, SUBJ{1}, sesdirs{1}, trigger{itrg});
% % w = what(examplefreqpath);
% % load(fullfile(examplefreqpath, w.mat{1}))
%
% switch analysistype{iat}
%     case 'low'
% %         load('/home/niels/MIBmeg/freqhigh/s1_TK/230611/stim/TK230611_mib2_400_type2event1_tsss_totalpow_freq.mat')
%         load('TK230611_mib7_400_type2event1_tsss_totalpow_freq')
%         CUTLO  = 3;         % low cutoff for plots
%         CUTHI  = 35;       % high cutoff for plots
%         ZLIM=[-0.10 0.10];
% %         ZLIM=[-0.40 0.40];
%         YTICKS = [10 20 30];
%     case 'high'
%         CUTLO  = 36;         % low cutoff for plots
%         CUTHI  = 140;       % high cutoff for plots
%         freq.freq = 36:2:150;
% %         ZLIM=[-0.025 0.025];
%         ZLIM=[-0.1 0.1];
% %         YTICKS = [40 80 120];
%         YTICKS = [50 75 100 125];
%     case 'full'
%         CUTLO  = 5;         % low cutoff for plots
%         CUTHI  = 140;       % high cutoff for plots
%         bandoi = [50 100];
%         ZLIM=[-0.12 0.12];
% %         ZLIM=[-0.1 0.1; -0.06 0.06; -0.06 0.06; -0.06 0.06 ]; %single subj
% end
% frind = find((freq.freq>=CUTLO) & (freq.freq<=CUTHI));
% faxis = freq.freq(frind);

% % sois = {
% %     'occ-globresp' 'frontal_trbeta-l' 'frontal_trbeta-r'   'frontal_trbeta'...
% %         'frontal_sustbeta-l' 'frontal_sustbeta-r'          'frontal_sustbeta' }; % 1, 4 and 7 most important
% sois = {
%     'occ-globresp' 'frontal_trbeta-l' 'frontal_trbeta-r'   'frontal_trbeta'...
%         'frontal_sustbeta-l' 'frontal_sustbeta-r'          'frontal_sustbeta'...  % 1, 4 and 7 most important
%         'motor' 'motor-l' 'motor-r'}; % 8,9 10
%
%     %     'occ-globresp_l'
%     %     'occ-globresp_r'
%     %     'centr-l'
%     %     'centr-r'
%     %     'motor'
%     %};
%
%     thepooling = sois{8}
%     plotTOPOblank({thepooling},1)

% soi_optimal_occ % get the sensors of interest (soi) list
load('sensorselection.mat')
% SOINsel = { 'occ-l'    'occ-r'  'occ-spec' 'frontal'}; %  'motor-l' 'motor-r' %         SOIN = {'occ-l', 'occ-r', 'occ-spec', 'frontal'};
SOINsel = { 'motor'    'occipital'}; %  'motor-l' 'motor-r' %         SOIN = {'occ-l', 'occ-r', 'occ-spec', 'frontal'};

%% average over subj
% soinameind = find(strcmp(  'occ-spec',  SOIN));
% sensind = match_str(chlabel,ft_channelselection(SOI{soinameind},chlabel));

% respavg
% 1     2    3      4        5          6       7       8       9       10
% subj nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2) pupil(2/3)

% respavgpool
% 1     2     3         4         5        6       7       8          |9       
% nchan nfreq ntimebins pharma(2) motor(2) diff(2) stim(2) choice(2)  |pupil(2/3)



if size(respavg,1) > 1 %average subj
    
    respavgpool = squeeze(mean(respavg)); %average over subj
    respavgpool(:,:,:,3,:,:,:) = mean(respavgpool(:,:,:,1:2,:,:,:), 4); %avg over motor regime
    respavgpool(:,:,3,:,:,:,:) = mean(respavgpool(:,:,1:2,:,:,:,:), 3); %avg over pharma
    respavgpool(:,:,4,:,:,:,:) = respavgpool(:,:,1,:,:,:,:) - respavgpool(:,:,2,:,:,:,:) ; %pharma diff
    respavgpool(:,:,:,4,:,:,:) = respavgpool(:,:,:,1,:,:,:) - respavgpool(:,:,:,2,:,:,:) ; %motor ipsi - contra
    respavgpool(:,:,:,:,4,:,:) = respavgpool(:,:,:,:,1,:,:) - respavgpool(:,:,:,:,2,:,:) ; %easy-hard
    respavgpool(:,:,:,:,:,4,:) = respavgpool(:,:,:,:,:,1,:) - respavgpool(:,:,:,:,:,2,:) ; %stim left-right    
    respavgpool(:,:,:,:,:,:,4) = respavgpool(:,:,:,:,:,:,1) - respavgpool(:,:,:,:,:,:,2) ; %choice left-right
else
    respavgpool = squeeze(respavg); %or just squeeze
end

% respavgpool(:,:,:,:,:,3) = respavgpool(:,:,:,:,:,1) - respavgpool(:,:,:,:,:,2);

% %if strfind(trigger{itrg}, 'stim')
%     % load('/home/niels/MIBmeg/freqlow/respavg/stim/reactiontimes.mat')
% %     load('/Volumes/Time Machine Backups/data/MIBmeg/freqlow/respavg/stim/reactiontimes.mat')
%     load('reactiontimes.mat') 
%     %     meanrt = squeeze(mean(cellfun(@mean,rt))) /1000; % for ievents and itypes
%     medianrt = squeeze(mean(cellfun(@median,rt))) /1000; % for ievents and itypes dims sub sess type event TAKE MEAN MEDIAN over subj!
% %     %     FE median analysis:
% %     rtcat=squeeze(rt(:,3,:,:)); %only take concat over sessions
% %     clear rt
% %     for itype = 1:3 %:ncond
% %         for ievent = 1:2 % off and on
% %             rt{itype,ievent}=vertcat(rtcat{:,itype,ievent});
% %         end
% %     end
% %     medianrt = cellfun(@median,rt) /1000; % median
% % end

%% freqstatistics on TFR
% nsub=11;
% respavg(:,:,:,:,:,:,3) = respavg(:,:,:,:,:,:,2) - respavg(:,:,:,:,:,:,1);
% respavgstat = squeeze(mean(respavg(:,:,sensind,:,:,:,:),3));
%
% cfg = [];
% cfg.frequency        = 'all';
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
% cfg.correctm         = 'cluster';
% %     cfg.correctm         = 'no';
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
% freq2stats = freq;
% freq2stats.dimord = 'subj_freq_time';
% freq2stats.label = {'custompooling'};
%
% stat=[];
% for ises=1:3 %1:3
%     for itype=1:3 %1:3 %1:3 %1:3
%
%         for ievent=3 %1:2 %3 %test d and r against 0
%             freq2stats.powspctrm = squeeze(respavgstat(:,ises,:,:,itype,ievent)); %DIMS: subj sess poi nfreq toi type event ltr
%             freq2statszero = freq2stats; %create zero freq to test against
%             freq2statszero.powspctrm = zeros(size(freq2stats.powspctrm));
%             stat{ises,itype,ievent} = ft_freqstatistics(cfg, freq2stats, freq2statszero);
%         end
%         % test disapp vs reapp
%         %                         cfg2 = [];
%         %                         cfg2.frequency        = 'all';
%         %                         cfg2.method           = 'stats';
%         %                         cfg2.tail             = 0;
%         %                         cfg2.alpha            = 0.05;
%         %                         cfg2.statistic = 'paired-ttest';   %paired-ttest  ttest2
%         %                         cfg2.design = cfg.design(2,:);
%
%         %                         cfg.correctm         = 'no';
%
%         %     freq2stats1 = freq2stats;
%         %     freq2stats2 = freq2stats;
%         %     freq2stats1.powspctrm = squeeze(respavglat(:,ises,ipoi,:, itoi,itype,1,iltr)); %DIMS: subj sess poi nfreq toi type event ltr
%         %     freq2stats2.powspctrm = squeeze(respavglat(:,ises,ipoi,:, itoi,itype,2,iltr)); %DIMS: subj sess poi nfreq toi type event ltr
%         %     stat{1,ises,ipoi,itoi,itype,3,iltr} = ft_freqstatistics(cfg, freq2stats1, freq2stats2);
%     end
% end


%% 
close all
load('sensorselection.mat')
% SOINsel = { 'occ-l'    'occ-r'  'occ-spec' 'frontal'}; %  'motor-l' 'motor-r' %         SOIN = {'occ-l', 'occ-r', 'occ-spec', 'frontal'};
% SOINind = find(strcmp(SOINsel{1}, {chans.group}));
sensind = chans(1).sens; %occ

resamplefs = 100;
TIMLO = -0.5;
TIMHI = 1.5;
taxis = TIMLO:1/resamplefs:TIMHI;

% dat = squeeze(mean(subjrespvar(sensind,:,:,:,3,3,3)));
dat = squeeze(mean(respavgpool(sensind,:,:,:,3,3,3)));
dat(:,:,3) = squeeze(mean(dat,3));
figure; hold on
COL = ['g', 'k']
for ipharm=1:2
    for imotor=1:2 %2%:2
        plot(taxis, dat(:,ipharm, imotor), COL(ipharm))
    end
end
        legend({'drug', 'placebo'})
        



















%% plotting: 2AFC TFR drug vs placebo: motor and visual cortex
close all
SAV=1;
set(0,'DefaultAxesFontsize',8)


% XLIM = [-0.3 1.3; -1.5 1.5 ]; % For fig 2
XLIM = [-1.5 0.25; -1.5 1.5 ]; % For fig 2
% XLIM = [-1.5 1.5]; % For fig S1
ZLIM = [-0.05 0.05]

showstats = 0;

TYP = '2afc'
for isoi = 1% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
%     sensind = match_str(chlabel,ft_channelselection(SOI{soinameind},chlabel));
    SOINind = find(strcmp(SOINsel{isoi}, {chans.group}));
    sensind = chans(SOINind).sens;
    
    for imotor = 1:4 % [1,2,4] % 
        
        figure;    iplot=0; hold on
        load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
        
        for idiff=[1,2,4]
            for ipharm = 1:4 % 1:2
                
                for istim = 3 % off and on
                    for ichoice = 3 %1:2% %:ncond
                        iplot = iplot+1; subplot(3,4,iplot); hold on
                        
                        dum = squeeze(mean(respavgpool(sensind,frind,:,ipharm, imotor, idiff, istim, ichoice))); %average over sensors
                        
                        scale = ZLIM; %{inorm}(1,:);
                        
                        %implement stats
                        if showstats
                            nsalpha=0.2;
                            pval = squeeze(stat{ises,itype,ievent}.prob);
                            pval = pval(frind,:);
                            thrcluster = zeros(size(pval));
                            thrcluster = thrcluster(frind,:);
                            thrcluster = thrcluster+nsalpha; % set non significant alpha
                            thrcluster(find(pval<0.05)) = 1;
                            % Transform cdat-values to have a 0-64 range, dependent on clim
                            % (think of it as the data having an exact range of min=clim(1) to max=(clim2), convert this range to 0-64)
                            %       tmpcdat = (tmpcdat + -clim(1)) * (64 / (-clim(1) + clim(2)));
                            tmpcdat = (dum + -ZLIM(1)) * (256 / (-ZLIM(1) + ZLIM(2)));
                            
                            % ind->rgb->hsv ||change saturation values||  hsv->rgb ->  plot
                            rgbcdat = ind2rgb(uint8(floor(tmpcdat)), colormap);
                            hsvcdat = rgb2hsv(rgbcdat);
                            hsvcdat(:,:,2) = hsvcdat(:,:,2) .* thrcluster;
                            dum = hsv2rgb(hsvcdat);
                        end
                        
                        imagesc(taxis,faxis,dum,scale);
                        
                        title({sprintf('%s %s %slocked [%g]', motor_conds{imotor}, diff_conds{idiff}, trigger, ZLIM(2)) ;
                            sprintf('%s Sess: %s', SOINsel{isoi}, pharm_conds{ipharm} )}); % 'FontWeight','bold'
                        hold on
                        yaxis = [CUTLO,CUTHI];
                        plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                        %                 if strfind(trigger{itrg}, 'stim')
                        %                     plot([medianrt(ises,itype,ievent),medianrt(ises,itype,ievent)],yaxis,'--k',[medianrt(ises,itype,ievent),medianrt(ises,itype,ievent)],yaxis,'--k'); %concat over sessions
                        %                 elseif strfind(trigger{itrg}, 'resp')
                        %                     plot([0-medianrt(ises,itype,ievent),0-medianrt(ises,itype,ievent)],yaxis,'--k',[0-medianrt(ises,itype,ievent),0-medianrt(ises,itype,ievent)],yaxis,'--k'); %concat over sessions
                        %                 end
                        
                        xlim(XLIM(1,:))
                        ylim([CUTLO 150])
                        YTICKS = [0:25:200];
                        set(gca,'Box','off','XTick',[-1 -0.5 0 0.5 1],...    [-0.5,0,0.5,1,1.5,2,2.5]
                            'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                            'TickDir','out', 'FontSize', 8);
                        LAB=1;
                        if LAB %&& (iplot == 9 || iplot == 19)
                            xlabel('Time (s)');
                            ylabel('Frequency (Hz)');
                            zlabel('Response (%)');
%                             h=colorbar;
                        end
                    end
                end
            end %itype
        end % ievent
        if SAV
            %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
            outpath = fullfile(PREOUT, 'poolings');
            
            warning off; mkdir(outpath); warning on
            outfile = sprintf( '%s%sTFRpooling_%slocked_%s_%s_%sfreq_%s', outpath, filesep, trigger, TYP, chans(SOINind).group, analysistype, motor_conds{imotor}  );
            disp(outfile)
%             print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
            print('-dpng', outfile) %, '-cmyk'
        end;
    end %ises
end


%% Plot single sessions averaged across conditions

close all
SAV=1;

set(0,'DefaultAxesFontsize',8)

XLIM = [-0.3 1.3; -1.5 1.5 ]; % For fig 2
% XLIM = [-1.5 1.5]; % For fig S1
ZLIM = [-0.5 0.5]

showstats = 0;

for isoi=2
    SOINind = find(strcmp(SOINsel{isoi}, {chans.group}));
    sensind = chans(SOINind).sens;
    for ipharm = 1:2
        for imotor = 1:2
            figure;    iplot=0; hold on
            load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
            
            for isub = 1:length(SUBJ)
                iplot = iplot+1; subplot(3,4,iplot); hold on

                dum = squeeze(mean(respavg(isub,sensind,frind,:,ipharm, imotor, 3, 3, 3))); %average over sensors
                
                scale = ZLIM; %{inorm}(1,:);
                %implement stats
                if showstats
                    nsalpha=0.2;
                    pval = squeeze(stat{ises,itype,ievent}.prob);
                    pval = pval(frind,:);
                    thrcluster = zeros(size(pval));
                    thrcluster = thrcluster(frind,:);
                    thrcluster = thrcluster+nsalpha; % set non significant alpha
                    thrcluster(find(pval<0.05)) = 1;
                    % Transform cdat-values to have a 0-64 range, dependent on clim
                    % (think of it as the data having an exact range of min=clim(1) to max=(clim2), convert this range to 0-64)
                    %       tmpcdat = (tmpcdat + -clim(1)) * (64 / (-clim(1) + clim(2)));
                    tmpcdat = (dum + -ZLIM(1)) * (256 / (-ZLIM(1) + ZLIM(2)));
                    
                    % ind->rgb->hsv ||change saturation values||  hsv->rgb ->  plot
                    rgbcdat = ind2rgb(uint8(floor(tmpcdat)), colormap);
                    hsvcdat = rgb2hsv(rgbcdat);
                    hsvcdat(:,:,2) = hsvcdat(:,:,2) .* thrcluster;
                    dum = hsv2rgb(hsvcdat);
                end
                
                imagesc(taxis,faxis,dum,scale);
                
                title({sprintf('Subj %s %s %s [%g]', SUBJ{isub}, motor_conds{imotor}, diff_conds{idiff}, ZLIM(2)) ;
                    sprintf('%s Sess: %s %slocked ', SOINsel{isoi}, pharm_conds{ipharm}, trigger )}); % 'FontWeight','bold'
                hold on
                yaxis = [CUTLO,CUTHI];
                plot([0,0],yaxis,'k',[0,0],yaxis,'k');
                %                 if strfind(trigger{itrg}, 'stim')
                %                     plot([medianrt(ises,itype,ievent),medianrt(ises,itype,ievent)],yaxis,'--k',[medianrt(ises,itype,ievent),medianrt(ises,itype,ievent)],yaxis,'--k'); %concat over sessions
                %                 elseif strfind(trigger{itrg}, 'resp')
                %                     plot([0-medianrt(ises,itype,ievent),0-medianrt(ises,itype,ievent)],yaxis,'--k',[0-medianrt(ises,itype,ievent),0-medianrt(ises,itype,ievent)],yaxis,'--k'); %concat over sessions
                %                 end
                
                xlim(XLIM(1,:))
                ylim([CUTLO 150])
                YTICKS = [0:25:200];
                set(gca,'Box','off','XTick',[-1 -0.5 0 0.5 1],...    [-0.5,0,0.5,1,1.5,2,2.5]
                    'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
                    'TickDir','out', 'FontSize', 8);
                LAB=1;
                if LAB %&& (iplot == 9 || iplot == 19)
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    zlabel('Response (%)');
                    %                     h=colorbar;
                end
            end
            if SAV
                %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                outpath = fullfile(PREOUT, 'poolings');
                
                warning off; mkdir(outpath); warning on
                outfile = sprintf( '%s%sTFRpooling_%slocked_%s_%s_%sfreq_%s_%s', outpath, filesep, trigger, TYP, chans(SOINind).group, analysistype, motor_conds{imotor}, pharm_conds{ipharm}  );
                disp(outfile)
%                 print('-dpng', outfile) %, '-cmyk'
                print('-dpdf', '-adobecs', '-painter', outfile) %, '-cmyk'
            end;
        end
    end
end

%% multiplot
close all
set(0,'DefaultAxesFontsize',12)
SAV=1
freq.dimord = 'chan_freq_time';

cfg=[];
load('colormap170613.mat'); %colormap(cmap);  %cmap = get(gcf, 'Colormap')
cfg.colormap = cmap;
cfg.layout = 'neuromag306cmb.lay';       %neuromag306all neuromag306mag neuromag306planar
cfg.zlim = [-0.2 0.2];

for ises = 3 %1:3
    for itype = 1%:2 %1:2% %:ncond
        for ievent = 1:2 % off and on
            figure;    iplot=0; hold on
            freq.powspctrm = squeeze(respavgpool(ises, :,frind,:,itype, ievent));
            ft_multiplotTFR(cfg, freq);
            title({sprintf('%s %slocked ', CON{itype,ievent}, trigger{1}) ; sprintf('Sess: %s', sesnames{ises} )}); % 'FontWeight','bold'
            if SAV
                %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                outpath = fullfile(PREOUT, 'multiplot');
                
                warning off; mkdir(outpath); warning on
                outfile = sprintf( '%s%smultiplot_%slocked_%s_%sfreq_%s', outpath, filesep, trigger{itrg}, CON{itype,ievent}, analysistype{iat}, sesnames{ises}  );
                disp(outfile)
                orient landscape
                %                 orient portrait
                %         print('-depsc', '-adobecs', '-painter', outfile) %, '-cmyk'
                print('-dpdf', '-adobecs', '-painter', outfile) %, '-cmyk'
                %                 print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
                %         print('-dpdf', outfile)
            end;
        end
    end
end

% multiplot, separately for left and right pressing subjects

% subtr=[];
% poi = {'occ-globresp' 'motor'}; %poolings of interest
% % poi = {'occ-globresp' 'occ-globresp'}; %poolings of interest
% % poi = { 'frontal_trbeta' 'frontal_sustbeta'}; %poolings of interest
% motl = find(strcmp('motor-l',  SOIN));
% motr = find(strcmp('motor-r',  SOIN));
% % % % % % % subj nr       1          2          3          4          5          6          7          8          9          10
% arealat = {'LatOcctotargetloc', 'LatMotortoBP', 'LatOcctoBP', 'LatMotortotargetloc'}; % legends
% subtr(:,1,:) = [motr,motl; motr,motl; motr,motl; motr,motl; motl,motr; motl,motr; motr,motl; motr,motl; motl,motr; motr,motl; motl,motr]; %contra - ipsi motor rel. to target loc
% subtr(:,2,:) = [motr,motl; motr,motl; motr,motl; motr,motl; motl,motr; motl,motr; motl,motr; motl,motr; motl,motr; motr,motl; motl,motr]; %contra - ipsi motor rel. to resp hand, dims subj poolingind area.!!!
SAV=1;
handleg = {'left'; 'right'};  %lefthanders, righthanders
LRhanders = { [1:4, 10]; [5:9,11]}; %lefthanders, righthanders

set(0,'DefaultAxesFontsize',12)

freq.dimord = 'chan_freq_time';

cfg=[];
load('colormap170613.mat'); %colormap(cmap);  %cmap = get(gcf, 'Colormap')
cfg.colormap = cmap;
cfg.layout = 'neuromag306cmb.lay';       %neuromag306all neuromag306mag neuromag306planar
cfg.zlim = [-0.2 0.2];

close all
for ihand=1:2 % left and right handers
    respavgpoolhand = squeeze(mean(respavg( LRhanders{ihand},:,:,:,:,:,:))); %average over subj
    for ises = 3 %1:3
        for itype = 1%:2 %1:2% %:ncond
            for ievent = 1:2 % off and on
                figure;    iplot=0; hold on
                
                freq.powspctrm = squeeze(respavgpoolhand(ises, :,frind,:,itype, ievent));
                ft_multiplotTFR(cfg, freq);
                title({sprintf('%s %slocked %s-handers', CON{itype,ievent}, trigger{1}, handleg{ihand}) ; sprintf('Sess: %s', sesnames{ises} )}); % 'FontWeight','bold'
                if SAV
                    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                    outpath = fullfile(PREOUT, 'multiplot');
                    
                    warning off; mkdir(outpath); warning on
                    outfile = sprintf( '%s%smultiplot_%slocked_%s_%sfreq_%s_%s-handers', outpath, filesep, trigger{itrg}, CON{itype,ievent}, analysistype{iat}, sesnames{ises}, handleg{ihand}  );
                    disp(outfile)
                    orient landscape
                    %                     orient portrait
                    %         print('-depsc', '-adobecs', '-painter', outfile) %, '-cmyk'
                    print('-dpdf', '-adobecs', '-painter', outfile) %, '-cmyk'
                    %                     print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
                    %         print('-dpdf', outfile)
                end;
            end
        end
    end
end




% multiplot, split up for subjects with target left and right

% subtr=[];
% poi = {'occ-globresp' 'motor'}; %poolings of interest
% % poi = {'occ-globresp' 'occ-globresp'}; %poolings of interest
% % poi = { 'frontal_trbeta' 'frontal_sustbeta'}; %poolings of interest
% motl = find(strcmp('motor-l',  SOIN));
% motr = find(strcmp('motor-r',  SOIN));
% % subj nr       1          2          3          4          5          6          7          8          9          10      11
% arealat = {'LatOcctotargetloc', 'LatMotortoBP', 'LatOcctoBP', 'LatMotortotargetloc'}; % legends
% subtr(:,1,:) = [motr,motl; motr,motl; motr,motl; motr,motl; motl,motr; motl,motr; motr,motl; motr,motl; motl,motr; motr,motl; motl,motr]; %contra - ipsi motor rel. to target loc
% subtr(:,2,:) = [motr,motl; motr,motl; motr,motl; motr,motl; motl,motr; motl,motr; motl,motr; motl,motr; motl,motr; motr,motl; motl,motr]; %contra - ipsi motor rel. to resp hand, dims subj poolingind area.!!!
SAV=1;
tlocleg = {'left'; 'right'};  %lefthanders, righthanders
LRtlocsubjects = { [1:4,7,8,10]; [5,6,9,11]}; %left side, right side

set(0,'DefaultAxesFontsize',12)

freq.dimord = 'chan_freq_time';

cfg=[];
load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
cfg.colormap = cmap;
cfg.layout = 'neuromag306cmb.lay';       %neuromag306all neuromag306mag neuromag306planar
cfg.zlim = [-0.2 0.2];

close all
for iloc=1:2 % left and right handers
    respavgpoolhand = squeeze(mean(respavg( LRtlocsubjects{iloc} ,:,:,:,:,:,:))); %average over subj
    for ises = 3 %1:3
        for itype = 1%:2 %1:2% %:ncond
            for ievent = 1:2 % off and on
                figure;    iplot=0; hold on
                
                freq.powspctrm = squeeze(respavgpoolhand(ises, :,frind,:,itype, ievent));
                ft_multiplotTFR(cfg, freq);
                title({sprintf('%s %slocked %s-side target', CON{itype,ievent}, trigger{1}, tlocleg{iloc}) ; sprintf('Sess: %s', sesnames{ises} )}); % 'FontWeight','bold'
                if SAV
                    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
                    outpath = fullfile(PREOUT, 'multiplot');
                    
                    warning off; mkdir(outpath); warning on
                    outfile = sprintf( '%s%smultiplot_%slocked_%s_%sfreq_%s_%s-sidetarget', outpath, filesep, trigger{itrg}, CON{itype,ievent}, analysistype{iat}, sesnames{ises}, tlocleg{iloc}  );
                    disp(outfile)
                    orient landscape
                    %                     orient portrait
                    %         print('-depsc', '-adobecs', '-painter', outfile) %, '-cmyk'
                    print('-dpdf', '-adobecs', '-painter', outfile) %, '-cmyk'
                    %                     print('-dpng', '-adobecs', '-painter', outfile) %, '-cmyk'
                    %         print('-dpdf', outfile)
                end;
            end
        end
    end
end







%% Plot single subj
%
% SAV=1;
% % close all
%
% for isoi = 1:length(sois)
%     soinameind = find(strcmp(  sois{isoi},  SOIN));
%     sensind = match_str(chlabel,ft_channelselection(SOI{soinameind},chlabel));
%
%     for ises = 1:2 % 1:3
%
%
%         for itype = 1 %:3 %:ncond
%             for ievent = 1:2 % off and on
%                 figure;    iplot=0;
%                 load('ColorMapTFRs');    colormap(ColorMapTFRs);
%
%                 for isub=1:length(SUBJ)
%
%                     iplot = iplot+1;
%                     dum = squeeze(mean(respavg(isub,ises, sensind,frind,:,itype, ievent))); %average over sensors
%                     % dum = squeeze(respavgpool(ises, 1,frind,:,itype,ievent));%take pooling if already averaged
%                     subplot(4,3,iplot);
%                     scale = ZLIM; %{inorm}(1,:);
%                     imagesc(taxis,faxis,dum,scale);
%                     title({sprintf('%s %s, [%g]', SUBJ{isub}, CON{itype,ievent}, ZLIM(2)) ; sprintf('%s Sess: %s', SOIN{soinameind}, sesnames{ises} )}); % 'FontWeight','bold'
%                     hold on
%                     yaxis = [CUTLO,CUTHI];
%                     plot([0,0],yaxis,'k',[0,0],yaxis,'k');
%                     if strfind(trigger{itrg}, 'stim')
%                         plot([medianrt(ises,itype,ievent),medianrt(ises,itype,ievent)],yaxis,'--k',[medianrt(ises,itype,ievent),medianrt(ises,itype,ievent)],yaxis,'--k'); %concat over sessions
%                     elseif strfind(trigger{itrg}, 'resp')
%                         plot([0-medianrt(ises,itype,ievent),0-medianrt(ises,itype,ievent)],yaxis,'--k',[0-medianrt(ises,itype,ievent),0-medianrt(ises,itype,ievent)],yaxis,'--k'); %concat over sessions
%                     end
%                     set(gca,'Box',BOX,'XTick',[-2:0.5:2],...    [-0.5,0,0.5,1,1.5,2,2.5]
%                         'YDir','normal','YTick',[0:10:CUTHI],...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
%                         'TickDir','out', 'FontSize', 8);
%                     if LAB %&& (iplot == 9 || iplot == 19)
%                         xlabel('Time (s)');
%                         ylabel('Frequency (Hz)');
%                         zlabel('Response (%)');
%                         colorbar
%                     end
%                 end %itype
%                 if SAV
%                     outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
%                     warning off; mkdir(outpath); warning on
%                     outfile = sprintf( '%s%sTFRpooling_%slock_singlesub_%s_%sfreq%s', outpath, filesep, trigger{itrg}, sois{isoi}, analysistype{iat}, CON{itype,ievent} );
%                     disp(outfile)
%                     orient landscape
%                     print('-dpsc2', outfile)
%                     print('-dpdf', outfile)
%                 end;
%
%             end % ievent
%         end  %isub
%     end %ises
%
% end


% %% plot
% figure; iplot=0;
%
% for itype=1:2 %1:3
%     for ievent=1:2 %test d and r against 0
%         iplot = iplot + 1;        subplot(2,2,iplot)
%         dum = squeeze(mean(respavgstat(:,frind,:,itype,ievent)));
%         imagesc(taxis, faxis, dum,scale);
% %         alpha(0.3)
%         hold on
%
%         pval = squeeze(stat{itype,ievent}.prob);
%
%         thrcluster = zeros(size(squeeze(stat{itype,ievent}.prob)));
%         thrcluster = thrcluster+0.3;
%         thrcluster(find(pval<0.05)) = 1;
%         alpha(thrcluster)
%
%         set(gca,'Box',BOX,'XTick',[-1 0 0.5 1],...    [-0.5,0,0.5,1,1.5,2,2.5]
%             'YDir','normal','YTick',YTICKS,...                %                 'YDir','normal','YTick',[15,20,50,100,150,200],...
%             'TickDir','out', 'FontSize', 8);
%         if LAB %&& (iplot == 9 || iplot == 19)
%             xlabel('Time (s)');
%             ylabel('Frequency (Hz)');
%             zlabel('Response (%)');
%             colorbar
%         end
%
% %
% %         for itoi = 1:61
% %             for ifrind = 1:length(frind)
% %                 if pval(ifrind, itoi) > THR
% %                     dum(ifrind,itoi) = NaN;
% %                     if groupstats
% %                         sigclusind(ifrind,itoi,ievent) = 0;
% %                         sigclus(:,ifrind,itoi,:,ievent) = NaN;
% %                     end
% %                 end
% %             end
%
% %         h = imagesc(taxis, faxis, dum, scale);
% %         set(h, 'alphadata', ~isnan(dum));
% %         hold off
%
%     end
% end


