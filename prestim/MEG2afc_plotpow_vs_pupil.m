% function MEG2afc_plotpow_vs_pupil()

PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/export/baseline/';

cd(PREIN)
% sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
sesnames = ['B' 'D' 'A' 'C'];

% example = dir('*.csv');
% textscan(example(1).name, ',', 0, 0, [0 0 1 13] )
% fprintf(fid, 'subj_idx, correct, RT, session_nr, run_nr, drug, simon, choice, diff, baseline_pupil, front_alpha, front_beta, occ_beta\n');
col=[];
col.subj = 1;
col.correct = 2;
col.RT = 3;
col.session_nr = 4;
col.run_nr = 5;
col.drug = 6;
col.simon = 7;
col.choice = 8;
col.diff = 9;
col.baseline_pupil = 10;
col.front_alpha = 11;
col.front_beta = 12;
col.occ_beta = 13;
% col.occ_gamma = 14;
col.pupil_bin0_4= 14;

csvdat = nan(21, 2,2, 1000, 15); % dimord subj, pharm, motor, trials, cols
container = [];
for isub = 1:22
    for ises = 1:4
        
%         % ATX = 1, placebo = 0
%         if ises==1, ipharm = 1;   imotor=1;
%         elseif ises==2, ipharm = 1;   imotor=2;
%         elseif ises==3, ipharm = 2;   imotor=1;
%         elseif ises==4, ipharm = 2;   imotor=2;
%         end

% vice versa!
        if ises==1, ipharm = 2;   imotor=1;
        elseif ises==2, ipharm = 2;   imotor=2;
        elseif ises==3, ipharm = 1;   imotor=1;
        elseif ises==4, ipharm = 1;   imotor=2;
        end

        csvname = sprintf('nk%d_%s_run*_trl_freq.csv', isub, sesnames(ises));
        csv = dir(csvname);
        if isempty(csv)
            warning('%s not found\n', csvname)
            continue
        end
        fprintf('%s\n', csv.name)
        temp = csvread(csv.name, 1, 0);
        ntrials = size(temp,1);
        csvdat( isub, ipharm, imotor, 1:ntrials, :) = temp;
        container = [container; temp];
            
    end
end
csvdat = csvdat([1:5,7:9,11:end], :, :, :, :); % remove 6 and 10
nsub = size(csvdat,1);

% export container
% write header and matrix to file
trl_outfile = 'data_alltrials.csv';
fid = fopen(fullfile(PREIN, trl_outfile), 'w') ;
fprintf(fid, 'subj_idx,response,rt,session_nr,run_nr,drug,simon,choice,diff,baseline_pupil,front_alpha,front_beta,occ_beta,pupilbin0_4,occ_alphabin\n');
fclose(fid);
dlmwrite(fullfile(PREIN, trl_outfile), container, '-append', 'precision', 10)


%% bin by power and plot pupil
meanpup = mean(csvdat(:,:,:,:,col.baseline_pupil),2); %drug
meanpup = nanmean(meanpup,3); % motor
meanpup = nanmean(meanpup,4); % trials

dat_bins = [];
for isub = 1:nsub
    for ipharm = 1:2        
        subjdat = squeeze(csvdat(isub, ipharm, :,:,:));
        subjdat = squeeze(cat(2, subjdat(1,:,:), subjdat(2,:,:))); % concat motor
        subjdat = subjdat(~isnan(subjdat(:,1)), :); %remove nans
        
        % do binning by meg
        prctl_size = 10;
        prctils = 0:prctl_size:100;
        nbins = length(prctils)-1;
        
        for icol = 1:3 % frontalpha col 11, frontbeta 12, occbeta 13
            colois = 11:13;
            colnames = {'occipital8-13Hz', 'occipital13-18Hz', 'occipital17-30Hz'} ;
            bin_edges = prctile(subjdat(:,colois(icol)), prctils );
            bin_ind = discretize(subjdat(:,colois(icol)), bin_edges, 1:nbins);
            
            %         pupdat = zscore(subjdat(:, 10));
            %         pupdat = ( subjdat(:, 10) - mean(subjdat(:, 10)) ) ./ mean(subjdat(:, 10)) * 100;
            pupdat = subjdat(:, 10);
            pupdat = (pupdat - meanpup(isub)) ./ meanpup(isub) * 100; % take all data as baseline
            %         pupdat = (pupdat - mean(pupdat)) ./ mean(pupdat);
            for ibin = 1:nbins
                dat_bins(isub, ipharm, ibin, icol) = mean( pupdat(bin_ind == ibin) ); % pupil col 10
            end
        end
    end
end

% %% plotting separate for drug and plac
% close all
% f = figure; hold on
% f.Position =[   2262         398         560         420];
% drug_conds = {'drug', 'placebo'};
% for ipharm=1:2
%     subplot(2,2,ipharm);
%     
%     dat = squeeze(dat_bins(:,ipharm,:));
% 
%     e = errorbar(1:nbins, mean(dat), std(dat)/sqrt(nsub), 'o', 'MarkerSize',10, 'MarkerEdgeColor','w','MarkerFaceColor','k');
% 	e.Color = 'k';
% %     lsline
%     axis square; box on
%     [r,p] = corr([1:nbins]', mean(dat)');
%     title(sprintf('%s r=%g, p=%g', drug_conds{ipharm}, r,p ))
%     xlim([0 nbins+1])
%     xlabel('pow')
%     ylabel('baseline pupil size (psc)')
% end


%% plotting drug and plac together
SAV = 1
close all
f = figure; hold on
f.Position =[  2262         181         700         600];
drug_conds = {'placebo', 'drug'};
% linecol = cbrewer('qual', 'Set1',3);
linecol = [0 0 1; 1 0 0];
r=[]; p=[];
for icol=1:3
    subplot(1,3,icol); hold on
    for ipharm=1:2
        
        dat = squeeze(dat_bins(:,ipharm,:,icol));
        
        e = errorbar(1:nbins, mean(dat), std(dat)/sqrt(nsub), ...
            'o', 'MarkerSize',12, 'MarkerEdgeColor','w','MarkerFaceColor',linecol(ipharm,:));
        e.Color = 'k';
        box on
        [r(ipharm),p(ipharm)] = corr([1:nbins]', mean(dat)');
        xlim([0 nbins+1])
        ylim([-8 8])
        xlabel(['MEG power ' upper(colnames{icol})])
        ylabel('baseline pupil size (psc)')
    end
    title(sprintf('%s r = %1.2f, p = %1.4f\n%s r = %1.2f, p = %1.4f', drug_conds{1}, r(1),p(1), drug_conds{2}, r(2),p(2) ))
    ref = refline(0,0);
    ref.Color = 'k';
    ref.LineStyle = '--';
    legend(drug_conds)
    legend boxoff
    ax=gca;
    ax.FontSize = 12;
end
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/plots';
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'megvspupil');
    
    mkdir(outpath)
    %     outfile = fullfile(outpath, sprintf( 'alphagroupvspower_%s_%s_%s_%s', respavg.sens.leg{isoi}, respavg.behav_conds{icond}, respavg.sdt_conds{istim, iresp}, respavg.freqband{iband}));
    outfile = fullfile(outpath, sprintf('binnedbymeg') );
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises


%% Correlate shift in meg and pupil across subjects
colois = 11:13;
% dimord subj, pharm, motor, trials, cols
dat = squeeze(cat(4, csvdat(:,:,1,:,:),  csvdat(:,:,2,:,:))); % concat motor
% dat(:,:,:,colois) = log(dat(:,:,:,colois));
dat = squeeze(nanmean(dat,3));  % avg over trials


% close all
f = figure; hold on
f.Position = [2278         651         800         400];
r=[]; p=[];

for icol = 1:3
    subplot(1,3,icol); hold on
    axis square; box on

%     megdat = (dat(:,1,colois(icol)) - dat(:,2,colois(icol)))  ./ dat(:,2,colois(icol)) * 100;
%         megdat = (dat(:,1,colois(icol)) - dat(:,2,colois(icol))) ;
        megdat = dat(:,2,colois(icol));
    
%         pupdat = (dat(:,1,10) - dat(:,2,10))  ./ dat(:,2,10) * 100;
%         pupdat = (dat(:,1,10) - dat(:,2,10));
        pupdat = dat(:,2,10);
    scatter( megdat, pupdat)
    [r(icol),p(icol)] = corr(megdat, pupdat);
end

%% plot RT and acc for drug and plac
% dimord subj, pharm, motor, trials, cols
SAV=1;
% behavdat = squeeze(cat(4, csvdat(:,:,1,:,:), csvdat(:,:,2,:,:))); % concat motor
% accdat = behavdat(:,:,:,col.correct);
% rtdat = behavdat(:,:,:,col.RT);
% behavdat_all = [];
% for isub=1:nsub
%     for idrug=1:2
%         accdat_subj = squeeze(accdat(isub, idrug,:));
%         accdat_subj = accdat_subj(~isnan(accdat_subj));
%         behavdat_all(isub, idrug, 1) = length(find(accdat_subj)) / length(accdat_subj) * 100;
%         
%         rtdat_subj = squeeze(rtdat(isub, idrug,:));
%         rtdat_subj = rtdat_subj(~isnan(rtdat_subj));
%         behavdat_all(isub, idrug, 2) = mean(rtdat_subj);
%         behavdat_all(isub, idrug, 3) = std(rtdat_subj);
%         
%     end
% end

behavdat = squeeze(cat(4, csvdat(:,:,1,:,:), csvdat(:,:,2,:,:))); % concat motor
% accdat = behavdat(:,:,:,col.correct);
% rtdat = behavdat(:,:,:,col.RT);
behavdat_all = [];
for isub=1:nsub
    for idrug=1:2
        accdat_subj = squeeze(behavdat(isub, idrug,:,:));
        accdat_subj = accdat_subj(~isnan(accdat_subj(:,1)),:);
        for idiff = 1:2
            trl_ind = accdat_subj(:,col.diff) == idiff-1; % get easy and hard
            
            behavdat_all(isub, idrug, mod(idiff,2)+1, 1) = ...
                length( find(accdat_subj(:,col.correct)==1 & trl_ind) ) ...
                / length(find(trl_ind)) * 100;
            
            behavdat_all(isub, idrug, mod(idiff,2)+1, 2) = mean( accdat_subj(trl_ind, col.RT) );
            behavdat_all(isub, idrug, mod(idiff,2)+1, 3) = std( accdat_subj(trl_ind, col.RT) );
            
            
%             rtdat_subj = squeeze(rtdat(isub, idrug,:));
%             rtdat_subj = rtdat_subj(~isnan(rtdat_subj));
%             behavdat_all(isub, idrug, 2) = mean(rtdat_subj);
%             behavdat_all(isub, idrug, 3) = std(rtdat_subj);
        end
    end
end

% behavdat_all = (behavdat_all - behavdat_all(:,1,:)) ./  behavdat_all(:,1,:);
% plot 
close all
measleg = {'Accuracy', 'Reaction time', 'RT variability' };
behavleg = {'Correct (%)', 'RT (s)', 'RT (sd)'};
diff_leg = {'Weak', 'Strong'};
f = figure;
f.Position = [ 2042         534         750         250];
for imeas=1:3
    subplot(1,3,imeas);
%     p = permtest(behavdat_all(:,1,:,imeas), behavdat_all(:,2,:,imeas));
    b=barweb( squeeze(mean(behavdat_all(:,:,:,imeas)))', squeeze(std(behavdat_all(:,:,:,imeas)))' / sqrt(nsub), 0.75, ...
        diff_leg, measleg{imeas}, 'Sensory evidence', behavleg{imeas} ); % sprintf('p = %g', p) 
    b.bars(1).FaceColor = [0 0 1];
    b.bars(1).EdgeColor = [0 0 1];
    b.bars(2).FaceColor = [1 0 0];
    b.bars(2).EdgeColor = [1 0 0];
    ax=gca;
    ax.FontSize=16;
    if imeas==1
        ylim([70 85])
        ax.YTick = [0:5:100];
    elseif imeas==2
        ylim([0.85 1.05])
        ax.YTick = [0:0.05:2];
    elseif imeas==3
        ylim([0.25 0.4])
        ax.YTick = [0:0.05:2];
        legend({'Placebo', 'ATX'}); legend boxoff
    end
    
end
if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'behavior_modelfree'));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf') %'-png',  '-pdf', , '-transparent'
    print(outfile, '-dpdf') %'-png',  '-pdf', , '-transparent'
%     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises


%% plot RT and acc for drug against plac high and low alpha
% dimord subj, pharm, motor, trials, cols
behavdat = squeeze(cat(4, csvdat(:,:,1,:,:), csvdat(:,:,2,:,:))); % concat motor

colnames = {'front_alpha' 'front_beta' 'occ_beta' 'baseline_pupil'};
coloi = 4; % bin this

accdat = []; rtdat = [];
for isub=1:nsub
    for idrug=1:2
        dat = squeeze(behavdat(isub,idrug,:,:));
%         dat = dat(~isnan(dat(:,1)),:);

        bin_edges = prctile(dat(:,col.(colnames{coloi})), 0:20:100 );
        nbins = length(bin_edges)-1;
        bin_alpha = discretize(dat(:,col.(colnames{coloi})), bin_edges, 1:nbins);
        for ibin= 1:nbins
            %accuracy:
%             trialsoi = dat(bin_alpha == ibin, col.correct);
            trialsoi = dat(bin_alpha == ibin & dat(:,col.diff) > -1 , col.correct);
            accdat(isub, ibin, idrug) =  length(find(trialsoi == 0)) / length(trialsoi); % get correctness for low alpha trials
            
            % RT
%             trialsoi = dat(bin_alpha == ibin, col.RT);
            trialsoi = dat(bin_alpha == ibin & dat(:,col.diff) > -1 , col.RT); % & trialsoi == 1
            rtdat(isub, ibin, idrug) =  mean(trialsoi); % get rt for low alpha trials

%             % gamma modulation
%             trialsoi = dat(bin_alpha == ibin & dat(:,col.diff) > -1 , col.occ_gamma);
%             rtdat(isub, ibin, idrug) =  mean(trialsoi); % get rt for low alpha trials
%             
            
        end            
       
    end
end

% plot
SAV = 1
close all
f = figure; hold on
iplot=0;
f.Position =[  2262         181         700         600];
drug_conds = {'drug', 'placebo'}; %% CHECK!!
linecol = cbrewer('qual', 'Set1',3);
for idrug = 1:2
    iplot = iplot + 1;
    subplot(2,2,iplot)
    barweb( mean(accdat(:,:,idrug)), std(accdat(:,:,idrug)) / sqrt(nsub), 0.75  )
%     ylim([0.75 0.8])
    [h,p] = ttest(accdat(:,1,idrug), accdat(:,2,idrug) );
    title(drug_conds{idrug})
    ylabel('Accuracy')
    % RT:
    iplot = iplot + 1;
    subplot(2,2,iplot)
    barweb( mean(rtdat(:,:,idrug)), std(rtdat(:,:,idrug)) / sqrt(nsub), 0.75 )
    ylim([0.85 1])
    [h,p] = ttest(rtdat(:,1,idrug), rtdat(:,2,idrug) );
    title(drug_conds{idrug})
    ylabel('RT (s)')
end
% legend({'low', 'med', 'high'})


%% check if occipital variance is bigger for drug than placebo
%  csvdat( isub, ipharm, imotor, 1:ntrials, :)
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/plots';

close all
SAV=1;

alphadat = squeeze(cat(4, csvdat(:,:,1,:,col.front_alpha), csvdat(:,:,2,:,col.front_alpha))); % concat motor
alphadat = log(alphadat);
alphavar = nanvar(alphadat, 0, 3); % non-normalized variance
% normalized into psc
% alphavar(:,:,2) = (alphavar - mean(alphavar,2)) ./ mean(alphavar,2) * 100;
% wrt placebo
alphavar(:,:,2) = (alphavar - alphavar(:,1)) ./ alphavar(:,1) * 100;

f = figure;
f.Position = [ 943   679   500   326]

for inorm = 1:2 % normalized
    subplot(1,2,inorm)
        if inorm==1
            b=barweb(mean(alphavar(:,:,inorm)), std(alphavar(:,:,inorm)) / sqrt(nsub), 0.75);
            %     b=barweb(mean(alphavar(:,2,inorm)), std(alphavar(:,2,inorm)) / sqrt(nsub), 0.75);
            b.bars(1).FaceColor = [0 0 1];
            b.bars(1).EdgeColor = [0 0 1];
            b.bars(2).FaceColor = [1 0 0];
            b.bars(2).EdgeColor = [1 0 0];
        else
            b=barweb(mean(alphavar(:,2,inorm)), std(alphavar(:,2,inorm)) / sqrt(nsub), 0.3);
            %     b=barweb(mean(alphavar(:,2,inorm)), std(alphavar(:,2,inorm)) / sqrt(nsub), 0.75);
            b.bars(1).FaceColor = [1 0 0];
            b.bars(1).EdgeColor = [1 0 0];
        end
    
    p =permtest(alphavar(:,1,inorm), alphavar(:,2,inorm))
    title(sprintf('Prestim alpha variability\np = %1.4f', p))
    ax=gca;
    ax.FontSize = 16;
    if inorm==1
        legend({'Placebo', 'ATX'}); legend boxoff
        ylabel('Normalized variance')
        ylim([0.4 0.85])
        ax.YTick = [0:0.1:1];
    else
%         ylabel('Normalized variance (% wrt grand average)')
        ylabel(sprintf('Normalized variance\n(%% from placebo)'))
        ylim([-10 7.5])
        ax.YTick = [-10:5:10];

    end
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'alphavar'));
    %                 end
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf') %'-png',  '-pdf',
%     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises


%% Determine two sub-groups: low and high baseline alpha/pupil.

%  csvdat( isub, ipharm, imotor, 1:ntrials, :)
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/plots';

% close all
SAV=1;

alphadat = squeeze(cat(4, csvdat(:,:,1,:,col.front_alpha), csvdat(:,:,2,:,col.front_alpha))); % concat motor
alphadat = log(alphadat);

alphadat_plac = squeeze(nanmean(alphadat(:,2,:),3));
alphadat_ATX = squeeze(nanmean(alphadat(:,1,:),3));
figure;
histogram(alphadat_plac, 20)
xlabel('Log alpha power')
ylabel('N subjects')

% lowgroup = SUBJ(alphadat_plac < median(alphadat_plac));
% highgroup = SUBJ(alphadat_plac >= median(alphadat_plac));
lowgroup = alphadat_plac < median(alphadat_plac);
highgroup = alphadat_plac >= median(alphadat_plac);

% p = permtest(alphadat_plac, alphadat_ATX)
% p = permtest((alphadat_ATX- alphadat_plac) ./ alphadat_plac)

% make mat for JW with plac en atx 1 en 2
% csvdat dimord csvdat( isub, ipharm, imotor, 1:ntrials, :) 
jwdat = squeeze(nanmean(csvdat(:,:,1:2,:,col.front_alpha), 4)); % dimord subj drug ses
jwdat = log(jwdat);
outmat = [jwdat(:,:,1) jwdat(:,:,2)]; % plac 1drug1  plac2 drug2 
% copy paste in xcel

figure; scatter(outmat(:,1), outmat(:,3))
corr(outmat(:,1), outmat(:,3))

%% correlate placebo 1 vs 2: baseline state stable trait?
SAV = 1;
alphadat_plac = squeeze( csvdat(:,1,:,:,col.front_alpha));  % take plac data ses 1 and 2
alphadat_plac = log(alphadat_plac);

% subj_oi = [1,3:19];
% subj_oi = [3:19];
subj_oi = [1:19];
alphadat_plac = nanmean(alphadat_plac,3);
alphadat_plac = alphadat_plac(subj_oi,:); % NK2 lacks plac sess 2

close all
figure; hold on; axis square; box on
s = scatter(alphadat_plac(:,1), alphadat_plac(:,2), 'filled');
l = lsline;
l.Color = 'k'

text(alphadat_plac(:,1), alphadat_plac(:,2), respavg.SUBJ(subj_oi))

s.MarkerFaceColor = 'k';
s.MarkerEdgeColor = 'w';
s.LineWidth = 2;
s.SizeData = 450
[p,r] = permuteCorr(alphadat_plac(:,1), alphadat_plac(:,2))
title(sprintf('Placebo 1 vs 2: r = %1.2f, p = %g', r,p))
xlabel('Placebo 1 occipital alpha power (log)')
ylabel('Placebo 2 occipital alpha power (log)')
xlim([-63 -58])
ylim([-63 -58])

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'plac1vs2corr'));
    %                 end
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf') %'-png',  '-pdf',
%     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

