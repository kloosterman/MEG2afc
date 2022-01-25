% read in para's from HDDM csv and plot
% MEG2afc_plot_hddmparas

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5bins_pupil.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_10bins_pupil.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_10bins_occalpha.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3binspupil_newbinning.csv';

% csvfile =
% '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3binsalpha_newbinning.csv';% % works well
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_alpha_3bins33.csv'; %5quantiles

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_alpha_3bins402040.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5binsalpha_newbinning.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_1bin.csv';

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5binsalpha_newbins5.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5binsalpha_newbins5_2.csv';

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_nobinning_incdiff.csv'; % used in SPR talk

csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_drug_diff2.csv'; % < 0.2 rt dropped, 5 rt bins

temp = csvread(csvfile, 1,1);
ddmleg = {'a_plac', 'a_drug', 'nd_plac', 'nd_drug', ...
    ... %     'v_plac_easy', 'v_drug_easy', 'v_plac_hard', 'v_drug_hard' }
    'v_atx_easy', 'v_plac_easy', 'v_atx_hard', 'v_plac_hard' }
paraleg = {'Boundary separation', 'Non-decision time', 'Drift rate'};

% dims now a drug 1-5, a plac 1-5
nbins = 1;
ndiff = 2;
nsub = 19;
ndrug = 2;
npara = 3;
% ddmdat = nan(nsub, nbins, ndrug, ndiff, npara);
ddmdat = temp;

% ddmdat = ddmdat(lowgroup,:,:,:);
% ddmdat = ddmdat(highgroup,:,:,:);

bodyweights =  [  74    90    81    55    55    62    66    68    56    72    71   NaN    60    72    63    63    67    84    64 ]';

ddmdat = [ddmdat bodyweights]

%% plotting bars for drug and plac, 3 DDM params
SAV = 1;
% close all
f = figure;
f.Position = [2043         535         749         250];
subplot(1,3,2)
dum = ddmdat(:,1:2);
% dum = (ddmdat(:,1:2) - ddmdat(:,1)) ./ ddmdat(:,1) * 100;
b = barweb(mean(dum), std(dum) / sqrt(nsub), 0.75 ); % bound
b.bars(1).FaceColor = [0 0 1];b.bars(1).EdgeColor = [0 0 1]; b.bars(2).FaceColor = [1 0 0]; b.bars(2).EdgeColor = [1 0 0];
p = permtest(ddmdat(:,1), ddmdat(:,2));
title(sprintf('%s\np = %g', paraleg{1}, p))
% ylim([1.25 1.55]);
ax=gca;
ax.XTickLabel = '';
ax.FontSize = 14;

subplot(1,3,3)
dum = ddmdat(:,3:4);
% dum = (ddmdat(:,3:4) - ddmdat(:,3)) ./ ddmdat(:,3) * 100;
b = barweb(mean(dum), std(dum) / sqrt(nsub), 0.75 ); % nd
b.bars(1).FaceColor = [0 0 1];b.bars(1).EdgeColor = [0 0 1]; b.bars(2).FaceColor = [1 0 0]; b.bars(2).EdgeColor = [1 0 0];p = permtest(ddmdat(:,3), ddmdat(:,4));
title(sprintf('%s\np = %g', paraleg{2}, p))
ylim([0.43 0.53]);
ax=gca;
ax.XTickLabel = '';
ax.FontSize = 14;

subplot(1,3,1)
dum=[];
% % dum(:,:,1) = ddmdat(:, [5 7]); % 'v_plac_easy', 'v_drug_easy', so dimord= subj diff drug 
% % dum(:,:,2) = ddmdat(:, [6 8]); % 'v_plac_hard', 'v_drug_hard'
% % dum(:,:,1) = ddmdat(:, [7 5]); % 'v_plac_easy', 'v_drug_easy', so dimord= subj diff drug 
% % dum(:,:,2) = ddmdat(:, [8 6]); % 'v_plac_hard', 'v_drug_hard'
% dum(:,:,1) = ddmdat(:, [6 5]); % 'v_plac_easy', 'v_drug_easy', so dimord= subj diff drug 
% dum(:,:,2) = ddmdat(:, [8 7]); % 'v_plac_hard', 'v_drug_hard'
% % dum = (dum - dum(:,:,1)) ./ dum(:,:,1) * 100;

dum(:,1,2) = ddmdat(:, 5); % 'v_drug_easy' so dimord= subj diff drug 
dum(:,1,1) = ddmdat(:, 6); % 'v_plac_easy' so dimord= subj diff drug 
dum(:,2,2) = ddmdat(:, 7); % 'v_drug_hard' so dimord= subj diff drug 
dum(:,2,1) = ddmdat(:, 8); % 'v_plac_hard' so dimord= subj diff drug 

p1 = permtest(dum(:,1,1), dum(:,1,2)); 
p2 = permtest(dum(:,2,1), dum(:,2,2));

b = barweb( squeeze(mean(dum)), squeeze(std(dum)) / sqrt(nsub), 0.75 ); % drift
b.bars(1).FaceColor = [0 0 1];b.bars(1).EdgeColor = [0 0 1]; b.bars(2).FaceColor = [1 0 0]; b.bars(2).EdgeColor = [1 0 0];title(sprintf('%s\np = %g, p = %g', paraleg{3}, p1, p2))
ylim([0.7 1.4]);
legend({ 'Placebo', 'ATX'}, 'Location', 'NorthEast'); legend boxoff
ax=gca;
ax.XTickLabel = {'Easy' , 'Hard'};
ax.XLabel.String = 'Sensory evidence';
ax.FontSize = 14;
ylabel('Parameter estimate (a.u.)')

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'ddm_nobins'));
    %                 end
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% correlate para's with weight
%%'a_plac'    'a_drug'    'nd_plac'    'nd_drug'    'v_plac_easy'    'v_plac_hard'    'v_drug_easy'    'v_drug_hard'
close all
dum = [];
% % dum(:,:,1) = ddmdat(:, [7 5]); % 'v_drug_easy', 'v_plac_easy',  so dimord= subj drug diff
% % dum(:,:,2) = ddmdat(:, [8 6]); % 'v_drug_hard' 'v_plac_hard', 
% dum(:,:,1) = ddmdat(:, [5 6]); %  'v_plac_easy', 'v_drug_easy',  so dimord= subj drug diff
% dum(:,:,2) = ddmdat(:, [7 8]); % 'v_plac_hard', 'v_drug_hard' 

dum(:,1,2) = ddmdat(:, 5); % 'v_drug_easy' so dimord= subj diff drug 
dum(:,1,1) = ddmdat(:, 6); % 'v_plac_easy' so dimord= subj diff drug 
dum(:,2,2) = ddmdat(:, 7); % 'v_drug_hard' so dimord= subj diff drug 
dum(:,2,1) = ddmdat(:, 8); % 'v_plac_hard' so dimord= subj diff drug 

dum(:,3,:) = mean(dum,2); % avg over drug
dum(:,:,3) = mean(dum,3); % avg over diff

dum(:,4,:) = dum(:,2,:) - dum(:,1,:); % drug-plac drift
dum(:,:,4) = dum(:,:,1) - dum(:,:,2); % easy-hard drift

% % drop high drift subj1 and FJ
% dum = dum([2:11, 13:end],:,:); 
% weight = ddmdat([2:11, 13:end], 9);

dum = dum([1:11, 13:end],:,:); 
weight = ddmdat([1:11, 13:end], 9);

diffleg = { 'Easy' 'Hard' 'diffavg' 'Easy-Hard'};
drugleg = {'Plac' 'ATX'  'drugavg' 'ATX-plac'};
f = figure; iplot = 0;
f.Position = [          2132          60         1100        1100];
for idrug = 1:4
    for idiff = 1:4
        iplot = iplot + 1;
        subplot(4,4,iplot)
        sc = scatter( weight, dum(:,idiff,idrug), 'o', 'filled', 'k');
        sc.SizeData = 70;
        axis square; box on
        [r,p] = corr( weight, dum(:,idiff,idrug), 'type', 'Spearman');
        title(sprintf('%s %s r = %g', drugleg{idrug}, diffleg{idiff},  r ))
        if p < 0.05; lsline; end
        ax = gca;
        ax.FontSize = 12;
        xlabel('Body weight (kg)')
        ylabel('Drift rate')
    end
end

%% plotting bars for ATX wrt plac, 3 DDM params
SAV = 1;
% close all
f = figure;
f.Position = [ 2150         642         600         400];

dum=[];
% dum(:,:,1) = ddmdat(:, [5 7]); % 'v_plac_easy', 'v_drug_easy', so dimord= subj diff drug 
% dum(:,:,2) = ddmdat(:, [6 8]); % 'v_plac_hard', 'v_drug_hard'
dum(:,:,1) = ddmdat(:, [7 5]); % 'v_plac_easy', 'v_drug_easy', so dimord= subj diff drug 
dum(:,:,2) = ddmdat(:, [8 6]); % 'v_plac_hard', 'v_drug_hard'
% dum = (dum - dum(:,:,1)) ./ dum(:,:,1) * 100;
dum = (dum - dum(:,:,1)) ;

ddmpsc = [];
ddmpsc(:,1:2) = dum(:,:,2); % drift

% dum = ddmdat(:,1:2);
% dum = (ddmdat(:,1:2) - ddmdat(:,1)) ./ ddmdat(:,1) * 100;
dum = (ddmdat(:,1:2) - ddmdat(:,1)) ;

ddmpsc(:,3) = dum(:,2);

% NDT
% dum = (ddmdat(:,3:4) - ddmdat(:,3)) ./ ddmdat(:,3) * 100;
dum = (ddmdat(:,3:4) - ddmdat(:,3)) ;
ddmpsc(:,4) = dum(:,2);

% stats vs 0
p=[];
for i=1:4
    p(i) = permtest(ddmpsc(:,i));
end

b = barweb(mean(ddmpsc), std(ddmpsc) / sqrt(nsub), 0.75 ); % bound
b.bars(1).FaceColor = [1 0.5 0 ]; b.bars(1).EdgeColor = b.bars(1).FaceColor;
b.bars(2).FaceColor = [0.9 0.1 0 ]; b.bars(2).EdgeColor = b.bars(2).FaceColor;
b.bars(3).EdgeColor = b.bars(3).FaceColor;
b.bars(4).EdgeColor = b.bars(4).FaceColor;
% p = permtest(ddmdat(:,1), ddmdat(:,2));
% title(sprintf('%s\np = %g', paraleg{1}, p))
% ylim([1.25 1.55]);
ax=gca;
ax.XTickLabel = '';
ax.FontSize = 14;
ylabel('Parameter change from placebo (%)')
legend({'Drift, weak evidence', 'Drift, strong evidence', 'Boundary separation', 'Non-decision time'}, 'Location', 'eastoutside')
legend boxoff
title(p)

% subplot(1,3,3)
% % dum = ddmdat(:,3:4);
% b = barweb(mean(dum), std(dum) / sqrt(nsub), 0.75 ); % nd
% b.bars(1).FaceColor = [0 0 1];b.bars(1).EdgeColor = [0 0 1]; b.bars(2).FaceColor = [1 0 0]; b.bars(2).EdgeColor = [1 0 0];p = permtest(ddmdat(:,3), ddmdat(:,4));
% title(sprintf('%s\np = %g', paraleg{2}, p))
% % ylim([0.43 0.53]);
% ax=gca;
% ax.XTickLabel = '';
% ax.FontSize = 14;
% 
% subplot(1,3,1)
% 
% p1 = permtest(dum(:,1,1), dum(:,1,2));
% p2 = permtest(dum(:,2,1), dum(:,2,2));
% 
% b = barweb( squeeze(mean(dum)), squeeze(std(dum)) / sqrt(nsub), 0.75 ); % drift
% b.bars(1).FaceColor = [0 0 1];b.bars(1).EdgeColor = [0 0 1]; b.bars(2).FaceColor = [1 0 0]; b.bars(2).EdgeColor = [1 0 0];title(sprintf('%s\np = %g, p = %g', paraleg{3}, p1, p2))
% ylim([0.7 1.25]);
% legend({'Placebo', 'ATX'}, 'Location', 'NorthWest'); legend boxoff
% ax=gca;
% ax.XTickLabel = {'Weak' , 'Strong'};
% ax.XLabel.String = 'Sensory evidence';
% ax.FontSize = 14;

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'ddm_nobins_fromplac'));
    %                 end
    disp(outfile)
%     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises


%%

% ddmdat = reshape(temp, nsub,nbins,ndrug, ndiff, npara); %dimord subj bin# drug para
% ddmdat(:,:,3,:) = mean(ddmdat,3);
% nsub = size(ddmdat,1);

% remove subjects?
% ddmdat = ddmdat(2:end,:,:,:); % drop NK1, extreme drift
% ddmdat = ddmdat([2:7, 9:end],:,:,:); % drop NK1, extreme drift

% % normalize drug and plac per bin
% for ibin=1:nbins
%     ddmdat2(:,ibin,1:2,:) = (ddmdat(:,ibin,1:2,:) - mean(ddmdat(:,ibin,3,:),2)); ...
% %         ./ mean(ddmdat(:,ibin,3,:),2) ;
% end

% normalize data per subject by dividing by largest drift across bins (Rajagovindan)
nsub = size(ddmdat,1);
ddmdat2 = [];
for isub=1:nsub
    for ipara=1:3 % drift etc
        %         ddmdat2(isub,:,:,ipara) = ddmdat(isub,:,:,ipara) ./ max(ddmdat(isub,:,2,ipara)) ;
        ddmdat2(isub,:,:,ipara) = (ddmdat(isub,:,:,ipara) - max(ddmdat(isub,:,1,ipara))) ./ max(ddmdat(isub,:,1,ipara)) * 100 ;
    end
end

% psc wrt ...
ddmdat_norm = []; range =[];
for idrug = 1:2
    for ipara = 1:3
        ddmdat_norm(:,:,idrug,ipara) = (ddmdat(:,:,idrug,ipara) - mean(ddmdat(:,:,1,ipara),2)) ...
            ./ mean(ddmdat(:,:,1,ipara),2) *100;
        
        %         range(:,idrug, ipara) = max(ddmdat_norm(:,:,idrug,ipara), [],2) - min(ddmdat_norm(:,:,idrug,ipara), [],2);
        
    end
end

PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/plots';
% plotting
SAV=1;

% addpath( '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/stats')
paraleg = {'Boundary separation', 'Non-decision time', 'Drift rate'};
drug_colors = {'b', 'r'};
close all
f = figure;
f.Position = [ 1998         466        1069         358];
for ipara = 1:3
    subplot(1,3,ipara); hold on; box on; axis square;
    for idrug = 1:2 % 1:2
        
        dum = ddmdat_norm(:,:,idrug,ipara);
        
        %                 dum = ddmdat2(:,:,idrug,ipara);
        %         dum = ddmdat(:,:,idrug,ipara);
        
        
        e = errorbar(1:nbins, mean(dum), std(dum) / sqrt(size(dum,1)), 'o', ...
            'MarkerFaceColor', drug_colors{idrug});
        %         e = errorbar(1:nbins, mean(ddmdat(:,:,idrug,ipara)), std(ddmdat(:,:,idrug,ipara)) / sqrt(size(ddmdat,1)), 'o', ...
        %              'MarkerFaceColor', drug_colors{idrug});
        e.Color = drug_colors{idrug};
        e.MarkerEdgeColor = [1 1 1];
        e.MarkerSize = 12;
        title(paraleg{ipara})
        %         legend boxoff
        xlim([0 nbins+1])
        ax=gca;
        ax.XTick = 1:nbins;
        ax.XTickLabel = {'Low' 'Medium' 'High'};
        ax.FontSize = 12;
        xlabel('Prestimulus occipital alpha power')
        %         xlabel('Baseline pupil bin')
        ylabel(sprintf('Parameter estimate\n (%% wrt grand average)'))
        %         ylabel(sprintf('Parameter estimate\n (%% wrt max placebo drift)'))
        %         ylabel(sprintf('Parameter estimate\n (raw DDM output)'))
        
        fitorder = 1; % linear
        if ipara == 3
            fitorder = 2; % quadratic
        end
        p = polyfit(1:nbins, mean(dum), fitorder);
        stepsize = 100;
        y = polyval(p, 1:nbins/stepsize:nbins);
        pl = plot(1:nbins/stepsize:nbins, y);
        pl.Color = e.MarkerFaceColor;
        %         pl.DisplayName = respavg.behav_conds{icond};
        pl.LineWidth = 1;
    end
    if ipara==2
        legend({ 'placebo', 'placebo fit', 'drug', 'drug fit'}, 'Location', 'northwest')
        legend boxoff
    end
    %     r = refline(0,0);
    %     r.Color = 'k';
    %     r.LineStyle = '--';
    
end
% test drug vs plac
% permtest(ddmdat_norm(:,1,1,3), ddmdat_norm(:,1,2,3))
% permtest(ddmdat_norm(:,2,1,3), ddmdat_norm(:,2,2,3))
% permtest(ddmdat_norm(:,3,1,3), ddmdat_norm(:,3,2,3))

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'binnedgain'));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% plot bars of non-binned DDM results
SAV=1;

% addpath( '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/stats')
paraleg = {'Boundary separation', 'Non-decision time', 'Drift rate'};
drug_colors = {'b', 'r'};
close all
f = figure;
f.Position = [ 1998         466        1069         358];
for ipara = 1:3
    subplot(1,3,ipara);
    
    dum = squeeze(ddmdat_norm(:,:,1:2,ipara));
    %     dum = squeeze(ddmdat(:,:,1:2,ipara));
    
    p = permtest(dum(:,1), dum(:,2) );
    b=barweb(mean(dum), std(dum) / sqrt(nsub), 0.75, ...
        [], sprintf('p = %g', p), paraleg{ipara} );
    b.bars(1).FaceColor = [0 0 1];
    b.bars(1).EdgeColor = [0 0 1];
    b.bars(2).FaceColor = [1 0 0];
    b.bars(2).EdgeColor = [1 0 0];
    ax=gca;
    ax.FontSize=16;
    %     if imeas==1
    %         ylim([50 100])
    %         ax.YTick = [0:25:100];
    %         legend({'placebo', 'drug'}); legend boxoff
    %     elseif imeas==2
    %         ylim([0.5 1.25])
    %         ax.YTick = [0:0.25:2];
    %     elseif imeas==3
    %         ylim([0 0.5])
    %         ax.YTick = [0:0.1:2];
    %     end
end




%% polyfit 2nd order on single subjects, get quadratic fit (gain)
% close all
SAV=1;
bin_fits = [];
for isub = 1:nsub
    for ipara = 1:3
        for idrug = 1:2
            dum = (ddmdat(isub,:,idrug,ipara) - mean(ddmdat(isub,:,3,ipara),2)) ...
                ./ mean(ddmdat(isub,:,3,ipara),2);
            %             dum = ddmdat(isub,:,idrug,ipara);
            bin_fits(isub,idrug,ipara,:) = polyfit(1:nbins, squeeze(dum), 2);
            
            %             bin_fits(isub,idrug,ipara,:) = polyfit(1:nbins, squeeze(ddmdat(isub,:,idrug,ipara)), 2);
        end
    end
end
% bin_fits = (bin_fits - mean(bin_fits, 2)) ./ mean(bin_fits, 2) * 100;

% plot
f = figure; iplot=0;
coeff_leg = {'Quadratic' 'Linear' 'Intercept'};
f.Position = [        2324          74         600         250];
for ipara = 3
    for ic = 1:3 % fit coeff
        iplot = iplot+1;
        subplot(1,3,iplot);
        b = barweb( mean(bin_fits(:,:,ipara,ic)), std(bin_fits(:,:,ipara,ic))/sqrt(nsub), 0.75 );
        b.bars(1).FaceColor = [0 0 1];
        b.bars(1).EdgeColor = [0 0 1];
        b.bars(2).FaceColor = [1 0 0];
        b.bars(2).EdgeColor = [1 0 0];
        
        %         [h,p] = ttest(bin_fits(:,1,ipara,ic), bin_fits(:,2,ipara,ic))
        p = permtest(bin_fits(:,1,ipara,ic), bin_fits(:,2,ipara,ic));
        title(sprintf('%s, p = %1.4f', coeff_leg{ic}, p))
        if ic==1
            legend({ 'placebo', 'drug'}, 'Location', 'northwest')
            legend boxoff
            ylabel('Coeff. weight')
        end
        if ic==2
        end
    end
end

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'polyfits'));
    %                 end
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    %     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

%% partial out the linear from the quadratic term
ipara = 3; resid=[];

for idrug=1:2
    quadrdat = bin_fits(:,idrug,ipara,1);
    lineardat = bin_fits(:,idrug,ipara,2);
    
    [~,~,resid(:,idrug)] = regress(quadrdat, lineardat );
    figure;
    scatter(quadrdat, lineardat)
end
f=figure;
barweb(mean(resid) , std(resid)/sqrt(nsub))




