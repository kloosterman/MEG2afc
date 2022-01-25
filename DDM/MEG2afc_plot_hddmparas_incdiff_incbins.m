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

csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3bins_incdiff.csv';

temp = csvread(csvfile, 1,1);
% ddmleg = {'a_plac', 'a_drug', 'nd_plac', 'nd_drug', ...
%     'v_plac_easy', 'v_plac_hard', 'v_drug_easy', 'v_drug_hard' }
paraleg = {'Boundary separation', 'Non-decision time', 'Drift rate'};

% dims now a drug 1-5, a plac 1-5
nbins = 3;
ndiff = 2;
nsub = 19;
ndrug = 2;
npara = 3;
% ddmdat = nan(nsub, nbins, ndrug, ndiff, npara);
ddmdat = temp; % dimord of cols: diff drug bin(3)

%% plotting
SAV = 1;
close all
f = figure;
f.Position = [  2104         762        1276         245];
% subplot(1,3,1)
% dum = ddmdat(:,1:2);
% % dum = (ddmdat(:,1:2) - ddmdat(:,1)) ./ ddmdat(:,1) * 100;
% b = barweb(mean(dum), std(dum) / sqrt(nsub), 0.75 ); % bound
% b.bars(1).FaceColor = [0 0 1];
% b.bars(1).EdgeColor = [0 0 1];
% b.bars(2).FaceColor = [1 0 0];
% b.bars(2).EdgeColor = [1 0 0];
% ylabel('Parameter estimate (a.u.)')
% p = permtest(ddmdat(:,1), ddmdat(:,2));
% title(sprintf('%s\np = %g', paraleg{1}, p))
% % ylim([1.25 1.6]);
% 
% subplot(1,3,2)
% dum = ddmdat(:,3:4);
% % dum = (ddmdat(:,3:4) - ddmdat(:,3)) ./ ddmdat(:,3) * 100;
% b = barweb(mean(dum), std(dum) / sqrt(nsub), 0.75 ); % nd
% b.bars(1).FaceColor = [0 0 1];
% b.bars(1).EdgeColor = [0 0 1];
% b.bars(2).FaceColor = [1 0 0];
% b.bars(2).EdgeColor = [1 0 0];
% p = permtest(ddmdat(:,3), ddmdat(:,4));
% title(sprintf('%s\np = %g', paraleg{2}, p))
% % ylim([0.4 0.55]);

subplot(1,4,3)
dum=[];
dum(:,1:3,1,1) = ddmdat(:, 13:15); % bin1-3, easy plac, dimord= subj bin diff drug 
dum(:,1:3,1,2) = ddmdat(:, 16:18); % bin1-3, easy drug, 
dum(:,1:3,2,1) = ddmdat(:, 19:21); % bin1-3, hard plac, 
dum(:,1:3,2,2) = ddmdat(:, 22:24); % bin1-3, hard drug, 

for idiff=1:2
    dum(:,:,idiff,:) = (dum(:,:,idiff,:) - mean(dum(:,:,idiff,1), 2)) ./ mean(dum(:,:,idiff,1), 2) * 100;
end
% p1 = permtest(dum(:,1,1), dum(:,1,2));
% p2 = permtest(dum(:,2,1), dum(:,2,2));
drug_colors = {'b', 'r'};
for idiff=1:2
        subplot(1,4,idiff+2); hold on
    for idrug = 1:2 % 1:2
        e = errorbar(1:nbins, squeeze(mean(dum(:,:,idiff,idrug))), squeeze(std(dum(:,:,idiff,idrug))) / sqrt(nsub), 'o', ...
            'MarkerFaceColor', drug_colors{idrug});
        e.Color = drug_colors{idrug};
        e.MarkerEdgeColor = [1 1 1];
        e.MarkerSize = 12;
    end
    title(paraleg{ipara})
    %         legend boxoff
    xlim([0 nbins+1])
    ax=gca;
    ax.XTick = 1:nbins;
    ax.XTickLabel = {'Low' 'Medium' 'High'};
    ax.FontSize = 12;
    xlabel('Prestimulus occipital alpha power')
    %         xlabel('Baseline pupil bin')
    ylabel(sprintf('Parameter estimate\n (%% from placebo)'))

end
% title(sprintf('%s\np1 = %g, p2 = %g', paraleg{3}, p1, p2))

if SAV
    %             outpath = fullfile(PREOUT, sesnames{ises}, 'poolings'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'gain');
    mkdir(outpath)
    outfile = fullfile(outpath, sprintf( 'ddm_binned'));
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




