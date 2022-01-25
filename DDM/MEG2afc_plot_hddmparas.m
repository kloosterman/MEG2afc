% read in para's from HDDM csv and plot
% MEG2afc_plot_hddmparas

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5bins_pupil.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_10bins_pupil.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_10bins_occalpha.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3binspupil_newbinning.csv';

csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3binsalpha_newbinning.csv';% % works well

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3bins_diff_collapsed2.csv';% 

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_alpha_3bins33.csv'; %5quantiles
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_3bins_diff_collapsed2.csv'; 

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_alpha_3bins402040.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5binsalpha_newbinning.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_1bin.csv';

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5binsalpha_newbins5.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_5binsalpha_newbins5_2.csv';

% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_nobinning_incdiff.csv';
% csvfile = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/params_nobinning_diffcollapsed.csv';

temp = csvread(csvfile, 1,1);
% dims now a drug 1-5, a plac 1-5
nbins = 3;
ddmdat = reshape(temp, 19,nbins,2,3); %dimord subj bin# drug para
ddmdat(:,:,3,:) = mean(ddmdat,3);
nsub = size(ddmdat,1);

% low and high alpha power groups
% ddmdat = ddmdat(lowgroup,:,:,:);
ddmdat = ddmdat(highgroup,:,:,:);

% remove subjects?
% ddmdat = ddmdat(2:end,:,:,:); % drop NK1, extreme drift
% ddmdat = ddmdat([2:7, 9:end],:,:,:); % drop NK1, extreme drift
% ddmdat = ddmdat([1:8, 10:end],:,:,:); % drop NK11, extreme paras
% ddmdat = ddmdat([1:7, 9:end],:,:,:); % drop NK9, extreme drift

% % normalize drug and plac per bin
% for ibin=1:nbins
%     ddmdat2(:,ibin,1:2,:) = (ddmdat(:,ibin,1:2,:) - mean(ddmdat(:,ibin,3,:),2)); ...
% %         ./ mean(ddmdat(:,ibin,3,:),2) ;
% end

% normalize data per subject by dividing by largest drift across bins in placebo (Rajagovindan)
nsub = size(ddmdat,1);
ddmdat2 = []; 
for isub=1:nsub
    for ipara=1:3 % drift etc
        %         ddmdat2(isub,:,:,ipara) = ddmdat(isub,:,:,ipara) ./ max(ddmdat(isub,:,2,ipara)) ;
        ddmdat2(isub,:,:,ipara) = (ddmdat(isub,:,:,ipara) - max(ddmdat(isub,:,1,ipara))) ./ max(ddmdat(isub,:,1,ipara)) * 100 ;
    end
end

%% psc wrt ...
ddmdat_norm = []; range =[];
for idrug = 1:2
    for ipara = 1:3
         ddmdat_norm(:,:,idrug,ipara) = (ddmdat(:,:,idrug,ipara) - mean(ddmdat(:,:,1,ipara),2)) ...
            ./ mean(ddmdat(:,:,1,ipara),2) *100;
        
%         range(:,idrug, ipara) = max(ddmdat_norm(:,:,idrug,ipara), [],2) - min(ddmdat_norm(:,:,idrug,ipara), [],2);
        
    end
end

%% test whether data is best fitted by linear or inverted U (tent) function 
% order = [];  % test if quadratic fit is justified
shapes = [1 2 3; 1 2 1]'; %linear or inverted U
shapes = shapes - mean(shapes); % center around 0
gof = [];
for idrug = 1:2
    for ipara= 3 % drift etc
        %         order(idrug,ipara) = sequential_regression(1:nbins, ddmdat2(:,:,idrug,ipara));
        for iref = 1:2 %1:2 % linear, ushape
            for isub = 1:nsub
                dum = ddmdat_norm(isub,:,idrug,ipara)'; % raw DDM output                
                % regress
                regressors = shapes(:,iref);
                regressors(:,2) = 1;
                betas = regress(dum, regressors );
                % create estimated model response:                
                modelresp = betas(1)*regressors(:,1) + betas(2);
%                 figure; plot([dum modelresp]); legend({'dum', 'modelresp'})
                gof(isub, idrug,ipara,iref) = goodnessOfFit( modelresp, dum, 'NRMSE' ); % measured data is ref
            end
        end
    end
end
squeeze(mean(gof(:,1:2,3,:)))
squeeze(gof(:,2,3,:))
squeeze(mean(gof(:,2,3,:)))

%%
PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/plots';
% plotting
SAV=1;

% addpath( '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/stats')
paraleg = {'Boundary separation', 'Non-decision time', 'Drift rate'};
drug_colors = {'b', 'r'};
% close all
f = figure;
f.Position = [ 1998         466        1069         358];
for ipara = 3
    subplot(1,3,ipara); hold on; box on; axis square;
    for idrug = 1:2
        
%                 dum = ddmdat_norm(:,:,idrug,ipara);
        
        dum = ddmdat2(:,:,idrug,ipara);
        
        % !!!! to make excitability go from low to high
        dum = fliplr(dum)
        
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
        %         xlabel('Prestimulus occipital alpha power')
        xlabel('Neural excitability')
        %         xlabel('Baseline pupil bin')
        %         ylabel(sprintf('Parameter estimate\n (%% wrt grand average)'))
        %         ylabel(sprintf('Parameter estimate\n (%% from placebo max)'))
        ylabel(sprintf('Parameter estimate\n (%% from placebo)'))
        %         ylabel(sprintf('Parameter estimate\n (raw DDM output)'))
        
        fitorder = sequential_regression(1:nbins, dum);
%         fitorder = 1; % linear
%         if ipara == 3
%             fitorder = 2; % quadratic
% %             ylim([-23 -5])
%         end
        
        p = polyfit(1:nbins, mean(dum), fitorder);
        stepsize = 100;
        y = polyval(p, 1:nbins/stepsize:nbins);
        pl = plot(1:nbins/stepsize:nbins, y);
        pl.Color = e.MarkerFaceColor;
        %         pl.DisplayName = respavg.behav_conds{icond};
        pl.LineWidth = 1;
    end
    if ipara==3
        legend({ 'placebo', 'placebo fit', 'drug', 'drug fit'}, 'Location', 'south')
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


%% polyfit 2nd order on single subjects, get quadratic fit (gain), test across subj
% close all
SAV=1;
bin_fits = [];
for isub = 1:nsub
    for ipara = 1:3
        for idrug = 1:2
            %             dum = (ddmdat(isub,:,idrug,ipara) - mean(ddmdat(isub,:,1,ipara),2)) ...
            %                 ./ mean(ddmdat(isub,:,1,ipara),2);
            %             dum = ddmdat_norm(isub,:,idrug,ipara);
            dum = ddmdat2(isub,:,idrug,ipara);
            dum = fliplr(dum);

            bin_fits(isub,idrug,ipara,:) = polyfit(1:nbins, squeeze(dum), order);
            
            %             bin_fits(isub,idrug,ipara,:) = polyfit(1:nbins, squeeze(ddmdat(isub,:,idrug,ipara)), 2);
        end
    end
end

% bin_fits = (bin_fits - mean(bin_fits, 2)) ./ mean(bin_fits, 2) * 100;

% bin_fits = (bin_fits - bin_fits(:,1,:,:)) ./ bin_fits(:,1,:,:) * 100; %psc wrt placebo

% bin_fits = (bin_fits - bin_fits(:,1,:,:)) ./ abs(bin_fits(:,1,:,:)) *
% 100; %psc wrt placebo NOT GOOD B/C values can get negative!

% bin_fits = (bin_fits - bin_fits(:,1,:,:));

% plot
f = figure; iplot=0;
coeff_leg = {'Quadratic' 'Linear' 'Intercept'};
f.Position = [        2324          74         600         250];
for ipara = 3
    for ic = 1:3 % fit coeff
        iplot = iplot+1;
        subplot(1,3,iplot); 
        b = barweb( mean(bin_fits(:,:,ipara,ic)), std(bin_fits(:,:,ipara,ic)) / sqrt(nsub), 0.75 );
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
    %     export_fig(outfile, '-pdf', '-transparent') %'-png',  '-pdf',
    print(outfile, '-dpdf')

%     export_fig(outfile, '-png', '-transparent') %'-png',  '-pdf',
    cd(outpath)
end %ises

% %% partial out the linear from the quadratic term
% ipara = 3; resid=[];
% 
% for idrug=1:2
%     quadrdat = bin_fits(:,idrug,ipara,1);
%     lineardat = bin_fits(:,idrug,ipara,2);
%     
%     [~,~,resid(:,idrug)] = regress(quadrdat, lineardat );
%     figure;
%     scatter(quadrdat, lineardat)
% end
% f=figure;
% barweb(mean(resid) , std(resid)/sqrt(nsub))
% 
% 
% 
% 
% 
% 
%% plot bars of non-binned DDM results
SAV=1;

% addpath( '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/custom_tools/stats')
paraleg = {'Boundary separation', 'Non-decision time', 'Drift rate'};
drug_colors = {'b', 'r'};
% close all
f = figure;
f.Position = [ 1998         466        1069         358];
for ipara = 1:3
    subplot(1,3,ipara); 
    
%     dum = squeeze(ddmdat_norm(:,:,1:2,ipara));
    dum = squeeze(ddmdat(:,:,1:2,ipara));
%     dum = squeeze(ddmdat2(:,:,1:2,ipara));

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

