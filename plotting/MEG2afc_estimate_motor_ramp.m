% function [ output_args ] = MEG2afc_estimate_motor_ramp( input_args )
% MEG2afc_estimate_motor_ramp plots time course of motor activity for RT
% bins, for high and low pupil, drug vs placebo
%   Detailed explanation goes here
% 
addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
ft_defaults

clear all; close all

trigger = 'resp'

MEG2afc_load_respavg

%% compute lateralization
NKsensorselection
% pooling order:
%     temp(:,1,:,:) = mean(respdat(:,occind,  :,:),2); %     temp(:,2,:,:) = mean(respdat(:,motorind,:,:),2);%     temp(:,3,:,:) = mean(respdat(:,leftoccind,  :,:),2); %     temp(:,4,:,:) = mean(respdat(:,rightoccind,:,:),2);%     temp(:,5,:,:) = mean(respdat(:,leftmotorind,  :,:),2); %     temp(:,6,:,:) = mean(respdat(:,rightmotorind,:,:),2);

respavg_bands=[];
respavg_bands(:,:,1,:, :,:,:,:, :,:) = mean(respavg(:,:,beta_band_ind,:, :,:,:,:, :,:),3);
respavg_bands(:,:,2,:, :,:,:,:, :,:) = mean(respavg(:,:,gamma_band_ind,:, :,:,:,:, :,:),3);
     
% Lateralization 10 dims
contra=[]; ipsi=[]; respavglat=[];
sois = [3 4 ; 5 6]; % LR occ then motor
%occ wrt target loc
contra(:,1,:,:,: ,:,:,:,:) = respavg_bands(:,sois(1,1),:,:,: ,:,:,2,:,:); %contra1: left sensors, right bp
contra(:,2,:,:,: ,:,:,:,:) = respavg_bands(:,sois(1,2),:,:,: ,:,:,1,:,:); % contra2: right sensors, left bp
ipsi(:,1,:,:,: ,:,:,:,:) = respavg_bands(:,sois(1,1),:,:,: ,:,:,1,:,:); %ipsi1
ipsi(:,2,:,:,: ,:,:,:,:) = respavg_bands(:,sois(1,2),:,:,: ,:,:,2,:,:); % ipsi2
respavglat(:,1,:,:,:, :,:,:,:) =  mean(contra,2) - mean(ipsi,2); % dims: subj_soi_bandoi_toi_npharm_nmotor_ndiff_stimorresp_npup_nRT

%motor wrt button press
contra=[]; ipsi=[];
contra(:,1,:,:,: ,:,:,:,:) = respavg_bands(:,sois(2,1),:,:,: ,:,:,:,2,:); %contra1: left sensors, target right 
contra(:,2,:,:,: ,:,:,:,:) = respavg_bands(:,sois(2,2),:,:,: ,:,:,:,1,:); % contra2: right sensors, target left
ipsi(:,1,:,:,: ,:,:,:,:) = respavg_bands(:,sois(2,1),:,:,: ,:,:,:,1,:); %ipsi1
ipsi(:,2,:,:,: ,:,:,:,:) = respavg_bands(:,sois(2,2),:,:,: ,:,:,:,2,:); % ipsi2
respavglat(:,2,:,:,:, :,:,:,:) =  mean(contra,2) - mean(ipsi,2); % dims: subj_soi_bandoi_toi_npharm_nmotor_ndiff_stimorresp_npup_nRT

% respavg_bands(:,:,1,:, :,:,:,:, :,:) = mean(respavg(:,:,beta_band_ind,:, :,:,:,:, :,:),3);
% respavg_bands(:,:,2,:, :,:,:,:, :,:) = mean(respavg(:,:,gamma_band_ind,:, :,:,:,:, :,:),3);
     
% % Lateralization 11 dimensions
% respavg_bands=[];
% respavg_bands(:,:,1,:, :,:,:,:, :,:,:) = mean(respavg(:,:,beta_band_ind,:, :,:,:,:, :,:,:),3);
% respavg_bands(:,:,2,:, :,:,:,:, :,:,:) = mean(respavg(:,:,gamma_band_ind,:, :,:,:,:, :,:,:),3);
% contra=[]; ipsi=[]; respavglat=[];
% sois = [3 4 ; 5 6]; % LR occ then motor
% %occ wrt target loc
% contra(:,1,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(1,1),:,:,: ,:,:,2,:,:,:); %contra1: left sensors, right bp
% contra(:,2,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(1,2),:,:,: ,:,:,1,:,:,:); % contra2: right sensors, left bp
% ipsi(:,1,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(1,1),:,:,: ,:,:,1,:,:,:); %ipsi1
% ipsi(:,2,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(1,2),:,:,: ,:,:,2,:,:,:); % ipsi2
% respavglat(:,1,:,:,:, :,:,:,:,:) =  mean(contra,2) - mean(ipsi,2); % dims: subj_soi_bandoi_toi_npharm_nmotor_ndiff_stimorresp_npup_nRT
% 
% %motor wrt button press
% contra=[]; ipsi=[];
% contra(:,1,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(2,1),:,:,: ,:,:,:,2,:,:); %contra1: left sensors, target right 
% contra(:,2,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(2,2),:,:,: ,:,:,:,1,:,:); % contra2: right sensors, target left
% ipsi(:,1,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(2,1),:,:,: ,:,:,:,1,:,:); %ipsi1
% ipsi(:,2,:,:,: ,:,:,:,:,:) = respavg_bands(:,sois(2,2),:,:,: ,:,:,:,2,:,:); % ipsi2
% respavglat(:,2,:,:,:, :,:,:,:,:) =  mean(contra,2) - mean(ipsi,2); % dims: subj_soi_bandoi_toi_npharm_nmotor_ndiff_stimorresp_npup_nRT

%% plot modulation time courses, Todo add diff
SAV =1;

% close all
imotor=3; idiff=3; istim=3; iresp=3; ipup=3; %irt= 4 ; icor=3;
color_code = {'r-' 'r--' 'b-' 'b--'};
figure; iplot=0;
set(gcf, 'Position', [0 0 375*2 210*3])
for isoi=1:2
    for iband=1:2
        iplot=iplot+1; subplot(2,2,iplot); hold on
        h=[]; iline=0;
        for idrug=1:2
            for idiff=1:2
                %             dum = squeeze(mean(respavg_bands(:,isoi,iband,:,idrug, imotor, idiff, istim, iresp, ipup, irt, icor))); %average over subj
                %             dum = squeeze(mean(respavg_bands(:,isoi,iband,:,idrug, imotor, idiff, istim, iresp, ipup, irt))); %average over subj
                dum = squeeze(mean(respavg_bands(:,isoi,iband,:,idrug, imotor, idiff, istim, iresp, ipup))); %average over subj
                h(idrug) = plot(taxis,dum, color_code{idrug}, 'Linewidth', 2);
                title({sprintf('%s %s %slocked',  sens.leg{isoi}, bands{iband}, trigger) ;
                    sprintf('%s  %s' )}); % 'FontWeight','bold'
                set(gca, 'Tickdir', 'out')
                xlim([TIMLO TIMHI])
                ylabel('Modulation (%)')
                xlabel(['Time from ' trigger ' (s)'])
            end
        end
        plot([0,0],get(gca,'ylim'), 'k');
        plot([TIMLO TIMHI], [0,0], '--k');
        legend(h, pharm_conds)

    end
end
if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%stc_modulation_%slocked', outpath, filesep, trigger);
    disp(outfile)
    export_fig(outfile, '-eps', '-transparent',  '-depsc') %'-png',
end;


%% plot lateralization time courses

SAV =1;

close all
imotor=3; istimorresp=3; ipup=3; %irt= 4 ; % icor=3;
color_code = {'-' '--'};
figure; iplot=0;
set(gcf, 'Position', [0 0 375*2 210*3])
colors= {'--r' '-r'; '--b' '-b'};
for isoi=1:2
    for iband=1:2
        iplot=iplot+1; subplot(2,2,iplot); hold on
        h=[];        iline=0;
        for idrug=1:2
            for idiff=1:2
                iline=iline+1;
                dum = squeeze(mean(respavglat(:,isoi,iband,:,idrug, imotor, idiff, istimorresp, ipup))); %average over subj
                h(iline) = plot(taxis,dum, colors{idrug,idiff}, 'Linewidth', 2);
                title({sprintf('%s %s %slocked',  sens.leg{isoi}, bands{iband}, trigger) ;
                    sprintf('%s  %s' )}); % 'FontWeight','bold'
                set(gca, 'Tickdir', 'out')
                xlim([TIMLO TIMHI])
                if isoi==1
                    ylabel('Lateralization wrt target loc (%)')
                else
                    ylabel('Lateralization wrt button press (%)')
                end
                xlabel(['Time from ' trigger ' (s)'])
                leg_entries{iline} = [pharm_conds{idrug} ' ' diff_conds{idiff}]
            end
        end
        plot([0,0],get(gca,'ylim'), 'k');
        plot([TIMLO TIMHI], [0,0], '--k');
        if iplot==1
            legend(h, leg_entries, 'location', 'northwest')
        end
    end
end
if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%stc_latr_%slocked', outpath, filesep, trigger);
    disp(outfile)
    export_fig(outfile, '-pdf', '-transparent',  '-depsc') %'-png',
end;

%% fit line to each subject, test drug against plac

SAV =0;

close all
imotor=3; idiff=1:2; istimorresp=3; ipup=3; irt= 4 ; % icor=3;
if strcmp(trigger, 'resp')
    tois = [-0.4 0; -0.5 -0.2]; %beta, gamma toi for fitting
else
    tois = [0.2 0.4; 0 0.4]; %beta, gamma toi for fitting
end

allp=[];
for isoi=1:2
    for iband=2 %1:2
        tind = taxis>=tois(iband,1) & taxis<=tois(iband,2);
        for idrug=1:2
            for idiff=1:2
                figure; iplot=0;
                set(gcf, 'Position', [0 0 375*2 210*3])
                for isub = 1:length(SUBJ)
                    %             for isub = [1:7, 9:16] % throw out NK11, 20 and 21
                    dum = squeeze(respavglat(isub,isoi,iband,tind,idrug, imotor, idiff, istimorresp, ipup)); %average over subj
                    p = polyfit(taxis(tind)', dum, 1);
                    %                 allp(isub, isoi, iband, idrug, :) = p;
                    allp(isub, idrug, idiff) = p(1);
                    % plot line and data
                    iplot=iplot+1; subplot(5,4,iplot); hold on
                    regline = polyval(p, taxis(tind));
                    plot(taxis(tind), dum); plot(taxis(tind), regline)
                    title(sprintf('%s %s %s', SUBJ{isub}, pharm_conds{idrug}, diff_conds{idiff}))
                end
            end
            if SAV
                outpath = fullfile(PREOUT, 'poolings');
                warning off; mkdir(outpath); warning on
                outfile = sprintf( '%s%sfits_per_sub_%slocked_%s_%s', outpath, filesep, trigger, pharm_conds{idrug}, bands{iband});
                disp(outfile)
                %                     export_fig(outfile, '-eps', '-transparent',  '-depsc') %'-png',
                export_fig(outfile, '-pdf', '-transparent',  '-depsc') %'-png',
            end;
        end
        dat = allp - repmat(mean(allp,2),1,2);
%         p = randtest1d(dat(:,1), dat(:,2), 0, 1000)
        figure; barwebNK(squeeze(mean(allp)), squeeze(std(dat)/ sqrt(isub)), 0.5, pharm_conds(1:2))
        shading flat
%         title(sprintf('%s %s-locked, p = %g',bands{iband}, trigger,p))
        title(sprintf('%s %s-locked, p = ',bands{iband}, trigger))
        ylabel(sprintf('Lateralization slope %g - %g s %s', tois(iband,1), tois(iband,2), trigger ))
        legend(diff_conds(1:2)); box off
        SAV=1
        if SAV
            outpath = fullfile(PREOUT, 'poolings');
            warning off; mkdir(outpath); warning on
            outfile = sprintf( '%s%sbar_slopes_test_%slocked_%s_%s', outpath, filesep, trigger, pharm_conds{idrug}, bands{iband});
            disp(outfile)
            %                     export_fig(outfile, '-eps', '-transparent',  '-depsc') %'-png',
            export_fig(outfile, '-pdf', '-transparent',  '-depsc') %'-png',
        end;

    end
end
%     if SAV
%         outpath = fullfile(PREOUT, 'poolings');
%         warning off; mkdir(outpath); warning on
%         outfile = sprintf( '%s%stc_latr_%slocked', outpath, filesep, trigger);
%         disp(outfile)
%         export_fig(outfile, '-eps', '-transparent',  '-depsc') %'-png',
%     end;
    
    
    
    
