%% Plot topo's: Gamma in visual cortex Stimlocked
close all

SAV= 1;

itrig = 2;

if itrig == 1  
%     bandoi = [50 70];  %gamma stim response: select occ sensors
%     toi = [0.25 0.5]; % stim
%     ZLIM    = [-0.03,0.03];
    bandoi = [4 8];  %gamma stim response: select occ sensors
    toi = [0.2 0.8]; % stim
%     
%     bandoi = [40 80];  %gamma stim response: select occ sensors
%     toi = [0.25 0.5]; % stim
%     ZLIM    = [-0.15,0.15];
%     
%     bandoi = [15 25];  %beta band occipital drug effect
%     toi = [0.25 1]; % stim
%     ZLIM    = [-0.1,0.1];
    
%     bandoi = [49 61];  %gamma band frontal drug effect
%     toi = [0.25 0.75]; % resp
%     ZLIM    = [-0.05,0.05];
    
%     bandoi = [7 14];  %gamma band frontal drug effect
%     toi = [0.25 0.5]; % stim
%     ZLIM    = [-0.02,0.02];

else
%     bandoi = [15 40];  % drug effect in beta
%     toi = [-0.4 0.25]; % resp
%     ZLIM    = [-0.05,0.05];
    
%     bandoi = [49 61];  %gamma band frontal drug effect
%     toi = [-0.4 0.2]; % resp
%     ZLIM    = [-0.025,0.025];
    
    bandoi = [30 70];  %gamma band lateralization effect
    toi = [-0.2 0]; % resp
%     ZLIM    = [-0.025,0.025];
    ZLIM    = [-0.15,0.15];

%     bandoi = [15 40];  % motor beta resp1-resp2
%     toi = [-0.5 0]; % resp
%     ZLIM    = [-0.1,0.1];    
% %     
%     bandoi = [50 80];  %gamma band motor response?
%     toi = [-0.25 0.25]; % resp
%     ZLIM    = [-0.05,0.05];    

end

NKsensorselection

cfg = [];
% cfg.layout = '/home/mpib/kloosterman/MATLAB/toolbox/fieldtrip-20150803/template/CTF275.lay';       % neuromag306cmb neuromag306all neuromag306mag neuromag306planar
cfg.layout = 'CTF275.lay';       %CTF275_RH.lay neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat
% cfg.layout = 'CTF275_helmet.mat';       % neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat
cfg.comment = 'no'; %sprintf('zlim: [%g %g]', ZLIM(1),ZLIM(2)); %no xlim
cfg.marker = 'on';
cfg.shading = 'flat';
cfg.style = 'straight'; %both
% cfg.interpolation =  'v4'; %'linear','cubic','nearest','v4' (default = 'v4') see GRIDDATA
cfg.markersymbol = '.';
cfg.markersize = 1;
freq.dimord = 'chan_freq_time';

for imotor = 3 %1:4 % [1,2,4] %
    
    figure;    iplot=0; hold on
    set(gcf, 'Position', [0 -200 375*3 210*4])

    load('colormap170613.mat'); colormap(cmap);  %cmap = get(gcf, 'Colormap')
    
    for idiff= [1,2,4]
        for ipharm = [1,2,4] % 1:4 % 1:2
            
            ZLIM    = [-0.01,0.01];
            if idiff == 4;  ZLIM = [-0.05 0.05];  end
            if ipharm == 4;  ZLIM = [-0.05 0.05];  end

            for itrig = 2  
                for ifreq = 2
                    for ilatr = 2
                        iplot = iplot+1; subplot(4,3,iplot); hold on
                        %                     iplot = iplot+1; subplot(1,1,iplot); hold on
                        
                        %                     respavg(isub,:,1:length(frind{ifreq}),1:length(tind{itrig}), 1:2, 1:3, itrig, ifreq, 2)
                        if ilatr == 1
                            freq.powspctrm = squeeze(mean(respavg(:,:,:,:, ipharm,idiff, itrig, ifreq, ilatr))); % avg over subj
                            freq.label=chlabel;
                        else
                            freq.powspctrm = squeeze(mean(respavg(:,LR_subtract_mat(:,2),:,:, ipharm,idiff, itrig, ifreq, ilatr))); % avg over subj
                            freq.label=chlabel(LR_subtract_mat(:,2));
                        end

                        freq.freq = faxis{ifreq};
                        freq.time = taxis{itrig};
                        %                     freq.mask = ~isnan(freq.powspctrm);
                        %                     cfg.maskparameter = 'mask';
                        
                        cfg.xlim = toi;
                        cfg.ylim = bandoi;
                        cfg.zlim = ZLIM;
                        
                        cfg.highlight          = 'off';
                        cfg.highlightsymbol    = 'o';
                        cfg.markersize         = 1;
                        %                     cfg.markersize         = 5;
                        cfg.highlightfontsize = 12;
                        % cfg.markercolor        = [1 0 0];
                        %                     cfg.highlightcolor        = [1 0 0];
                        cfg.highlightcolor        = [0 0 0];
                        cfg.highlightsize      = 6;
                        %                     cfg.highlightchannel = occind;
                        %                     cfg.highlightchannel = frontalind;
                        %                     cfg.highlightchannel = motorind;
                        %                     cfg.highlightchannel = [];
                        
                        %                     cfg.style = 'straight_imsat';
                        
                        cfg.channel = ft_channelselection(LR_subtract_mat(:,2), chlabel);
                        cfg.interpolatenan = 'yes';
                        test = ft_topoplotTFR(cfg,freq);
                        
                        title({sprintf('%s %s [%g]', diff_conds{idiff}, trigger_leg{itrig}, ZLIM(2)) ; % motor_conds{imotor},
                            sprintf('%g-%g s, %d-%d Hz, %s', toi, bandoi, pharm_conds{ipharm} )}); % 'FontWeight','bold'
                        
                    end
                end
            end
        end
    end
        if SAV
            outpath = fullfile(PREOUT, 'topos'); 
            warning off; mkdir(outpath); warning on
           
            outfile = sprintf( '%s%stopo_%slocked_freq_%g-%gms_%d-%dHz', outpath, filesep, trigger_leg{itrig}, ...
                toi*1000, bandoi  );
            disp(outfile)

            %             orient portrait
            orient landscape
%             print('-dpdf', outfile)     
            print('-dpng', outfile)     
%             export_fig(outfile, '-eps') %'-png', ,  '-depsc'  '-transparent'

        end;
        
end


%% select X sensors showing gamma
SAV=1;
close all
ZLIM = [-0.02 0.02];
nsens = 10; % select n sensors:
idrug=3; imotor =3; idiff=3; istim=3; iresp=4;

% cfg2.latency = [0.25 0.5]; % for gamma after stim onset
% cfg2.frequency = [40 80];
cfg2.latency = toi;
cfg2.frequency = bandoi;
% bandoi = [40 60];  %gamma band lateralization effect
% toi = [-0.3 0.1]; % resp

cfg2.avgovertime= 'yes';
cfg2.avgoverfreq= 'yes';

% freq.powspctrm = squeeze(respavgpool(:,:,:, idrug, imotor, idiff, istim, iresp));

freq.powspctrm = squeeze(mean(respavg(:,:,:,:, idrug,imotor,idiff, istim, iresp))); % avg over subj

                    
freq2 = ft_selectdata(cfg2,freq);

[~,maxinds] = sort(freq2.powspctrm, 'descend');
[~,mininds] = sort(freq2.powspctrm, 'ascend');


cfg.zlim = ZLIM;
cfg.highlight          = 'labels'; %on
cfg.highlightsymbol    = 'o';
cfg.markersize         = 5;
cfg.highlightfontsize = 6;
% cfg.markercolor        = [1 0 0];
cfg.highlightcolor        = [0 0 1];
%                     cfg.highlightcolor        = [0 0 1];
cfg.highlightsize      = 10;
% cfg.highlightchannel = chlabel(maxinds(1:nsens)); %motorind
cfg.highlightchannel = chlabel([mininds(1:nsens) maxinds(1:nsens)]); %motorind
% cfg.highlightchannel = chlabel(mininds(1:nsens) ); %motorind

figure; colormap(cmap)
ft_topoplotTFR(cfg,freq2);

title({sprintf('%s %s %s [%g]', motor_conds{imotor}, diff_conds{idiff}, trigger, ZLIM(2)) ;
    sprintf('%g-%g s, %d-%d Hz, %s %d sensors', toi, bandoi, pharm_conds{ipharm}, nsens )}); % 'FontWeight','bold'

if SAV
    %                 outpath = fullfile(PREOUT, sesnames{ises}, 'topos'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'topos');
    warning off; mkdir(outpath); warning on
    
    outfile = sprintf( '%s%stopoGAMMA_%slocked_%s_freq_%g-%gms_%d-%dHz', outpath, filesep, trigger, ...
        motor_conds{imotor}, toi*1000, bandoi  );
    disp(outfile)
    
    %             orient portrait
    orient landscape
    %             print('-depsc2', outfile)
    print('-dpdf', outfile)
    %                 print('-dpng', outfile)
    %             print('-depsc', '-adobecs', '-painter', outfile) %, '-cmyk'
    
    %                 print('-dpng', outfile)
end;



%% select X sensors showing beta over motor cortex (resp1-resp2)
SAV=1;
close all
ZLIM = [-0.1 0.1];
nsens = 19; % select n sensors:
idrug=3; imotor =3; idiff=3; istim=3; iresp=4;

% cfg2.latency = [0 1];
cfg2.latency = [-0.5 0];
cfg2.frequency = [15 40];
cfg2.avgovertime= 'yes';
cfg2.avgoverfreq= 'yes';
freq.powspctrm = squeeze(respavgpool(:,:,:, idrug, imotor, idiff, istim, iresp));

freq2 = ft_selectdata(cfg2,freq);

[~,maxinds] = sort(freq2.powspctrm, 'descend');
[~,mininds] = sort(freq2.powspctrm, 'ascend');


cfg.zlim = ZLIM;
cfg.highlight          = 'labels'; %on
cfg.highlightsymbol    = 'o';
cfg.markersize         = 5;
cfg.highlightfontsize = 6;
% cfg.markercolor        = [1 0 0];
cfg.highlightcolor        = [0 0 1];
%                     cfg.highlightcolor        = [0 0 1];
cfg.highlightsize      = 10;
% cfg.highlightchannel = chlabel([mininds(1:nsens) maxinds(1:nsens)]); %motorind
% cfg.highlightchannel = chlabel(mininds(1:nsens) ); %motorind
cfg.highlightchannel = chlabel(maxinds(1:nsens) ); %motorind
figure; colormap(cmap)
ft_topoplotTFR(cfg,freq2);

title({sprintf('%s %s %s [%g]', motor_conds{imotor}, diff_conds{idiff}, trigger, ZLIM(2)) ;
    sprintf('%g-%g s, %d-%d Hz, %s %d sensors', toi, bandoi, pharm_conds{ipharm}, nsens )}); % 'FontWeight','bold'

if SAV
    %                 outpath = fullfile(PREOUT, sesnames{ises}, 'topos'); %e.g. plots/poolings/pressoff/lowfreq
    outpath = fullfile(PREOUT, 'topos');
    warning off; mkdir(outpath); warning on
    
    outfile = sprintf( '%s%stopoGAMMA_%slocked_%s_freq_%g-%gms_%d-%dHz', outpath, filesep, trigger, ...
        motor_conds{imotor}, toi*1000, bandoi  );
    disp(outfile)
    
    %             orient portrait
    orient landscape
    %             print('-depsc2', outfile)
    print('-dpdf', outfile)
    %                 print('-dpng', outfile)
    %             print('-depsc', '-adobecs', '-painter', outfile) %, '-cmyk'
    
    %                 print('-dpng', outfile)
end;



