function plotTOPOblank(freq, SAV, soinames)
% plotTOPOblank plots the sensor poolings listed in soi
% give SOIN as argument to plot all

close all
if nargin < 2
    SAV=0;
end

% load('sensorselectionTOM.mat')
load('sensorselection.mat')

% if nargin == 2
%    soinames = {'occ-globresp'}; % 'occ-globresp'    'occ-globresp-l'    'occ-globresp-r'    'motor'    'motor-l'    'motor-r'    };
   soinames = chans.group;
% end
% load('plottopofreq.mat')
% soiUnilateral_cmb_CABMSI_optimal_occ

% for i=1:length(soinames)
%     sois(i) = find(strcmp(soinames{i},  SOIN)); %find soi indices
% end

PREOUT = '~/plots/';
% SAV = 0;

cfg = [];
cfg.layout = 'CTF275.lay';       %neuromag306all neuromag306mag neuromag306planar
cfg.comment = 'no';
cfg.interactive='yes';
cfg.marker = 'on';
cfg.shading            = 'flat';
cfg.style = 'blank' ;
cfg.colorbar = 'no';
cfg.highlight          = 'on';
% cfg.highlight          = 'labels';
cfg.highlightsymbol    = 'o';
cfg.markersize         = 1;
cfg.highlightfontsize = 12;
% cfg.markercolor        = [1 0 0];
cfg.highlightcolor        = [1 0 0];
cfg.highlightsize      = 5;

figure;
iplot=0;
for isoi= 1:length(soinames) 
    if mod(isoi, 6) == 0
        if SAV
            outpath = fullfile(PREOUT); %e.g. plots/poolings/pressoff/lowfreq
            warning off; mkdir(outpath); warning on
            outfile = sprintf('%sPoolingtopo_%d', outpath, isoi );
            disp(outfile)
            orient landscape
            print('-dpdf', outfile)
%             print('-depsc2', outfile)
        end
        figure;
        iplot=0;
    end
    iplot=iplot+1;
    subplot(2,3, iplot)
    cfg.highlightchannel = chans(isoi).sens;

    ft_topoplotTFR(cfg,freq);
    title(gca, sprintf('%s', chans(isoi).group), 'Fontsize', 12);
end


if SAV
    outpath = fullfile(PREOUT); %e.g. plots/poolings/pressoff/lowfreq
    warning off; mkdir(outpath); warning on
    outfile = sprintf('%sPoolingtopo_%d', outpath, isoi );
    disp(outfile)
    orient landscape
    print('-dpdf', outfile)
%     print('-depsc2', outfile)
end