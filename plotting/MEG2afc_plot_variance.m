% MEG2afc_plot_variance
% take the respvar and plot bars for drug and plac

addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
ft_defaults

trigger = 'resp'

MEG2afc_load_respavg

%% plot bars
close all
figure;    iplot=0; hold on
set(gcf, 'Position', [0 -200 375*3 210*4])

TYP = '2afc'; ipup=3; idrug=1:2;
for isoi = 1%:2% 1:length(SOINsel) %8:10 %[4, 7] %1:length(sois)
    for idiff=3%1:4
%         for ipup= [3 1 2 4] %1:4
%             for idrug = 1:4 % [1:2,4] % 1:2
                dum = squeeze(mean(respvar(:, idrug,3,idiff,3,3,ipup)));
                sem = squeeze(std(respvar(:, idrug,3,idiff,3,3,ipup))) / sqrt(nsub);
                barweb(dum, sem)
                legend(pharm_conds(1:2))
%             end
%         end
    end
end

