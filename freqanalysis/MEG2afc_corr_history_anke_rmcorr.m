function MEG2afc_corr_history_anke_rmcorr(megdat)

rp = megdat.rp_lat_runs;
dimord = megdat.rp_lat_runs_dimord;

disp 'avg over diff'
rp(:,:,:,:,3,:) = mean(rp,5);

disp 'drug-plac'
rp(:,:,4,:,:,:) = rp(:,:,1,:,:,:) - rp(:,:,2,:,:,:);
size(rp)

disp 'select diff drug'
idiff = 3;
idrug = 4;
rpsel = squeeze(rp(:,:,idrug,:,idiff,:)); % subj_runs_ses_repprobvslat

% disp 'plot lat and rp over time'
% figure; hold on
% for ises=1:2
%     temp = squeeze(rpsel(:,:,ises,1));
%     plot(1:8,nanmean(temp))
% end

disp 'collapse ses'
rpsel = reshape(rpsel, 18, 16, 2); % subj_runs_repprobvslat

disp 'demean rp and lat within subj'
rpsel = rpsel - nanmean(rpsel,2);

disp 'plot regression lines single subjects'
figure; hold on
for isub = 1:18
  rpsub = squeeze(rpsel(isub,:,:));
  rpsub = rpsub(~isnan(rpsub(:,1)),:);
  subjfit(isub,:) = polyfit(rpsub(:,1), rpsub(:,2), 1 );
  subjfit(isub,2) = 0;
  temp = polyval(subjfit(isub,:), [-1 1] )
  plot([-1 1], temp, 'Linewidth', 1, 'Color', [0.5 0.5 0.5])
end
temp = polyval(mean(subjfit), [-1 1] )
plot([-1 1], temp, 'Linewidth', 4, 'Color', 'k')
axis square; box on
[h,p]=ttest(subjfit(:,1))
title(p)
  
disp 'lump runs together'
rpsel = reshape(rpsel, [], 2); % subj_runs_repprobvslat
%TODO why nan for lat and not for repprob?

disp 'remove nans'
rpsel = rpsel(~isnan(rpsel(:,1)),:);

disp 'corr and scatter'
[r, p] = corr(rpsel(:,1), rpsel(:,2))
figure; scatter(rpsel(:,1), rpsel(:,2))
title(r)


