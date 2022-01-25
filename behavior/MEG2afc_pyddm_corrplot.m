t = readtable('/Users/kloosterman/Downloads/atx_2.csv')
close all
corrdat = table2array(t(:,2:end-2))
figure
imagesc(corr(corrdat), [-1 1])
colorbar

c=[];
c(:,1) = mean([t.v1_0 t.v1_1],2) - mean([t.v0_0 t.v0_1],2); % atx-plac v
c(:,2) = t.a1 - t.a0;
c(:,3) = t.t1 - t.t0;
c(:,4) = t.b1 - t.b0; % drift bias
c(:,5) = t.z1 - t.z0; % starting point
c(:,6) = t.u; % collapsing bound

ct = array2table(c, 'VariableNames', {'v', 'a', 't', 'b', 'z', 'u'})

% figure
% imagesc(corrplot(c), [-1 1])
% colorbar
%%
% corrplot(ct, 'type', 'Spearman', 'testR', 'on')
corrplot(ct, 'type', 'Pearson', 'testR', 'on')
% saveas(gcf, '/Users/kloosterman/Dropbox/tardis_code/MATLAB/MEG_HH_analysis/behavior/ddm_corrplot.pdf')
saveas(gcf, '/Users/kloosterman/Dropbox/tardis_code/MATLAB/MEG_HH_analysis/behavior/ddm_corrplot.png')

%% corr pyddm with quantile fit drift rate
ddm_acc_perrun = squeeze(megdat.behavior.ddm_acc_perrun.v(:,9,4,3,3))
corr(ddm_acc_perrun, ct.v)
figure; scatter(ddm_acc_perrun, ct.v)

dprime = squeeze(megdat.behavior.dprime(:,9,4,3,3))
corr(dprime, ct.v)
figure; scatter(dprime, ct.v)
corr(dprime, ddm_acc_perrun)
figure; scatter(dprime, ddm_acc_perrun)

