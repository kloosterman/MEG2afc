function driftrates = MEG2afc_load_driftrates161116(PREOUT)

if nargin == 0
    PREOUT = '/Users/kloosterman/gridmaster2012/kloosterman/plots'; 
end

% Condition labels
%--------------------------------------------------------------------------
pharm_conds = {'atomox' 'placebo' 'drugscomb' 'atomox - placebo'};
motor_conds = {'pressipsi' 'presscontra' 'regimecomb' 'pressipsi - presscontra'};
diff_conds = {'easy' 'hard' 'diffcomb' 'easy - hard'};
pupilconds = {'pupilhigh' 'pupillow' 'pupilcomb' 'pupilhi-lo'};
stim_conds = {'stimleft' 'stimright' 'allstim' 'stimleft-right'};
choice_conds = {'choiceleft' 'choiceright' 'allchoice' 'choiceleft-right'};
rt_conds = {'slow' 'medium' 'fast' 'allrts'};
correct_conds = {'correct' 'error' 'corr+err'};
sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj


driftrates = [];
if ismac
    % temp = csvread('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/MEG_HH_analysis/behavior/ddmparas161116.csv', 1,0);
%     temp = csvread('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/MEG_HH_analysis/behavior/ddmparas291116.csv', 1,0);
    temp = csvread('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/MEG_HH_analysis/behavior/params_no_blink_trials070417.csv', 1,0);
else
%     temp = csvread('/home/mpib/kloosterman/MATLAB/MEG_HH_analysis/behavior/ddmparas291116.csv', 1,0);
    temp = csvread('/home/mpib/kloosterman/MATLAB/MEG_HH_analysis/behavior/params_no_blink_trials070417.csv', 1,0);
end
% driftrates(:,1,:) = temp(1:19,2:3); % dims subj drug diff
% driftrates(:,2,:) = temp(20:end,2:3); % assuming we have N = 19!!!

% % for 291116 data 
% % NK1 is 3 std's away from mean: outlier
% driftrates(:,1,:) = temp(2:19,2:3); % dims subj drug diff
% driftrates(:,2,:) = temp(21:end,2:3); % assuming we have N = 19!!!
% 
% driftrates(:,3,:) = mean(driftrates,2);
% driftrates(:,4,:) = driftrates(:,1,:) - driftrates(:,2,:);
% driftrates(:,:,3) = mean(driftrates,3); % diff
% driftrates(:,:,4) = driftrates(:,:,1) - driftrates(:,:,2);

% 070416 data: NK1 already out
driftrates(:,1,:) = temp(1:18,2:3); % dims subj drug diff
driftrates(:,2,:) = temp(19:end,2:3); % assuming we have N = 19!!!

driftrates(:,3,:) = mean(driftrates,2);
driftrates(:,4,:) = driftrates(:,1,:) - driftrates(:,2,:);
driftrates(:,:,3) = mean(driftrates,3); % diff
driftrates(:,:,4) = driftrates(:,:,1) - driftrates(:,:,2);



% outliercutoff = mean(driftrates(:,3,3)) + 3*std(driftrates(:,3,3)); % = 1.9242 w N=19
% 1.9586 for NK1, so official outlier!

% close all
SAV=1;
% peasy = randtest1d(driftrates(:,1, 1), driftrates(:,2, 1), 0, 1000);
% phard = randtest1d(driftrates(:,1, 2), driftrates(:,2, 2), 0, 1000);
% pall = randtest1d(driftrates(:,1, 3), driftrates(:,2, 3), 0, 1000)

[~,peasy] = ttest(driftrates(:,1, 1), driftrates(:,2, 1));
[~,phard] = ttest(driftrates(:,1, 2), driftrates(:,2, 2));
[~,pall] = ttest(driftrates(:,1, 3), driftrates(:,2, 3))
% 

if ~ismac
    return
end

figure;
barweb(squeeze(mean(driftrates(:,1:2, 1:2)))', squeeze(std(driftrates(:,1:2, 1:2)))'/sqrt(length(driftrates)))
title(sprintf('N = %d, Easy p = %g, Hard p = %g', size(driftrates,1), peasy, phard))
ylim([0 1.4])
set(gca, 'xticklabel', diff_conds(1:2))
legend(pharm_conds{1:2})
if SAV
    outpath = fullfile(PREOUT, 'poolings');
    warning off; mkdir(outpath); warning on
    outfile = sprintf( '%s%s_driftrates', outpath, filesep); % motor_conds{imotor}, pupilconds{ipup}
    disp(outfile)
    export_fig(outfile, '-png', '-painters', '-transparent', '-depsc', '-nocrop') %'-png',  '-pdf',
end

% % corrtype = 'Pearson';
% corrtype = 'Spearman';
