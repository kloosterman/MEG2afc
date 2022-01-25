function trial_vars = MEG2afc_var_acrosstrials(SUBJ)

% PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/export/13-17and19-25hz_resp';
PREIN = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/export/19-25hz_resp';

% sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
% sesdir_codes = ['B' 'D' 'A' 'C'];

trial_vars = nan(length(SUBJ), 2, 2); % dims subj drug motor

icol = 2; % beta1 or 2

cd(PREIN)
for isub = 1:length(SUBJ)
    seslist = dir([lower(SUBJ{isub}) '_*']);
    for ises=1:length(seslist)
        trialdat = load(seslist(ises).name);
        if isnan(var(trialdat(:,2)))
            disp(length(find(isnan(trialdat(:,2)))))
        end
        if strfind(seslist(ises).name, 'A')
            trial_vars(isub, 2, 1) = var(trialdat(:,icol));
        elseif strfind(seslist(ises).name, 'B')
            trial_vars(isub, 1, 1) = var(trialdat(:,icol));
        elseif strfind(seslist(ises).name, 'C')
            trial_vars(isub, 2, 2) = var(trialdat(:,icol));
        elseif strfind(seslist(ises).name, 'D')
            trial_vars(isub, 1, 2) = var(trialdat(:,icol));
        end
    end
            
end

trial_vars(:,:,3) = nanmean(trial_vars,3);

figure;
% bar(trial_vars(:,:,3))
barweb(mean(trial_vars(:,:,3)), std(trial_vars(:,:,3))/sqrt(length(SUBJ)) )


