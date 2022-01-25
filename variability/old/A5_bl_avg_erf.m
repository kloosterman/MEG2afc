function [] = A5_bl_avg_erf(subjects)
% split up the files into reflocked, stimlocked, resplocked and fblocked
% can be run in nojvm (or batch)

tic
addpath(genpath('~/code/MEG/NewPipeline'));
addpath(genpath('~/code/Tools'));
addpath('~/Documents/fieldtrip/');
ft_defaults; warning off;

if ~exist('subjects', 'var'),
    % get all those subjects who have a folder
    cd ~/Data/MEG-PL;
    s = dir('P*');
    s = {s(:).name};
    for i = 1:length(s), subjects(i) = str2num(s{i}(2:3)); end
    
    % determine which subject this job will work on
    for sj = unique(subjects),
        % if this subject and session is not complete
        if ~exist([sprintf('~/Data/MEG-PL/P%02d/MEG/Locked', sj) sprintf('/P%02d_A5started.log', sj)], 'file') && ...
                ~exist([sprintf('~/Data/MEG-PL/P%02d/MEG/Locked', sj) sprintf('/P%02d_A5finished.log', sj)], 'file') && ...
                exist([sprintf('~/Data/MEG-PL/P%02d/MEG/Locked', sj) sprintf('/P%02d_A4finished.log', sj)], 'file'),
            break
        end
    end
end

% compute and combine the planars after this operation
% prepare neighbours
cfgneigh             = [];
cfgneigh.method      = 'template';
cfgneigh.layout      = 'CTF275';
neighbours           = ft_prepare_neighbours(cfgneigh);

for sj = unique(subjects),
    clearvars -except sj subjects neighbours
    
    disp(['Analysing subject ' num2str(sj)]);
    subjectdata = subjectspecifics(sj);
    cd(subjectdata.lockdir);
    system(['touch ' sprintf('~/Data/MEG-PL/P%02d/MEG/Locked', sj) sprintf('/P%02d_A5started.log', sj)]);
    
    % ==================================================================
    % get baseline
    % ==================================================================
    
    load(sprintf('P%02d_ref_locked.mat', sj));
    timelock            = rmfield(timelock, 'cfg'); % remove old cfg
   
    % ==================================================================
    % REF
    % ==================================================================
  
    % select timewin
    cfg = [];
    cfg.latency = [-.3 .9];
    timelock = ft_selectdata(cfg, timelock);
    
    % save separate files
    split_Stim_Corr(sj, timelock, 'ref', neighbours);
    
    % ==================================================================
    % STIM
    % ==================================================================
    
    load(sprintf('P%02d_stim_locked.mat', sj));
    timelock            = rmfield(timelock, 'cfg'); % remove old cfg
    % disp('baseline correction...');
    % tic; timelock.trial       = bsxfun(percchange, timelock.trial, refbl); toc;
    
    % select timewin
    cfg = [];
    cfg.latency = [-.2 1];
    timelock = ft_selectdata(cfg, timelock);
    
    % save separate files
    split_Stim_Corr(sj, timelock, 'stim', neighbours);
    
    % ==================================================================
    % RESP
    % ==================================================================
    
    load(sprintf('P%02d_resp_locked.mat', sj));
    timelock            = rmfield(timelock, 'cfg'); % remove old cfg
    %disp('baseline correction...');
    %tic; timelock.trial       = bsxfun(percchange, timelock.trial, refbl); toc;
    
    % select timewin
    cfg = [];
    cfg.latency = [-.5 .5];
    timelock = ft_selectdata(cfg, timelock);
    
    % save separate files
    split_Stim_Corr(sj, timelock, 'resp', neighbours);
    
    % ==================================================================
    % FEEDBACK
    % ==================================================================
    
    load(sprintf('P%02d_fb_locked.mat', sj));
    timelock            = rmfield(timelock, 'cfg'); % remove old cfg
    %disp('baseline correction...');
    %tic; timelock.trial       = bsxfun(percchange, timelock.trial, refbl); toc;
    
    % select timewin
    cfg = [];
    cfg.latency = [-.3 1.5];
    timelock = ft_selectdata(cfg, timelock);

    % save separate files
    split_Stim_Corr(sj, timelock, 'fb', neighbours);
    
    % delete the temporary file
    delete([sprintf('~/Data/MEG-PL/P%02d/MEG/Locked', sj) sprintf('/P%02d_A5started.log', sj)]);
    system(['touch ' sprintf('~/Data/MEG-PL/P%02d/MEG/Locked', sj) sprintf('/P%02d_A5finished.log', sj)]);
    
end
end


function [] = split_Stim_Corr(sj, timelock, lock, neighbours)
% splits data into error and correct + strong vs weak

subjectdata             = subjectspecifics(sj);
trl                     = timelock.trialinfo;

cfg                     = [];
cfg.trials              = find(trl(:,4) == -1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_stimweak_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg                     = [];
cfg.trials              = find(trl(:,4) == -1 & trl(:,7) == -1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_stimweak_respweak_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,4) == -1 & trl(:,7) == 1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_stimweak_respstrong_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,4) == 1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_stimstrong_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,4) == 1 & trl(:,7) == 1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_stimstrong_respstrong_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,4) == 1 & trl(:,7) == -1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_stimstrong_respweak_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,7) == -1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_respweak_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,7) == 1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_respstrong_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,8) == 1);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_correct_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

cfg = [];
cfg.trials              = find(trl(:,8) == 0);
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_incorrect_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);

% take the mean of all trials
cfg = [];
cfg.trials              = find(~isnan(trl(:,4))); % all trials that have a proper response
cfg.outputfile          = [subjectdata.lockdir sprintf('/P%02d_%s_mean_lockbl.mat', sj, lock)];
select_avg_var(cfg, timelock, neighbours);
end

% subfunction that does the actual selecting of data and writing to disk
function [] = select_avg_var(cfg, data, neighbours)
% should contain cfg.trials and cfg.outputfile

% select the trials for this condition
cfg1            = [];
cfg1.trials     = cfg.trials;
data            = ft_selectdata(cfg1, data);

% % SPLIT ANYTHING THAT'S NOT MEG CHANNELS OFF
cfgSel             = [];
cfgSel.channel     = find(~strncmp(data.label, 'M', 1));
data_extrachans    = ft_selectdata(cfgSel, data);
cfgSel.channel     = find(strncmp(data.label, 'M', 1));
data               = ft_selectdata(cfgSel, data);

% compute planar gradiometers
cfg4                = [];
cfg4.feedback       = 'no';
cfg4.method         = 'template';
cfg4.planarmethod   = 'sincos';
cfg4.channel        = 'MEG'; % this will leave in the chans that are not in the neighbour struct
cfg4.neighbours     = neighbours;
data_planar         = ft_megplanar(cfg4, data);

% remove the single trial data
timelock_planar         = ft_timelockanalysis([], data_planar);
timelock_othersens      = ft_timelockanalysis([], data_extrachans);

% combine again
cfg3                    = [];
cfg3.combinemethod      = 'sum';
timelock_combined_avg   = ft_combineplanar(cfg3, timelock_planar);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the same thing, but now using VAR instead of AVG
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timelock_planar.avg     = timelock_planar.var;
cfg3                    = [];
cfg3.combinemethod      = 'sum';
timelock_combined_var   = ft_combineplanar(cfg3, timelock_planar);

% some fields are gone, reappend
timelock_combined_avg.var = timelock_combined_var.avg;
timelock_combined_avg.dof = timelock_planar.dof(1:length(timelock_combined_avg.label), :);

% write to file
cfgAppend               = [];
cfgAppend.outputfile    = cfg.outputfile;
ft_appendtimelock(cfgAppend, timelock_combined_avg, timelock_othersens);

% add var
timelock_combined_var.var = timelock_combined_var.avg;
timelock_combined_var.dof = timelock_planar.dof(1:length(timelock_combined_avg.label), :);
timelock_othersens.avg    = timelock_othersens.var;

cfgAppend               = [];
cfgAppend.outputfile    = regexprep(cfg.outputfile, '_lockbl.mat', '_var.mat');
ft_appendtimelock(cfgAppend, timelock_combined_var, timelock_othersens);

end
