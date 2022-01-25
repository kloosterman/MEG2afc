% MEG2afc_concat_runs_torquesetup
% Make jobs for each subject to concatenate runs and collect sessions

close all
clear all
cd ~

if ismac
    matlabpath = '/Users/kloosterman/gridmaster2012/MATLAB'; % '/mnt/homes/home022/nkloost1';
else
    matlabpath = '/home/mpib/kloosterman/MATLAB'; % '/mnt/homes/home022/nkloost1';
end
addpath(genpath(fullfile(matlabpath, 'MEG_HH_analysis')))
addpath(fullfile(matlabpath, 'fieldtrip-20141231'))
ft_defaults
addpath(fullfile(matlabpath, 'fieldtrip-20141231', 'qsub'))

trigger = 'resp'; %  'stim' resp

backend = 'torque' %torque local

% memreq = 6;
timreq = 8; % in mins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJ  = {
    'NK1' ...
    'NK3'  ...
    'NK4'   'NK5'   'NK7'   'NK8'   'NK9'  'NK11'   ...
    'NK12'       'NK13'          'NK15'        'NK16'     'NK17'        'NK18'     'NK19'    'NK20'  'NK21' }; %  'NK6' out,  'NK10' don't exist, 'NK2' incomplete
% SUBJ  = { 'NK15'}; %  'NK6' out  'NK14' no MRI
% sesdirs = { 'B' 'D' 'A' 'C'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj

TYP = '2afc';
if ismac
    PREIN = fullfile('/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/source', trigger);
else
    PREIN = fullfile('/home/mpib/kloosterman/projectdata/MEG2afc/source', trigger);
end
% % % PREOUT = fullfile(cfg.basepath, 'source', runcfg.trigger{1}, batch(irun).subj, batch(irun).type, filesep);


subjdirlist = {}; triggerlist = {};
for isub = 1:length(SUBJ)
    fprintf('Adding subject %s %s source trials . . . \n',SUBJ{isub}, TYP)
    subjdirlist{isub} = fullfile(PREIN, SUBJ{isub});
    triggerlist{isub} = trigger;
end

%Submit jobs
fprintf('Running MEG2afc_concat_source for %d subjects\n', length(subjdirlist))
setenv('TORQUEHOME', 'yes')  %    yes or ''
mkdir('~/qsub'); cd('~/qsub');

nnodes = length(SUBJ); % how many jobs?
stack = round(length(subjdirlist)/nnodes);
qsubcellfun(@MEG2afc_concat_source, subjdirlist, triggerlist, 'memreq', 4, 'timreq', timreq*60, ...
    'stack', stack, 'StopOnError', false, 'backend', backend, 'options', '-l nodes=1:ppn=1');  % , 'matlabcmd', 'matlab'
cd('~/qsub');
% % % % % % save(['argout_' datestr(now)], 'argout')
% cellfun(@MEG2afc_concat_runs, subjdirlist, triggerlist, 'UniformOutput', false);

cd ~

