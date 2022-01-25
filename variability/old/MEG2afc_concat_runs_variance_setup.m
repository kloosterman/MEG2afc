% MEG2afc_concat_runs_torquesetup
% Make jobs for each subject to concatenate runs and collect sessions

close all
clear all
cd ~

addpath(genpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/MEG_HH_analysis'))
addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/')
ft_defaults
addpath('/mnt/homes/home022/nkloost1/Documents/MATLAB/fieldtrip/qsub')

backend = 'local'; %torque local
memreq = 8;
timreq = 30; % in mins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBJ  = { 'NK1'  'NK3'   'NK4'   'NK5'   'NK7'   'NK8'   'NK9'  'NK11'   'NK12'   'NK13'   'NK14' ...
    'NK15'     'NK16'     'NK17'     'NK18'     'NK19'     'NK20'     'NK21' }; %  'NK6' out,  'NK10' don't exist, 'NK2' incomplete
% SUBJ  = { 'NK12'}; %  'NK6' out
% sesdirs = { 'B' 'D' 'A' 'C'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj

trigger = 'stim'; %  'stim' resp
PREIN = '/mnt/homes/home022/nkloost1/projectdata/2afc/preproc/';

TYP = '2afc';

subjdirlist = {}; triggerlist = {};
for isub = 1:length(SUBJ)
    fprintf('Adding subject %s %s data trials . . . \n',SUBJ{isub}, TYP)
    subjdirlist{isub} = fullfile(PREIN, SUBJ{isub});
    triggerlist{isub} = trigger;
end 

%Submit jobs
fprintf('Running MEG2afc_concat_runs_variance for %d subjects\n', length(subjdirlist))
switch backend
    case 'local'
        cellfun(@MEG2afc_concat_runs_variance, subjdirlist, triggerlist, 'UniformOutput', false);
    case {'torque'}
        setenv('TORQUEHOME', 'yes')  %    yes or ''
        mkdir('~/qsub'); cd('~/qsub');
        
        nnodes = length(SUBJ); % how many jobs?
        stack = round(length(subjdirlist)/nnodes);
        qsubcellfun(@MEG2afc_concat_runs_variance, subjdirlist, triggerlist, 'memreq', memreq*1024^3, 'timreq', timreq*60, ...
            'stack', stack, 'StopOnError', false, 'backend', backend, 'options', '-l nodes=1:ppn=1');
    otherwise
        error('Unknown backend, aborting . . .\n')
end

cd ~