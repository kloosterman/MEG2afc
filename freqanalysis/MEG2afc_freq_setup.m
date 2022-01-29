function MEG2afc_freq_setup()
% run from runMIBmeg_analysis

if ismac
%   basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/'; %yesno or 2afc
  backend = 'local';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
  backend = 'slurm';
%     backend = 'torque';
%         backend = 'local';
  %   compile = 'yes';
  compile = 'no';
end
timreq = 800; %in minutes per run
memreq = 6000; % in MB

% PREIN = fullfile(basepath, 'preproczap-plus');
% PREOUT = fullfile(basepath, 'freqzap-plus');

% PREOUT = fullfile(basepath, 'freqzap_log_BLperrun_'); % logscaled freq axis

linenoise_rem = 'zapline-plus'; % bandstop zapline zapline-plus DFT
PREIN = fullfile(basepath, sprintf('preproc%s', linenoise_rem));
PREOUT = fullfile(basepath, sprintf('freq%s', linenoise_rem));

mkdir(PREOUT)

overwrite = 1;

SUBJ= [1:5, 7:9, 11:21]; % all
% SUBJ= [1:5, 7:9, 11:13, 15:21]; % 14 already run
% SUBJ= 15;

sesdirs = {'A' 'B' 'C' 'D'};
sesnames = {'plac_ipsi', 'drug_ipsi', 'plac_contra', 'drug_contra'};

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
cfglist = {};
for isub = 1:length(SUBJ)
  for ises = 1:4
    cfg.PREIN = PREIN;
    cfg.SUBJ = sprintf('NK%d', SUBJ(isub));
    cfg.PREOUT = PREOUT;
    cfg.sesname = sesnames{ises};
    cfg.outfile = fullfile(PREOUT, sprintf('%s_%s_freq.mat', cfg.SUBJ, sesnames{ises}));
    if ~exist(cfg.outfile, 'file') || overwrite
      cfglist{end+1} = cfg;
    end
  end
end

cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running MEG2afc_freq for %d cfgs\n', length(cfglist))

if strcmp(backend, 'slurm')
%   options = '-D. --cpus-per-task=4'; % '--cpus 4' '-D. -c2'  --gres=gpu:1
  options = '-D. -c4 '; % '-D. -c2'  --gres=gpu:1
else
  options =  '-l nodes=1:ppn=4'; % torque %-q testing or gpu
end

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');
if strcmp(compile, 'yes')
  fun2run = qsubcompile(@MEG2afc_freq_BL_perrun, 'toolbox', {'signal', 'stats'}); %
  %   fun2run = qsubcompile({@MEG2afc_preproc @sortTrials_MEGhh_2afc @interpolate_blinks}, ...
  %       'executable', 'run_kloosterman_master_p10908_b18.sh'); % compiled function
else
  fun2run = @MEG2afc_freq_BL_perrun;
end

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist)
  return
end

qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', 1, ...
  'StopOnError', false, 'backend', backend, 'options', options);
