function MEG2afc_proc_eye_setup()
if ismac
  %     backend = 'parfor';
  backend = 'local';
%   backend = 'torque';
else
  %     backend = 'slurm';
  backend = 'torque';
  backend = 'local';
end

% computation settings
compile = 'yes';
stack = 1;
memreq = 25000; % in MB now! with 16e9 mem_available
timreq = 2*24*60; % days in minutes
overwrite = 0;


%% set up paths
if ismac
%   basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc';
  basepath =  '/Users/kloosterman/beegfs/projectdata/MEG2afc';
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc';
end
% cd(basepath)
% cd('preproc')

PREIN = fullfile(basepath, 'eye', 'edf');
PREOUT = fullfile(basepath, 'eye', 'preproc');
mkdir(PREOUT)

SUBJ= [1:5, 7:9, 11:21];

% sesdirs = {'drug_contra'	'drug_ipsi'	'plac_contra'	'plac_ipsi'};

cfglist = {};
cfg = []; 
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
for isub = SUBJ
  cfg.isub = isub;
  cfglist = [cfglist cfg];
end

fprintf('Running eye analysyis for %d cfgs\n', length(cfglist))
funhandle = str2func('MEG2afc_proc_eye');

% cfglist = cfglist(1)

switch backend
  case 'local'
    cellfun(funhandle, cfglist(:));
  case 'peer'
    peercellfun(funhandle, cfglist(:));
  case {'torque' 'qsublocal' 'slurm'}
    
    setenv('TORQUEHOME', 'yes')
    mkdir('~/qsub'); cd('~/qsub');
    if strcmp(backend, 'slurm')
      options = '-D. -c2 --gres=gpu:1';
    else
      options =  '-l nodes=1:ppn=2'; % torque %-q testing or gpu
    end
    switch compile
      case 'no'
%         nnodes = 60; % how many licenses available?
%         %                 stack = round(length(cfglist(:))/nnodes); % only used when not compiling
%         qsubcellfun(funhandle, cfglist(:), 'memreq', memreq, 'timreq', timreq*60, ...
%           'stack', stack, 'StopOnError', false, 'backend', backend, 'options', '--gres=gpu:1'); %-q gpu
        qsubcellfun(funhandle, cfglist, 'memreq', memreq, 'timreq', timreq*60, ...
          'stack', stack, 'StopOnError', false, 'backend', backend, 'options', options, ...
          'UniformOutput', false );

        
      case 'yes'
        compiledfun = qsubcompile(funhandle, 'toolbox', {'signal', 'stats'});
%         compiledfun = qsubcompile(funhandle, 'toolbox', ...
%           {'signal', 'stats'}, 'executable', ['run_' 'kloosterman_master_p17668_b1' '.sh']);
        
        qsubcellfun(compiledfun, cfglist(:), 'memreq', memreq, 'timreq', timreq*60, ...
          'stack', stack, 'StopOnError', false, 'backend', backend, 'options', options, ...
          'UniformOutput', false );
    end
    
    
  case 'parfor'
    parfor ibatch = 1:length(cfglist(:))
      MEG2afc_mse_analysis(cfglist{ibatch})
    end
  otherwise
    error('Unknown backend, aborting . . .\n')
end
