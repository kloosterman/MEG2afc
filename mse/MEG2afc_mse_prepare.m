function MEG2afc_mse_prepare()
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

% analysis settings
removeERP = 'yes';
removeERPmethod = 'subtract';
trigger = {'stim'};

%% set up paths
if ismac
%   basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc';
  basepath =  '/Users/kloosterman/beegfs/projectdata/MEG2afc';
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc';
end
cd(basepath)
cd('preproc')

SUBJ= dir('NK*');

% output
PREOUT = fullfile(basepath, 'mse');
mkdir(PREOUT)

sesdirs = {'drug_contra'	'drug_ipsi'	'plac_contra'	'plac_ipsi'};

cfglist = {};     ctr=0;
for isub = 1:length(SUBJ)
  cd(SUBJ(isub).name)
  
  for ises= 1:4
    try cd(sesdirs{ises})
    catch
      fprintf('! ! ! %s not found\n', sesdirs{ises})
      continue
    end
    runlist = dir('*data.mat');
    for irun = 1:length(runlist)   %1:length(runlist)
      
      %         PREIN = fullfile(basepath, 'preproc', SUBJ{isub}, sprintf('ses%d', ises));
      %         if ~exist(PREIN, 'dir')
      %             fprintf('%s not found\n', PREIN)
      %             continue
      %         end
      cfg = [];
      cfg.fileload = fullfile(pwd, runlist(irun).name);
      cfg.outputfile = fullfile(PREOUT, runlist(irun).name); % all in one folder
      cfg.trigger = trigger;
      cfg.removeERP = removeERP;
      cfg.removeERPmethod = removeERPmethod;
      cfg.SUBJ = SUBJ(isub).name;

      if ~overwrite && exist(cfg.outputfile, 'file')
        fprintf('%s exists!\n', cfg.outputfile)
      else
        fprintf('Adding %s . . . \n', cfg.outputfile)
        ctr = ctr + 1;
        cfglist{ctr} = cfg;
      end

    end
    cd ..
  end
  cd ..
end

fprintf('Running entropyanalysis for %d cfgs\n', length(cfglist))
funhandle = str2func('MEG2afc_mse_analysis');

cfglist = cfglist(1)

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


