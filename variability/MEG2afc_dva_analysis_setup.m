function MEG2afc_dva_analysis_setup()% run from runMIBmeg_analysisif ismac  basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc  backend = 'local';  compile = 'no';else  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc  %     backend = 'slurm';  backend = 'torque';%       backend = 'local';%   compile = 'yes';  compile = 'no';endtimreq = 100; %in minutes per runmemreq = 10000; % in MBPREIN = fullfile(basepath, 'preproc');PREOUT = fullfile(basepath, 'dvalock');mkdir(PREOUT)overwrite = 1;SUBJ= [1:5, 7:9, 11:21]; % all% SUBJ= [1:5, 7:9, 11:13, 15:21]; % 14 already run% SUBJ= 14;sesdirs = {'A' 'B' 'C' 'D'};sesnames = {'plac_ipsi', 'drug_ipsi', 'plac_contra', 'drug_contra'};%make cells for each subject, to analyze in parallelcfg = [];cfg.PREIN = PREIN;cfg.PREOUT = PREOUT;cfglist = {};for isub = 1:length(SUBJ)  for ises = 1:4    cfg.PREIN = PREIN;    cfg.SUBJ = sprintf('NK%d', SUBJ(isub));    cfg.PREOUT = PREOUT;    cfg.sesname = sesnames{ises};    %   if ~exist(cfg.outfile, 'file') || overwrite    cfglist = [cfglist cfg];    %   end  endendcfglist = cfglist(randsample(length(cfglist),length(cfglist)));fprintf('Running MEG2afc_dva_analysis for %d cfgs\n', length(cfglist))if strcmp(backend, 'slurm')  options = '-D. -c2'; % --gres=gpu:1else  options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpuendsetenv('TORQUEHOME', 'yes')mkdir('~/qsub'); cd('~/qsub');if strcmp(compile, 'yes')  fun2run = qsubcompile(@MEG2afc_dva_analysis, 'toolbox', {'signal', 'stats'}); %  %   fun2run = qsubcompile({@MEG2afc_preproc @sortTrials_MEGhh_2afc @interpolate_blinks}, ...  %       'executable', 'run_kloosterman_master_p10908_b18.sh'); % compiled functionelse  fun2run = @MEG2afc_dva_analysis;endif strcmp(backend, 'local')  cellfun(fun2run, cfglist, 'uni', false)  returnendqsubcellfun(fun2run, cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', 1, ...  'StopOnError', false, 'backend', backend, 'options', options);