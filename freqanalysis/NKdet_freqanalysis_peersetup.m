% MIBexp_freqanalysis_peersetup
% called by runMIBmeg_analysis

overwrite = runcfg.overwrite;
cfg = [];
cfg.output = 'pow';
% cfg.output = 'fourier';
cfg.channel = 'MEG';
cfg.keeptapers  = 'no';
cfg.pad = 7;
cfg.method = 'mtmconvol';
% cfg.precision = 'single';
cfg.runMIBmegcfg = runcfg;

cfg1={}; cfg2={}; inputfile={}; outputfile={}; ctr = 0; %make cells for each subject, to analyze in parallel

for iana=1:length(runcfg.freq.analysistype) %high low
  cfg.freqanalysistype = runcfg.freq.analysistype{iana};
  cfg.phaselocktype = runcfg.freq.phaselocktype{1};
  cfg.timreq = runcfg.freq.timreq;
  cfg.sourceloc = sourceloc;
  cfg.keeptrials  = 'yes'; % needed for fourier-output
  
  switch cfg.freqanalysistype
    case 'high'
      cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
      cfg.foi = 36:2:150;
      cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
      cfg.tapsmofrq = ones(length(cfg.foi),1) .* 8;
    case 'low'
      cfg.taper = 'hanning'; % low frequency-optimized analysis
      cfg.foi = 2:35;
      cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.5;
      cfg.tapsmofrq = ones(length(cfg.foi),1) .* 2;
    case 'full'
      cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
      cfg.foi = 5:2:150;
      cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;   % length of time window = 0.4 sec
      cfg.tapsmofrq = ones(length(cfg.foi),1) .* 5;
    case 'lowfine'
      cfg.taper        = 'hanning';
      cfg.foi          = 2:0.25:35;                      % 0.25 Hz steps
      cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.6;   % length of time window = 0.6 sec
      cfg.pad          = 'nextpow2';
      cfg.tapsmofrq = ones(length(cfg.foi),1) .* 0; % was 1 before
    case 'highfine'
      cfg.taper        = 'hanning';
      cfg.foi = 36:1:150;
      cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.6;   % length of time window = 0.6 sec
      cfg.pad          = 'nextpow2';
      cfg.tapsmofrq = ones(length(cfg.foi),1) .* 4;
    otherwise
      error('Unexpected analysisname. abort.');
  end %switch
  
  for itrg = 1:length(runcfg.trigger) %resp stim
    cfg.trigger=runcfg.trigger{itrg};
    switch cfg.trigger
      case 'resp'
        cfg.toi = -1.5:0.05:0.3;
      case 'stim'
        cfg.toi = -0.5:0.05:1.5; % 2 sec TFR's
      otherwise
        cfg.toi          = -0.3;                           % one timepoint
    end
    
    for i = 1:length(runcfg.batchlists)
      batch=[];
      eval(runcfg.batchlists{i}); %load in batchlist file, batch, PRE come out
      switch cfg.phaselocktype
        case 'evoked'
          %                     MIBexp_concatenatedata(batch, PRE, cfg.trigger, 1)
          batch(length(batch)).exp = length(batch) + 1; % catdata saved under exp+1
          batch = batch(length(batch)); %only analyze catdata
      end
      switch cfg.trigger %concat all trialsmib runs within subj
        case 'trialsmib'
          MIBexp_concatenatedata(batch, PRE, cfg.trigger, 1)
          batch(1).exp = length(batch) + 1; % catdata saved under exp+1
          batch = batch(1); %only analyze catdata
      end
      
      for irun = 1:length(batch) %for each run per subject
        %                 PREIN = ['/mnt/homes/home020/meindertsmat/data/MEG/preproc/' PRE];
        %                 PREOUT = ['/mnt/homes/home020/meindertsmat/data/MEG/freq/' cfg.freqanalysistype filesep PRE cfg.trigger filesep];
        
        %                 basepath = ['/mnt/homes/home022/nkloost1/projectdata/2afc/' ]; %yesno or 2afc
        if ismac
%           basepath = '/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
          basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
        else
%           basepath = '/home/mpib/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
          basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
        end
        %                 PREIN = fullfile(basepath, 'preproc', PRE);
        %                 PREOUT = fullfile(basepath, 'freq', cfg.trigger, PRE);
        %                 infile = sprintf('%s%s_%s%d_%s_data', PREIN, batch(irun).subj, batch(irun).type, batch(irun).exp, 'stim');
        
        PREIN = fullfile(basepath, 'preproc', batch(irun).subj, batch(irun).type, filesep);
        %                 PREOUT = fullfile(basepath, 'freq', cfg.freqanalysistype, cfg.trigger, batch(irun).subj, batch(irun).type, filesep);
        
        PREOUT = fullfile(basepath, 'freq', cfg.freqanalysistype, cfg.trigger, batch(irun).subj, batch(irun).type, filesep);
        infile = sprintf('%s%s_ses%d_%s_%s_run%d_stim_data', PREIN, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
          batch(irun).type, batch(irun).exp);
        
        
        %                 % load in preprocinfo and figure out if condition is there, % if yes put in queue
        %                 preprocinfofile = [infile '_preprocinfo.mat'];
        %                 infile = [infile '.mat' ];
        %                 try load(preprocinfofile)
        %                     switch cfg.trigger
        %                         case 'resp'
        %
        %
        %                     end
        cfg.trials = 'all'; %find(preprocinfo.trl(:,4) == itype & preprocinfo.trl(:,5) == ievent);
        %                     cfg.trials = 1:5; %find(preprocinfo.trl(:,4) == itype & preprocinfo.trl(:,5) == ievent);
        
        %                     outfile = sprintf('%s%s_%s%d_%d_%s_freq.mat', PREOUT, batch(irun).subj, batch(irun).type, ...
        %                         batch(irun).exp, cfg.t_ftimwin(1)*1000,  cfg.phaselocktype );
        outfile = sprintf('%s%s_ses%d_%s_%s_run%d_%s_%s_freq.mat', PREOUT, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
          batch(irun).type, batch(irun).exp, runcfg.trigger{itrg}, cfg.phaselocktype);
        
        if (~exist(outfile, 'file') && ~isempty(batch(irun).dataset)) || overwrite % add to the joblist if outf does not exist and not commented out in batch
          ctr = ctr + 1;
          cfg1{ctr}=cfg;
          cfg2{ctr}.vartrllength = 2; % for var trial length data, nans etc
          inputfile{ctr} = infile;
          outputfile{ctr} = outfile;
          disp(outfile)
        else
          [dummy, name]=fileparts(outfile);
          fprintf('%s exists! Skip it\n', name)
        end
        %                 catch ME
        %                     fprintf('%s not found\n', preprocinfofile)
        %                 end
      end % for irun = 1:length(batch)
    end %length(runcfg.batchlists)
  end %itrg
  
  %Submit jobs
  fprintf('Running MIBexp_freq_analysis for %d cfgs\n', length(cfg1))
  disp(cfg1{1})
  disp(cfg2{1})
  switch runcfg.freq.parallel
    case 'local'
      cellfun(@NKdet_freqanalysis, cfg1, cfg2, inputfile, outputfile);
    case 'peer'
      peercellfun(@NKdet_freqanalysis, cfg1, cfg2, inputfile, outputfile);
    case {'torque' 'slurm'}
      setenv('TORQUEHOME', 'yes')  %    yes or ''
      mkdir('~/qsub'); cd('~/qsub');
      
      memreq = 10000; % in MB
      %             memreq = 5000; % in MB
      if strcmp(runcfg.freq.parallel, 'slurm')
        %       options = '-D. -c2 --gres=gpu:1';
        options = '-D. -c2';
      else
        options =  '-l nodes=1:ppn=2'; % torque %-q testing or gpu
      end

      switch runcfg.freq.compile
        case 'no'
          nnodes = 60; % how many licenses?
%           stack = round(length(cfg1)/nnodes);
          stack = 1;
          qsubcellfun(@NKdet_freqanalysis, cfg1, cfg2, inputfile, outputfile, 'compile', 'no', ...
            'memreq', memreq, 'timreq', cfg.timreq*60, 'stack', stack, 'StopOnError', false, 'backend', runcfg.freq.parallel, 'options', options);
          %                     qsubcellfun(@NKdet_freqanalysis, cfg1, cfg2, inputfile, outputfile, 'compile', 'no', ...
          %                         'memreq', memreq, 'timreq', cfg.timreq*60, 'stack', stack, 'StopOnError', false, 'backend', runcfg.freq.parallel, 'options', '-l nodes=1:ppn=1 -q gpu');
        case 'yes'
          %                     compiledfun = qsubcompile(@NKdet_freqanalysis, 'toolbox', {'signal', 'stats'}); %
          
          compiledfun = qsubcompile(@NKdet_freqanalysis, 'toolbox', {'signal', 'stats'} ...
            , 'executable', ...
            'run_kloosterman_master_p47677_b1.sh');
          
          
          qsubcellfun(compiledfun, cfg1, cfg2, inputfile, outputfile, ...
            'memreq', memreq, 'timreq', cfg.timreq*60,'stack', 1, 'StopOnError', false, 'backend', runcfg.freq.parallel, 'options', options);
      end
    case 'parfor'
      parfor i = 1:length(cfg1)
        tic
        NKdet_freqanalysis(cfg1{i}, cfg2{i}, inputfile{i}, outputfile{i})
        toc
      end
    otherwise
      error('Unknown backend, aborting . . .\n')
  end
  
end % freqanalysistype

