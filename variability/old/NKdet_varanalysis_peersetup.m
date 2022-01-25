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
    cfg.phaselocktype = runcfg.freq.phaselocktype{iana};
    cfg.timreq = runcfg.freq.timreq;
    cfg.sourceloc = sourceloc;    

    switch cfg.freqanalysistype
        case 'high'
            cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
            cfg.keeptrials  = 'no';
            cfg.foi = 36:2:150;
            cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
            cfg.tapsmofrq = ones(length(cfg.foi),1) .* 8;
        case 'low'
            cfg.taper = 'hanning'; % low frequency-optimized analysis
            cfg.keeptrials  = 'yes'; % needed for fourier-output
%             cfg.keeptapers = 'yes'; % idem
            cfg.foi = 3:35;
            cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
            cfg.tapsmofrq = ones(length(cfg.foi),1) .* 4.5;
        case 'full'
            cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
            cfg.keeptrials  = 'yes';
            cfg.foi = 5:2:150;
            cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;   % length of time window = 0.4 sec
            cfg.tapsmofrq = ones(length(cfg.foi),1) .* 5; 
        otherwise
            error('Unexpected analysisname. abort.');
    end %switch
    
    for itrg = 1:length(runcfg.trigger) %resp stim
        cfg.trigger=runcfg.trigger{itrg};
        switch cfg.trigger
            case 'resp'
                cfg.toi = -1.5:0.05:0.3;
            case 'stim'
                 cfg.toi = -0.5:0.05:2.5; % 3 sec TFR's
            otherwise               
        end
        
        for i = 1:length(runcfg.batchlists)
            batch=[];
            eval(runcfg.batchlists{i}); %load in batchlist file, batch, PRE come out
            switch cfg.phaselocktype
                case 'evoked'
                    MIBexp_concatenatedata(batch, PRE, cfg.trigger, 1)
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
                
                basepath = ['/mnt/homes/home022/nkloost1/projectdata/2afc/' ]; %yesno or 2afc
%                 PREIN = fullfile(basepath, 'preproc', PRE);
%                 PREOUT = fullfile(basepath, 'freq', cfg.trigger, PRE);
%                 infile = sprintf('%s%s_%s%d_%s_data', PREIN, batch(irun).subj, batch(irun).type, batch(irun).exp, 'stim');

                PREIN = fullfile(basepath, 'preproc', batch(irun).subj, batch(irun).type, filesep);
                PREOUT = fullfile(basepath, 'freq', cfg.trigger, batch(irun).subj, batch(irun).type, filesep);
                infile = sprintf('%s%s_ses%d_%s_%s_run%d_stim_data', PREIN, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
                    batch(irun).type, batch(irun).exp);

                
                % load in preprocinfo and figure out if condition is there, % if yes put in queue
                preprocinfofile = [infile '_preprocinfo.mat'];
                infile = [infile '.mat' ];
                try load(preprocinfofile)
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
                    else
                        [dummy, name]=fileparts(outfile);
                        fprintf('%s exists! Skip it\n', name)
                    end
                catch ME
                    fprintf('%s not found\n', preprocinfofile)
                end
            end % for irun = 1:length(batch)
        end %length(runcfg.batchlists)
    end %itrg
    
    %Submit jobs
    fprintf('Running MIBexp_freq_analysis for %d cfgs\n', length(cfg1))
    disp(cfg1{1})
    disp(cfg2{1})
    switch runcfg.freq.parallel
        case 'local'
            cellfun(@NKdet_varanalysis, cfg1, cfg2, inputfile, outputfile);
        case {'torque'}
            setenv('TORQUEHOME', 'yes')  %    yes or ''
            mkdir('~/qsub'); cd('~/qsub');

            memreq = 8; % in GB
            switch runcfg.freq.compile
                case 'no'
                    nnodes = 24; % how many licenses?
                    stack = round(length(cfg1)/nnodes);
                    qsubcellfun(@NKdet_varanalysis, cfg1, cfg2, inputfile, outputfile, 'compile', 'no', ...
                        'memreq', memreq*1024^3, 'timreq', cfg.timreq*60, 'stack', stack, 'StopOnError', false, 'backend', runcfg.freq.parallel, 'options', '-l nodes=1:ppn=1');
                case 'yes'
                    compiledfun = qsubcompile(@NKdet_varanalysis, 'toolbox', {'signal', 'stats'});
                    qsubcellfun(compiledfun, cfg1, cfg2, inputfile, outputfile, ...
                        'memreq', memreq*1024^3, 'timreq', cfg.timreq*60,'stack', 1, 'StopOnError', false, 'backend', runcfg.freq.parallel, 'options', '-l nodes=1:ppn=1');
            end
        otherwise
            error('Unknown backend, aborting . . .\n')
    end
    
end % freqanalysistype

