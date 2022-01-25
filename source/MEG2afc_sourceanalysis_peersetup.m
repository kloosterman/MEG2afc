% MEG2afc_sourceanalysis_peersetup

cfg = [];
if ismac
    cfg.basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc';
else
    cfg.basepath = '/home/mpib/kloosterman/projectdata/MEG2afc';
end
cfg.trigger = runcfg.trigger;
cfg.trial_duration = [-0.4 2];
% cortical stimulus response, gamma
cfg.bsl_interval = [-0.2 0];  % always stimlocked
cfg.exp_interval = [-0.2 0];
cfg.tapsmofrq = 20;
cfg.foi = 50;
% fill in below:
cfg.hdmfile = []; 
cfg.sourcemfile  = [];
cfg.inputfile = [];
cfg.outputfile = [];
cfg.trials = [];

sourcecfg = {}; ctr = 0;
for i = 1:length(runcfg.batchlists)
    batch=[];   
    eval(runcfg.batchlists{i}); %load in batchlist file, batch, PRE come out
    
    for irun = 1:length(batch) %for each run per subject
        if isempty(batch(irun).dataset); continue; end % probably commented out

        PREIN = fullfile(cfg.basepath, 'preproc', batch(irun).subj, batch(irun).type, filesep);
        PREOUT = fullfile(cfg.basepath, 'source', runcfg.trigger{1}, batch(irun).subj, batch(irun).type, filesep);
        inputfile = sprintf('%s%s_ses%d_%s_%s_run%d_stim_data', PREIN, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
            batch(irun).type, batch(irun).exp);
        
        if ~exist([inputfile '_preprocinfo.mat'])
            disp('Preprocinfo not found. Problem?')
            continue
        end
        load([inputfile '_preprocinfo.mat']) % for debugging
        cfg.trials = 'all';
        
        outputfile = sprintf('%s%s_ses%d_%s_%s_run%d_%s_source_%dto%d_vs_%dto%d_foi%d_smo%d.mat', PREOUT, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
            batch(irun).type, batch(irun).exp, runcfg.trigger{1}, cfg.exp_interval*1000, cfg.bsl_interval*1000, cfg.foi, cfg.tapsmofrq);
        
        if (~exist(outputfile, 'file') && ~isempty(batch(irun).dataset)) || runcfg.overwrite % add to the joblist if outf does not exist and not commented out in batch
            ctr = ctr + 1;
            sourcecfg{ctr}=cfg;
            sourcecfg{ctr}.trials = cfg.trials;
            sourcecfg{ctr}.inputfile = [inputfile '.mat'];
            sourcecfg{ctr}.outputfile = outputfile;
            sourcecfg{ctr}.hdmfile = fullfile(cfg.basepath, 'MRI', 'NKdet', batch(irun).subj, [batch(irun).subj '_hdm.mat']);
            sourcecfg{ctr}.sourcemfile =  fullfile(cfg.basepath, 'MRI', 'NKdet', batch(irun).subj, [batch(irun).subj '_sourcemodel.mat']);
        else
            [dummy, name]=fileparts(outputfile);
            fprintf('%s exists! Skip it\n', name)
        end
        
    end
end
%Submit jobs
fprintf('Running source analysis for %d cfgs\n', length(sourcecfg))
disp(sourcecfg{1})
% switch runcfg.parallel
%     case 'local'
%         cellfun(@MEG2afc_sourceanalysis, sourcecfg);
%     case {'torque'}
        setenv('TORQUEHOME', 'yes')  %    yes or ''
        mkdir('~/qsub'); cd('~/qsub');
        
        memreq = 1; % in GB  memreq*1024^3  3 is added anyway
%         switch runcfg.compile
%             case 'no'
%                 nnodes = 24; % how many licenses?
%                 stack = round(length(sourcecfg)/nnodes);
%                 qsubcellfun(@MEG2afc_sourceanalysis, sourcecfg, 'compile', 'no', ...
%                     'memreq', memreq, 'timreq', runcfg.timreq*60, 'stack', stack, 'StopOnError', false, 'backend', runcfg.parallel, 'options', '-l nodes=1:ppn=1');
%             case 'yes'
%                 compiledfun = qsubcompile(@MEG2afc_sourceanalysis, 'toolbox', {'signal', 'stats'});
%                 qsubcellfun(compiledfun,sourcecfg, ...
%                     'memreq', memreq, 'timreq', runcfg.timreq*60,'stack', 1, 'StopOnError', false, 'backend', runcfg.parallel, 'options', '-l nodes=1:ppn=1');
%         end
        
        if strcmp(runcfg.compile, 'no')
            nnodes = 13; % how many licenses?
            stack = round(length(sourcecfg)/nnodes);
        else
            stack = 1;
        end
        qsubcellfun(@MEG2afc_sourceanalysis, sourcecfg, 'compile', runcfg.compile, ...
            'toolbox', {'signal', 'stats'}, 'memreq', memreq, 'timreq', runcfg.timreq*60, 'stack', stack, 'StopOnError', false, 'backend', runcfg.parallel, 'options', '-l nodes=1:ppn=1');
% end

