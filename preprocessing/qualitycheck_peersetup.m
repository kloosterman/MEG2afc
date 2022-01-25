cfg1 = {};
overwrite = runcfg.overwrite;

%make cells for each subject, to analyze in parallel
ctr = 1;
for isub = 1:length(runcfg.batchlists)
    eval(runcfg.batchlists{isub}); %load in batchlist file, batch, PREOUT and PREIN come out
    for irun=1:length(batch)
        if isempty(batch(irun).dataset) % probably commented out
            fprintf('Empty batchlist entry! Probably commented out?\n')
            continue
        end
        
        basepath = '/mnt/homes/home022/nkloost1/projectdata/2afc/'; %yesno or 2afc
        PREIN = fullfile(basepath, 'data', PRE);
        dataset = fullfile(PREIN, 'meg', batch(irun).dataset);
        [~, out] = fileparts(dataset);
        outfile = fullfile(qualityoutputdir, [out '.pdf']);
        
        if exist(outfile, 'file') && ~overwrite % if the matfile does not yet exist, then add to the joblist
            fprintf('%s exists! Skipping . . . \n', outfile)
            continue
        end
        
        cfg1{ctr}.dataset = dataset;
        cfg1{ctr}.plotunit = 720; % 12 minutes
        cfg1{ctr}.qualityoutputdir = qualityoutputdir;

        ctr =ctr+1;
    end
    clear batch
end

switch runcfg.parallel
    case 'local'
        cellfun(@meg_qualitycheck, cfg1);
    case 'peer'
        peercellfun(@meg_qualitycheck, cfg1);
    case {'torque' 'qsublocal' } %'local'}

        %         timreq = 12; %in minutes per run
        timreq = 2; %in minutes per run
        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        
        nnodes = 32; % how many licenses available?
        stack = round(length(cfg1)/nnodes); % only used when not compiling
        qsubcellfun(@meg_qualitycheck, cfg1, 'memreq', 1024^3, 'timreq', timreq*60, ...
            'stack', stack, 'StopOnError', true, 'backend', runcfg.parallel);
        
    otherwise
        error('Unknown backend, aborting . . .\n')
end
