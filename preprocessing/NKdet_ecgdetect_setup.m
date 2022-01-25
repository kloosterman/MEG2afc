% MIBexp_preproc_peersetup
% run from runMIBmeg_analysis

cfg1 = {};
cfg2 = {};
cfg3 = {};
outputfile = {};
overwrite = runcfg.overwrite;

%make cells for each subject, to analyze in parallel
ctr = 1;
for isub = 1:length(runcfg.batchlists)
    eval(runcfg.batchlists{isub}); %load in batchlist file, batch, PREOUT and PREIN come out
    for irun=1:length(batch)
        if isempty(batch(irun).dataset) % probably commented out
            continue
        end
%         basepath = '/mnt/homes/home022/nkloost1/projectdata/2afc/'; %yesno or 2afc
        basepathout = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/'; %yesno or 2afc
        basepath = '/Volumes/NielsDet/2AFC/'; %yesno or 2afc
        
        PREIN = fullfile(basepath, 'data', PRE);
        PREOUT = fullfile(basepathout, 'preproc', batch(irun).subj, batch(irun).type, filesep);

        for itrg = 1:length(runcfg.trigger)
            outfile = sprintf('%s%s_ses%d_%s_%s_run%d_%s_data', PREOUT, batch(irun).subj, batch(irun).sessionno, batch(irun).sessiondate, ...
                batch(irun).type, batch(irun).exp, runcfg.trigger{itrg});
            if ~exist([outfile '_heartbeats.mat'], 'file') || overwrite % if the matfile does not yet exist, then add to the joblist
                outputfile{ctr} = outfile;
                cfg1{ctr}.runcfg = runcfg;  %analysis specifics
                cfg1{ctr}.runcfg.PRE = PRE;  %analysis specifics
                cfg1{ctr}.runcfg.batch = batch(irun); %analysis specifics
                cfg1{ctr}.datafile= fullfile(PREIN, 'meg', batch(irun).dataset);
                cfg1{ctr}.etdatafile= fullfile(PREIN, 'eye', batch(irun).eyedat);
                cfg1{ctr}.headerfile = cfg1{ctr}.datafile;
                cfg1{ctr}.headerformat = ft_filetype(cfg1{ctr}.datafile); 
                cfg1{ctr}.channel = {'EEG059'};
                cfg1{ctr}.fsample = 1200;

                cfg1{ctr}.trialfun = ['sortTrials_MEGhh_2afc']; 
                
                cfg1{ctr}.trialdef.trg = runcfg.trigger{itrg}; %baseline, stim or resp
                cfg1{ctr}.dftfilter = 'yes';
%                 cfg1{ctr}.dftfreq = [50 100 150 200];  % line noise removal using discrete fourier transform
%                 cfg1{ctr}.dftfreq = [49:0.1:51, 99:0.1:101,  149:0.1:151, 199:0.1:201];  % line noise removal using discrete fourier transform
                cfg1{ctr}.padding   = 10;
                cfg1{ctr}.dftfreq = [49:0.1:51, 99.5:0.1:100.5 149.5:0.1:150.5 ];  % line noise removal using discrete fourier transform
                cfg1{ctr}.demean = 'yes';
                cfg1{ctr}.continuous   = 'yes';
                
                cfg1{ctr}.hpfilter      = 'yes';
                cfg1{ctr}.hpfreq        = 0.5;
                cfg1{ctr}.hpinstabilityfix = 'reduce'; % filter order from 6 to 5
                
                cfg1{ctr}.trialdef.begtim = -1.0;  % before stim onset
                cfg1{ctr}.trialdef.endtim = 1; % after report/confidence
                
                cfg1{ctr}.artfrej = runcfg.preproc.artf_rejection; %do artifact rejection? yes
                cfg1{ctr}.artf_feedback = runcfg.preproc.artf_feedback; %feedback for inspection automatic artf detection
                cfg1{ctr}.loadartf = runcfg.preproc.loadartf;%load from file?
                %automatic artf detection: cfg specified in resp scripts
                cfg1{ctr}.carthr = 5e-12;
                cfg1{ctr}.musclethr = batch(irun).musclethr;
                cfg1{ctr}.jumpthr = batch(irun).jumpthr;
                cfg1{ctr}.eogverthr = batch(irun).eogverthr;
                cfg1{ctr}.eoghorthr = batch(irun).eoghorthr;
                
% %                 cfg1{ctr}.musclethr = 7.5; %!!! fixed
%                 % visual artifact rejection parameters
%                 cfg2{ctr}.method   = 'summary'; % channel trial summary
%                 cfg2{ctr}.channel = 'MEG';
%                 cfg2{ctr}.alim     = 1e-10;
%                 cfg2{ctr}.megscale = 1;
%                 cfg2{ctr}.eogscale = 5e-8;
% %                 cfg2{ctr}.layout = '/home/niels/matlab/fieldtrip/template/layout/neuromag306all.lay';       %neuromag306cmb neuromag306all neuromag306mag neuromag306planar
%                 
%                 %resampling parameters
%                 cfg3{ctr}.resample = 'yes';
%                 cfg3{ctr}.fsample = 1200;
%                 cfg3{ctr}.resamplefs = 500;
%                 cfg3{ctr}.detrend = 'yes';
                
                ctr = ctr + 1;
            else
                fprintf('%s exists! Skip it\n', outfile)
            end   %%% if ~exist([outfile '.mat'], 'file')
        end
    end
    clear batch
end
fprintf('Running MIBexp_preproc for %d cfgs\n', length(cfg1))

switch runcfg.preproc.parallel
    case 'local'
        cellfun(@NKdet_ecgdetect, cfg1, outputfile);
    case 'peer'
        peercellfun(@NKdet_ecgdetect, cfg1, outputfile);
    case {'torque' 'qsublocal' } %'local'}
%             ntest=1; % for testing with qsub
%             outputfile=outputfile(1:ntest); cfg1=cfg1(1:ntest); cfg2=cfg2(1:ntest); cfg3=cfg3(1:ntest);

        timreq = 100; %in minutes per run
        memreq = 4; % in GB
%         timreq = 2; %in minutes per run
        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        switch runcfg.preproc.compile
            case 'no'
                nnodes = 10; % how many licenses available?
                stack = round(length(cfg1)/nnodes); % only used when not compiling 
%                 qsubcellfun(@MIBexp_preproc, cfg1, outputfile, 'memreq', memreq*1024^3, 'timreq', timreq*60, ...
%                     'stack', stack, 'StopOnError', true, 'backend', runcfg.preproc.parallel);

                qsubcellfun(@NKdet_ecgdetect, cfg1, outputfile, 'memreq', memreq*1024^3, 'timreq', timreq*60, ...
                    'stack', stack, 'StopOnError', true, 'backend', runcfg.preproc.parallel, 'options', '-l nodes=1:ppn=1');
                
            case 'yes'
                compiledfun = qsubcompile({@NKdet_ecgdetect @sortTrials_MIBexperiment}, 'toolbox', {'signal', 'stats'});
                qsubcellfun(compiledfun, cfg1, outputfile, 'memreq', 4*1024^3, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', true, 'backend', runcfg.preproc.parallel, 'options', '-l nodes=1:ppn=1');
        end
        
    case 'parfor'
        parfor icfg = 1:length(cfg1(:))
            NKdet_ecgdetect(cfg1{icfg}, outputfile{icfg})
        end

    otherwise
        error('Unknown backend, aborting . . .\n')
end

