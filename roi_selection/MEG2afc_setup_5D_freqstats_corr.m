% function MEG2afc_setup_5D_freqstats_corr()
% Make jobs for each condition to run freqanalysis

% restoredefaultpath
if ismac
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman'; % on the cluster
%     backend = 'none'; % local torque
    backend = 'localloop'; % local torque localloop
%     backend = 'parfor'; % local torque
else
    basepath = '/home/mpib/kloosterman'; % on the cluster
    addpath(fullfile(basepath, 'MATLAB',  'tools/qsub_tardis'))
    backend = 'torque'; % local torque
end
addpath(genpath(fullfile(basepath, 'MATLAB',   'MEG_HH_analysis')))
rmpath(genpath(fullfile(basepath, 'MATLAB', 'MEG_HH_analysis/.git/')))
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip-20161220')) %inc JJ edit ft_artifact_zvalue
ft_defaults
addpath(fullfile(basepath, 'MATLAB', 'tools', 'fieldtrip-20161220/external/spm8/')) %for spm_bwlabel


%%
% %%%%%
% MEG2afc_load_respavg % respavg not loaded when running on the cluster
% % %%%%%%
% % 
% % %%
% driftrates = MEG2afc_load_driftrates161116;
driftrates = MEG2afc_load_drifts_regression;

%% Control panel
testdv = {'modulation'; 'correlation'}; % test modulation or power-drift correlation
itest = 1;

datatype = { 'powermod' 'latr_wrt_resp' 'latr_wrt_choice'  'latr_wrt_stim' }; % use powermodulation or powerlateralization
%idat = 3; % loop across

% freqrange = 'lowhigh'; % low high or full
freqrange = 'full'; % low high or full

corrtype = 'Pearson';
% corrtype = 'Spearman';

%%%%%%%

compile = 'yes';

% memreq = 6000;
memreq = 12000;
% timreq = 30; % 7 mins per perm * 500 perms = 3500 min
timreq = 60; % 0.8 mins per perm * 500 perms = 500 min
% timreq = 15; % 7 mins per perm * 500 perms = 3500 min

overwrite_input = true; % freq that goes into stats
overwrite_output = true; % output stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MEG2afc_sensorselection

cfg = [];
cfg.channel          = 'all'; %{'MEG'}
% cfg.latency          = 'all';
cfg.frequency        = [5 100] %'all';
% cfg.frequency        = [2 40] %'all';
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 1000
cfg.avgoverfreq      = 'no';
% prepare_neighbours determines what sensors may form clusters
cfg0_neighb.method    = 'template'; %TODO do based on data (chans missing)
if ismac
    cfg0_neighb.template  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/ctf275_neighb.mat';
else
    cfg0_neighb.template  = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20161220/template/neighbours/ctf275_neighb.mat';
end
cfg.neighbours       = ft_prepare_neighbours(cfg0_neighb);
% ft_neighbourplot(cfg0_neighb, freq)

cfgcell = {};       ctr = 0;
for idat = 1:4
    data = [];
    switch datatype{idat}
        case 'powermod'
            sensoi = 1:length(chlabel);
            data.label = chlabel;
        otherwise
            sensoi = LR_subtract_mat(:,2); % right hemisphere sensors
            data.label = chlabel(sensoi); % only select right side sensors
    end

    outdir = fullfile(basepath, 'projectdata/MEG2afc/freq/stats_5D', testdv{itest}, datatype{idat}, freqrange, corrtype);
    
    mkdir(outdir)
    mkdir(fullfile(outdir, 'freq'))
    mkdir(fullfile(outdir, 'stat'))
    cd(outdir)

    for idiff = 3 %1:3 %3%:4
        for idrug = 2 %4%3:4
            for itrig = 1:2
                switch testdv{itest}
                    case 'modulation' % FIX this
                        cfg.statistic        = 'ft_statfun_depsamplesT';
                        design = zeros(2,2*nsub);
                        for i = 1:nsub
                            design(1,i) = i;
                            design(1,nsub+i) = i;
                        end
                        design(2,1:nsub)        = 1;
                        design(2,nsub+1:2*nsub) = 2;
                        cfg.design   = design;
                        cfg.uvar     = 1;
                        cfg.ivar     = 2;
                    case 'correlation'
                        cfg.statistic        = 'ft_statfun_correlationT';
                        cfg.type             = corrtype;
                        cfg.design   = squeeze(driftrates(:, idrug, idiff))';
                        cfg.uvar     = [];
                        cfg.ivar     = 1;
                        data.drifts = cfg.design;
                end
                cfg.latency          = XLIM(itrig,:);
%                 if itrig == 1;                    cfg.latency          = [-0.25 1.5];
%                 else
%                     cfg.latency          = [-1.5 0.25];
%                 end

                cfg.inputfile = fullfile(outdir, 'freq', sprintf('freq_idrug%d_idiff%d_itrig%d_itest%d_idat%d_regr.mat', idrug, idiff, itrig, itest, idat));
                if ~exist(cfg.inputfile, 'file') || (overwrite_input && ismac)

                    data.time = taxis{itrig};
                    data.dimord = 'subj_chan_freq_time';
                    
                    switch freqrange
                        case 'lowhigh'
                            data.powspctrm = squeeze(cat(3,  respavg(:,sensoi, 1:length(frind{1}), 1:length(tind{itrig}), idrug, idiff, itrig, 1, idat), ...
                                                             respavg(:,sensoi, 1:length(frind{2}), 1:length(tind{itrig}), idrug, idiff, itrig, 2, idat)  ));
                            data.freq = faxis_all;
                        case 'full'
                            data.powspctrm = squeeze( respavg(:,sensoi, 1:length(frind{3}), 1:length(tind{itrig}), idrug, idiff, itrig, 3, idat));
                            data.freq = faxis{3};
                        case 'low'
                            data.powspctrm = squeeze( respavg(:,sensoi, 1:length(frind{1}), 1:length(tind{itrig}), idrug, idiff, itrig, 1, idat));
                            data.freq = faxis{1};
                        case 'high'
                            data.powspctrm = squeeze( respavg(:,sensoi, 1:length(frind{2}), 1:length(tind{itrig}), idrug, idiff, itrig, 2, idat));
                            data.freq = faxis{2};
                    end
                    
                    fprintf('SAVING %s\n', cfg.inputfile);
                    save(cfg.inputfile, 'data')
%                     load(cfg.inputfile)
                end
                
                cfg.outputfile = fullfile(outdir, 'stat',sprintf('stat_idrug%d_idiff%d_itrig%d_itest%d_idat%d_regr.mat', idrug, idiff, itrig, itest, idat));
                if ~exist(cfg.outputfile, 'file') || overwrite_output % ismac?
                    ctr=ctr+1;
                    cfgcell{ctr} = cfg; % add to queue
                end
                
            end
        end
    end
end

% clear respavg

switch backend
    case 'local'
        cellfun(@ft_freqstatistics, cfgcell);
    case 'localloop' % if multiple inputfiles
        % cellfun(@(a,b) ft_freqstatistics(cfg,a,b), cfgcell);
        for ibatch = 1:length(cfgcell(:))
            cfg = cfgcell{ibatch};
            load(cfg.inputfile)
            datazero = data; %create zero freq to test against
            datazero.powspctrm = zeros(size(data.powspctrm));
            cfg = rmfield(cfg, 'inputfile'); 
            ft_freqstatistics(cfg, data, datazero)
        end


    case 'peer'
        peercellfun(@ft_freqstatistics, cfgcell);
    case {'torque' 'qsublocal'}
        
        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        switch compile
            case 'no'
                %                 nnodes = 30; % how many licenses available?
                %                 stack = round(length(cfgcell(:))/nnodes); % only used when not compiling
                
                qsubcellfun(@ft_freqstatistics, cfgcell, 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 3, 'StopOnError', true, 'backend', backend, 'options', '-l nodes=1:ppn=1');
                
            case 'yes'
                %                 compiledfun = qsubcompile({@ft_freqstatistics, @ft_statfun_depsamplesT, @ft_statfun_correlationT, @ft_statistics_montecarlo, @spm_bwlabel},...
                %                     'toolbox', {'stats'});
                compiledfun = qsubcompile({@ft_freqstatistics, @ft_statfun_depsamplesT, @ft_statfun_correlationT, @ft_statistics_montecarlo, @spm_bwlabel},...
                    'toolbox', {'stats'}, 'executable', 'run_kloosterman_master_p9950_b2.sh');
                qsubcellfun(compiledfun, cfgcell, 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', false, 'backend', backend, 'options', '-l nodes=1:ppn=1', ... % -q gpu
                    'UniformOutput', false );
        end
    case 'parfor'
%         parpool(3)
        parfor ibatch = 1:length(cfgcell(:))
            ft_freqstatistics(cfgcell{ibatch})
        end
    otherwise
        warning('Unknown backend, aborting . . .\n')
end










