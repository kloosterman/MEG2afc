function MEG2afc_setup_5D_freqstats_corr()
% Make jobs for each condition to run freqanalysis

restoredefaultpath
if ismac
    matlabpath = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/';
    basepath = '/Users/kloosterman/gridmaster2012/projectdata/'; % on the cluster
    parallel = 'parfor'; % local torque
else
    matlabpath = '/home/mpib/kloosterman/MATLAB/';
    basepath = '/home/mpib/kloosterman/projectdata/MEG2afc/'; % on the cluster
    addpath(fullfile(matlabpath, 'tools/qsub_tardis'))
    parallel = 'torque'; % local torque
end
addpath(genpath(fullfile(matlabpath,    'MEG_HH_analysis')))
addpath(fullfile(matlabpath, 'tools', 'fieldtrip-20161220')) %inc JJ edit ft_artifact_zvalue
ft_defaults
addpath(fullfile(matlabpath, 'tools', 'fieldtrip-20161220/external/spm8/')) %for spm_bwlabel


MEG2afc_load_respavg 


compile = 'no';

memreq = 8;
% timreq = 30; % 7 mins per perm * 500 perms = 3500 min
timreq = 350; % 10 mins per perm * 500 perms = 5000 min
% timreq = 15; % 7 mins per perm * 500 perms = 3500 min

testdv = {'modulation'; 'correlation'}; % test power modulation or power-drift correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = {};
ctr=1;
for idiff = 1:4 
    for idrug = 1:4 
        for itrig = 1:2
            for itest = 2%:2
                cfg{ctr}.idiff = idiff;
                cfg{ctr}.idrug = idrug;
                cfg{ctr}.itrig = itrig;
                cfg{ctr}.testdv = testdv{itest}; % mod or corr
                cfg{ctr}.chlabel = chlabel;
                cfg{ctr}.faxis_all = faxis_all;
                cfg{ctr}.frind = frind;
                cfg{ctr}.taxis = taxis;
                cfg{ctr}.PREIN = PREIN;
                ctr=ctr+1;
            end
        end
    end
end

switch parallel
    case 'local'
        cellfun(@MEG2afc_5D_freqstats_corr, cfg, respavg);
    case 'peer'
        peercellfun(@MEG2afc_5D_freqstats_corr, cfg);
    case {'torque' 'qsublocal'}

        setenv('TORQUEHOME', 'yes')
        mkdir('~/qsub'); cd('~/qsub');
        switch compile
            case 'no'
%                 nnodes = 30; % how many licenses available?
%                 stack = round(length(cfg(:))/nnodes); % only used when not compiling 

                qsubcellfun(@MEG2afc_5D_freqstats_corr, cfg, 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', true, 'backend', parallel, 'options', '-l nodes=1:ppn=1');
                
            case 'yes'
                compiledfun = qsubcompile({@MEG2afc_5D_freqstats_corr, @ft_statfun_depsamplesT, @ft_statfun_correlationT, @ft_statistics_montecarlo, @spm_bwlabel},...
                    'toolbox', {'stats'}); 
                qsubcellfun(compiledfun, cfg, 'memreq', memreq, 'timreq', timreq*60, ...
                    'stack', 1, 'StopOnError', false, 'backend', parallel, 'options', '-l nodes=1:ppn=1', ...
                    'UniformOutput', false );
        end
    case 'parfor'
        parpool(3)
        parfor ibatch = 1:length(cfg(:))
            MEG2afc_5D_freqstats_corr(cfg{ibatch}, squeeze(respavg(:,:, :, :, cfg{ibatch}.idrug, cfg{ibatch}.idiff, cfg{ibatch}.itrig, :)))
        end
    otherwise
        error('Unknown backend, aborting . . .\n')
end










