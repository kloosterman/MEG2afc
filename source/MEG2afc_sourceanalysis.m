function MEG2afc_sourceanalysis(sourcecfg)

% set up paths and variables
trigger = sourcecfg.trigger;
trial_duration = sourcecfg.trial_duration;
exp_interval = sourcecfg.exp_interval;
bsl_interval = sourcecfg.bsl_interval;
tapsmofrq = sourcecfg.tapsmofrq;
foi = sourcecfg.foi ;
hdmfile = sourcecfg.hdmfile ;
sourcemfile = sourcecfg.sourcemfile ;
basepath = sourcecfg.basepath;
inputfile = sourcecfg.inputfile;
outputfile = sourcecfg.outputfile;

if ismac
    matlabpath = '/Users/kloosterman/gridmaster2012/MATLAB'; % '/mnt/homes/home022/nkloost1';
else
    matlabpath = '/home/mpib/kloosterman/MATLAB'; % '/mnt/homes/home022/nkloost1';
end
% addpath(genpath(fullfile(matlabpath, 'MEG_HH_analysis')))
% addpath(fullfile(matlabpath, 'fieldtrip-20141231'))
% ft_defaults

mkdir(fileparts(outputfile));
fprintf('Loading %s\n', inputfile);
try load(inputfile);
    
    % Interpolate missing channels
    label_all = ft_senslabel( 'ctf275');
    match = match_str(label_all, data.label,1);
    cfg=[];
    cfg.method = 'template';
    cfg.template = 'ctf275_neighb.mat';
    neighbours = ft_prepare_neighbours(cfg);
    cfg=[];
    cfg.missingchannel = label_all(find(match==0));
    cfg.neighbours = neighbours;
    cfg.method = 'spline';
    data = ft_channelrepair(cfg, data);

    % Ensure all resulting trials are equal length
    cfg           = [];
    cfg.toilim    = trial_duration ; % eg [-0.4 0.4]
    cfg.minlength = 'maxperlen'; 
    data          = ft_redefinetrial(cfg, data);
    
    % add choice to trialinfo field
    if strfind(inputfile, '_ipsi_')
        data.trialinfo(:,9) = data.trialinfo(:,3);
    elseif strfind(inputfile, '_contra_') % press contra
        data.trialinfo(:,9) = mod(data.trialinfo(:,3),2) + 1;  % resp and choice are flipped
    end
    
    % Now, we 'cut' out the pre- and post-stimulus time windows:
    cfg        = [];
    cfg.toilim = bsl_interval;
    data_bsl   = ft_redefinetrial(cfg, data);

    % make resp-locked if necessary
    if strcmp( trigger, 'resp')
        cfg=[];
        cfg.offset = -round((data.trialinfo(:,5) / ft_findcfg(data.cfg, 'origfs')) * ft_findcfg(data.cfg, 'resamplefs'));
        data = ft_redefinetrial(cfg,data);
    end

    cfg        = [];
    cfg.toilim = exp_interval;
    data_exp   = ft_redefinetrial(cfg, data);
    
    % to not run every function twice, we can combine the two data structure for now.
    cfg      = [];
    data_cmb = ft_appenddata(cfg, data_bsl, data_exp);
    % give a number to each trial: 0 = baseline, 1 = experimental condition
    data_cmb.trialinfo = [zeros(length(data_bsl.trial), 1); ones(length(data_exp.trial), 1)];
    
    % Calculating the cross spectral density matrix
    cfg            = [];
    cfg.method     = 'mtmfft';
    cfg.output     = 'fourier'; % 'powandcsd'; takes longer
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq  = tapsmofrq;
    cfg.foi        = foi;
    freq_cmb       = ft_freqanalysis(cfg, data_cmb);
    
    % Now, we can separate the two conditions again:
    cfg                = [];
    cfg.trials         = freq_cmb.trialinfo == 0;
    freq_bsl           = ft_selectdata(cfg, freq_cmb);
    % remember the number of tapers per trial
    freq_bsl.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
    freq_bsl.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);
    
    cfg.trials         = freq_cmb.trialinfo == 1;
    freq_exp           = ft_selectdata(cfg, freq_cmb);
    % remember the number of tapers per trial
    freq_exp.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
    freq_exp.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);
    
    % Before computing the leadfields, we need to load again our source- and headmodels if they are not in memory anymore:
    load(hdmfile);
    load(sourcemfile)
    
    %% Since we already verified that sensors, head- and sourcemodel align up,
    % we can continue to computing the leadfield matrices by incorporating our just computed frequency data:
    cfg             = [];
    cfg.grid        = sourcemodel;
    cfg.vol         = hdm;
    cfg.channel     = {'MEG'};
    cfg.grad        = freq_cmb.grad;
    sourcemodel_lf  = ft_prepare_leadfield(cfg, freq_cmb);
    
    % Source analysis and contrasting conditions
    cfg                   = [];
    cfg.frequency         = freq_cmb.freq;
    cfg.grad              = freq_cmb.grad;
    cfg.method            = 'dics';
    cfg.keeptrials        = 'yes';
    cfg.grid              = sourcemodel_lf;
    cfg.vol               = hdm;
    cfg.keeptrials        = 'yes';
    cfg.dics.lambda       = '5%';
    cfg.dics.keepfilter   = 'yes';
    cfg.dics.fixedori     = 'yes';
    cfg.dics.realfilter   = 'yes';
    source                = ft_sourceanalysis(cfg, freq_cmb);
    
    templatedir  = fullfile(matlabpath, '/fieldtrip-20151005/template/sourcemodel');
    template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm'));

    % beam pre- and poststim by using the common filter        
    cfg.grid.filter  = source.avg.filter;
    for ichoice = 1:3
        selcfg = [];
        selcfg.trials = find(data.trialinfo(:,9) == ichoice);
        if isempty(selcfg.trials); selcfg.trials = 'all'; end
        
        source_bsl       = ft_sourceanalysis(cfg, ft_selectdata(selcfg, freq_bsl));
        source_exp       = ft_sourceanalysis(cfg, ft_selectdata(selcfg, freq_exp));
        
        source_diff = source_exp;
        source_diff.avg.pow = (source_exp.avg.pow ./ source_bsl.avg.pow) - 1;
        
        source_diff.pos = template.sourcemodel.pos;
        source_diff.dim = template.sourcemodel.dim;
        
        templatefile = fullfile(matlabpath, 'fieldtrip-20141231/external/spm8/templates/T1.nii');
        template_mri = ft_read_mri(templatefile);
        
        cfgint              = [];
        cfgint.voxelcoord   = 'no';
        cfgint.parameter    = 'avg.pow';
        cfgint.interpmethod = 'nearest';
        source_diff_int{ichoice} = ft_sourceinterpolate(cfgint, source_diff, template_mri);
    end
    
    %     %plotting
    %     close all
    %     cfg               = [];
    %     cfg.method        = 'surface';
    %     cfg.funparameter  = 'avg.pow';
    %     cfg.maskparameter = cfg.funparameter;
    %     cfg.funcolorlim   = [0.0 0.2];
    %     cfg.opacitylim    = [0.0 0.2];
    %     cfg.opacitymap    = 'rampup';
    %     ft_sourceplot(cfg,source_diff_int{3});
    
    save(outputfile, 'source_diff_int')

catch ME
    disp(getReport(ME))
    fid = fopen([basepath 'HHmeg2afc_source_errorlog.txt'], 'at');
    fprintf(fid,'%s\n%s\n%s\n\n\n', datestr(now), inputfile, getReport(ME));
    fclose('all');
end

