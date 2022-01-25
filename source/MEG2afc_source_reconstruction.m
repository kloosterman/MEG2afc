
%% incorporate our data.
% 30-70 Hz: ctr freq = 50 Hz, is 5 cycles
% Time intervals: -0.4 to -0.3, 0.1 to 0.2

cd('/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/preproc/NK1/plac_ipsi')
load('NK1_ses1_20141115_plac_ipsi_run1_stim_data.mat')

% % make resp-locked
% cfg=[];
% cfg.offset = -round((data.trialinfo(:,5) / ft_findcfg(data.cfg, 'origfs')) * ft_findcfg(data.cfg, 'resamplefs'));
% data = ft_redefinetrial(cfg,data);
                
cfg           = [];                                           
% cfg.toilim    = [-0.4 0.2];           
cfg.toilim    = [-0.4 0.4];           
cfg.minlength = 'maxperlen'; % this ensures all resulting trials are equal length
data          = ft_redefinetrial(cfg, data);

% Now, we 'cut' out the pre- and post-stimulus time windows:
cfg        = [];                                           
% cfg.toilim = [-0.4 -0.3];                       
cfg.toilim = [-0.4 -0.2];                       
data_pre   = ft_redefinetrial(cfg, data);
     
cfg.toilim = [0.2 0.4];                       
data_post   = ft_redefinetrial(cfg, data);

% to not run every function twice, we can combine the two data structure for now. 
cfg      = [];
data_cmb = ft_appenddata(cfg, data_pre, data_post);
% give a number to each trial: 0 = baseline, 1 = experimental condition
data_cmb.trialinfo = [zeros(length(data_pre.trial), 1); ones(length(data_post.trial), 1)];

%%NK edit: interpolate missing channels
label_all = ft_senslabel( 'ctf275');
match = match_str(label_all, data_cmb.label,1);
cfg=[];
cfg.method = 'template';
cfg.template = 'ctf275_neighb.mat';
neighbours = ft_prepare_neighbours(cfg);

cfg=[];
cfg.missingchannel = label_all(find(match==0));
cfg.neighbours = neighbours;
data_cmb = ft_channelrepair(cfg, data_cmb);

% Calculating the cross spectral density matrix
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
% cfg.output     = 'powandcsd'; %takes longer
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 20;
cfg.foi        = 50;
freq_cmb       = ft_freqanalysis(cfg, data_cmb);

% Now, we can separate the two conditions again:
cfg                = [];
cfg.trials         = freq_cmb.trialinfo == 0;
freq_pre           = ft_selectdata(cfg, freq_cmb);
% remember the number of tapers per trial
freq_pre.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
freq_pre.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);

cfg.trials         = freq_cmb.trialinfo == 1;
freq_post           = ft_selectdata(cfg, freq_cmb);
% remember the number of tapers per trial
freq_post.cumtapcnt = freq_cmb.cumtapcnt(cfg.trials);
freq_post.cumsumcnt = freq_cmb.cumsumcnt(cfg.trials);

% Before computing the leadfields, we need to load again our source- and headmodels if they are not in memory anymore:
load /Users/kloosterman/gridmaster2012/projectdata/MEG2afc/MRI/NKdet/NK1/NK1_hdm.mat;
load /Users/kloosterman/gridmaster2012/projectdata/MEG2afc/MRI/NKdet/NK1/NK1_sourcemodel.mat;

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

% beam pre- and poststim by using the common filter
cfg.grid.filter  = source.avg.filter;
source_pre       = ft_sourceanalysis(cfg, freq_pre);
source_post       = ft_sourceanalysis(cfg, freq_post);

source_diff = source_post;
% source_diff.avg.pow = (source_pre.avg.pow ./ source_post.avg.pow) - 1;
source_diff.avg.pow = (source_post.avg.pow ./ source_pre.avg.pow) - 1;

source_diff.pos = template.sourcemodel.pos;
source_diff.dim = template.sourcemodel.dim;

templatefile = '/Users/kloosterman/gridmaster2012/MATLAB/fieldtrip-20141231/external/spm8/templates/T1.nii';
template_mri = ft_read_mri(templatefile);

cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'avg.pow';
cfg.interpmethod = 'nearest';
source_diff_int  = ft_sourceinterpolate(cfg, source_diff, template_mri);

%% plot
        load('colormap170613.mat'); % colormap(cmap);  %cmap = get(gcf, 'Colormap')

% close all;
cfg               = [];
cfg.method        = 'ortho'; % 'slice',  'ortho', 'surface',
cfg.funparameter  = 'avg.pow';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = [0.0 0.2];
cfg.opacitylim    = [0.0 0.2]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap = cmap(129:256,:);
% cfg.camlight ='no';
ft_sourceplot(cfg,source_diff_int);

