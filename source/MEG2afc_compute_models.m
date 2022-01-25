%% MEG2afc_compute_models
% Take mri file and compute headmodel and sourcemodel for each subject and
% plot and save result

% SUBJ  = { 
% %     'NK1'  'NK3'   'NK4'   'NK5'   'NK7'   'NK8'   'NK9'  'NK11'  'NK12'   'NK13'   
%     'NK15' ...
%          'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
%     }; %   'NK14' has no MRI,  'NK6' out,  'NK10' don't exist, 'NK2' incomplete  'NK12' 'NK15' have bad ddm fits

%   SUBJ  = { 
% 'NK2'    }; %   'NK14' has no MRI,  'NK6' out,  'NK10' don't exist, 'NK2' incomplete  'NK12' 'NK15' have bad ddm fits

  SUBJ  = { 
'NK14'    }; %   'NK14' has no MRI,  'NK6' out,  'NK10' don't exist, 'NK2' incomplete  'NK12' 'NK15' have bad ddm fits
%   SUBJ  = { 
% 'NK2'    }; %   'NK14' has no MRI,  'NK6' out,  'NK10' don't exist, 'NK2' incomplete  'NK12' 'NK15' have bad ddm fits

% basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc'; /Users/kloosterman/beegfs/projectdata/MEG2afc/MRI
basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc';

for isub=1:length(SUBJ)
    
    cd(fullfile(basepath, 'MRI', SUBJ{isub}));
    mri_file = dir('*_V2.mri'); % V1 does not work?   
    
    mri = ft_read_mri(mri_file(end).name);
    cfg=[];
    segmentedmri = ft_volumesegment(cfg, mri);
    
    % add anatomical information to the segmentation
    segmentedmri.transform = mri.transform;
    segmentedmri.anatomy   = mri.anatomy;
    
    cfg              = [];
    cfg.funparameter = 'gray';
    ft_sourceplot(cfg,segmentedmri);
    title(SUBJ{isub})

    outfile = [SUBJ{isub} '_hdm_segmentedmri'];
    export_fig(outfile, '-png', '-transparent',  '-depsc') %'-png',

    cfg        = [];
    cfg.method = 'singleshell';
    hdm        = ft_prepare_headmodel(cfg, segmentedmri);
    
    save([SUBJ{isub} '_hdm'], 'hdm')
    
    %%
%     templatedir  = '/Users/kloosterman/gridmaster2012/MATLAB/fieldtrip-20151005/template/sourcemodel';
    templatedir  = '/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20170611/template/sourcemodel';
    template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm'));
    % inverse-warp the template grid to subject specific coordinates
    cfg                = [];
    cfg.grid.warpmni   = 'yes';
    cfg.grid.template  = template.sourcemodel;
    cfg.grid.nonlinear = 'yes'; % use non-linear normalization
    cfg.mri            = mri;
    sourcemodel        = ft_prepare_sourcemodel(cfg);
    
    save([SUBJ{isub} '_sourcemodel'], 'sourcemodel')
    %%
    close all
    hdm_cm = ft_convert_units(hdm, 'cm');
    
    figure; hold on     % plot all objects in one figure
    
    ft_plot_vol(hdm_cm, 'edgecolor', 'none')
    
    alpha 0.4           % make the surface transparent
    
    ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
    
    example_data = dir(fullfile(basepath, 'preproc', SUBJ{isub}, 'drug_contra', '*data.mat'));
    load(fullfile(basepath, 'preproc', SUBJ{isub}, 'drug_contra', example_data(1).name))
    
    ft_plot_sens(data.grad);
    title(SUBJ{isub})
    
    outfile = [SUBJ{isub} '_sourcemodel_mesh'];
    export_fig(outfile, '-png', '-transparent',  '-depsc') %'-png',

end
