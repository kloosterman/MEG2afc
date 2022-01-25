% MIBexp_artefact_EOG
fprintf('\n\nLooking for EOG artifacts . . .\n')
cfg     = [];
cfg.datafile = cfg1.datafile;
cfg.headerfile = cfg1.datafile;
cfg.headerformat = cfg1.headerformat;
cfg.trl = data.cfg.trl;
% cfg.padding = 1;
cfg.padding = 0;
cfg.continuous = 'yes';
% cutoff and padding
% select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
cfg.artfctdef.zvalue.channel = {'EEG058'};
cfg.artfctdef.zvalue.cutoff      = cfg1.eogverthr;
cfg.artfctdef.zvalue.trlpadding  = 0.5;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0.1;

cfg.artfctdef.zvalue.bpinstabilityfix = 'reduce';

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [1 15];
cfg.artfctdef.zvalue.bpfiltord  = 3;
cfg.artfctdef.zvalue.hilbert    = 'yes';
% feedback
cfg.artfctdef.zvalue.interactive = cfg1.artf_feedback; % 'yes';
[cfg, artfctdef.eog_ver.artifact] = ft_artifact_zvalue(cfg);

fprintf('%d EOG artifacts found\n', length(artfctdef.eog_ver.artifact))


% %%
% cfg     = [];
% % cfg.dataset = cfg1.datafile;
% cfg.headerfile = cfg1.datafile;
% cfg.datafile = cfg1.datafile;
% cfg.trl = data.cfg.trl;
% cfg.continuous = 'yes';
% cfg.artfctdef.eog.bpfilter   = 'yes';
% cfg.artfctdef.eog.bpfilttype = 'but';
% cfg.artfctdef.eog.bpfreq     = [1 15];
% cfg.artfctdef.eog.bpfiltord  = 4;
% cfg.artfctdef.eog.hilbert    = 'yes';
% 
% cfg.artfctdef.eog.channel      = {'EEG058'};
% cfg.artfctdef.eog.cutoff       = 4;
% cfg.artfctdef.eog.trlpadding   = 0.5;
% cfg.artfctdef.eog.fltpadding   = 0.1;
% cfg.artfctdef.eog.artpadding   = 0.1;
% 
% cfg.artfctdef.eog.bpinstabilityfix = 'reduce';
% cfg.artfctdef.eog.interactive = 'yes';
% [cfg, artifact] = ft_artifact_eog(cfg, data);
% 
