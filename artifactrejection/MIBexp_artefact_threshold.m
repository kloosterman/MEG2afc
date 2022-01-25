% MIBexp_artefact_threshold
fprintf('\n\nLooking for CAR artifacts . . .\n')
cfg     = [];
cfg.datafile = cfg1.datafile;
cfg.headerfile = cfg1.datafile;
cfg.headerformat = cfg1.headerformat;
cfg.trl = cfg1.trl;
cfg.continuous = 'yes';

cfg.artfctdef.threshold.channel = {'MEG'};
cfg.artfctdef.threshold.bpfilter = 'no';
% cfg.artfctdef.threshold.range = cfg1.carthr;
% cfg.artfctdef.threshold.range = 0.75e-11; % thomas setting
cfg.artfctdef.threshold.range = 5e-12;

[cfg, artfctdef.threshold.artifact ] = ft_artifact_threshold(cfg);

fprintf('%d car artifacts found\n', length(artfctdef.threshold.artifact))
