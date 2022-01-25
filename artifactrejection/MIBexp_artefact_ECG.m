% MIBexp_artefact_EOG
fprintf('\n\nLooking for ECG artifacts . . .\n')
  


cfg     = [];
cfg.datafile = cfg1.datafile;
cfg.headerfile = cfg1.datafile;
cfg.headerformat = cfg1.headerformat;
cfg.trl = ft_findcfg(data.cfg, 'trl');
cfg.continuous = 'yes';

cfg.artfctdef.ecg.channel = {'EEG059'};
cfg.artfctdef.ecg.pretim  = 0.05; %pre-artifact rejection-interval in seconds
cfg.artfctdef.ecg.psttim  = 0.3;  %post-artifact rejection-interval in seconds
cfg.artfctdef.ecg.method  = 'zvalue'; %peak-detection method
cfg.artfctdef.ecg.cutoff  = 2.5; %peak-threshold
cfg.artfctdef.ecg.inspect = 'EEG059'; %Nx1 list of channels which will be shown in a QRS-locked average

cfg.artfctdef.ecg.hpfreq = 0.5; %Nx1 list of channels which will be shown in a QRS-locked average
cfg.artfctdef.ecg.hpfiltord = 5; %Nx1 list of channels which will be shown in a QRS-locked average
cfg.artfctdef.ecg.hpfilttype  = 'but'; %Nx1 list of channels which will be shown in a QRS-locked average

% [cfg, artifact] = ft_artifact_ecg(cfg);
[cfg, artifact] = ft_artifact_ecg(cfg,data);

export_fig([outputfile 'Heartbeats'], '-png') %'-png',  '-pdf',
% export_fig([outputfile 'Heartbeats'], '-pdf') %'-png',  '-pdf',
close all

% 
% % cutoff and padding
% % select a set of channels on which to run the artifact detection (e.g. can be 'MEG')
% cfg.artfctdef.zvalue.channel = {'EEG059'};
% cfg.artfctdef.zvalue.cutoff      = cfg1.eoghorthr;
% cfg.artfctdef.zvalue.trlpadding  = 0;
% cfg.artfctdef.zvalue.artpadding  = 0.1;
% cfg.artfctdef.zvalue.fltpadding  = 0.2;
% 
% % algorithmic parameters
% cfg.artfctdef.zvalue.bpfilter   = 'yes';
% cfg.artfctdef.zvalue.bpfilttype = 'but';
% cfg.artfctdef.zvalue.bpfreq     = [1 15];
% cfg.artfctdef.zvalue.bpfiltord  = 4;
% cfg.artfctdef.zvalue.hilbert    = 'yes';
% % feedback
% cfg.artfctdef.zvalue.interactive = cfg1.artf_feedback;
% [cfg, artfctdef.eog_hor.artifact] = ft_artifact_zvalue(cfg);
% 
% fprintf('%d EOG artifacts found\n', length(artfctdef.eog_hor.artifact))
