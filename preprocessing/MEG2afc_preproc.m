function MEG2afc_preproc(cfg)% function MEG2afc_preproc(cfg1, cfg2, cfg3, outputfile)% MEG2afc preproc inc eyedo_only_heartbeats = 0;if ismac  edf2asc = '/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/custom_tools/eyelink/mac/edf2asc';else  edf2asc = '/home/mpib/kloosterman/MATLAB/tools/custom_tools/eyelink/linux/edf2asc';endPREIN = cfg.PREIN;PREOUT = cfg.PREOUT;subjno = cfg.subjno;ses = cfg.ses;irun = cfg.irun;linenoise_rem = cfg.linenoise_rem;% sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj% sesdir_codes = ['B' 'D' 'A' 'C'];if strcmp(ses, 'B'),     cond = 'drug_ipsi' ; %ipharm = 1;   imotor=1;elseif strcmp(ses, 'D'), cond = 'drug_contra' ; %ipharm = 1;   imotor=2;elseif strcmp(ses, 'A'), cond = 'plac_ipsi' ; %ipharm = 2;   imotor=1;elseif strcmp(ses, 'C'), cond = 'plac_contra' ; %ipharm = 2;   imotor=2;endSUBJ = sprintf('NK%d', subjno);megdir = fullfile(PREIN, SUBJ, ses, 'meg' );disp(megdir);cd(megdir);meglist = dir('NK*.ds');disp('put meg files in chronological order')datestring = {};for imeg = 1:length(meglist)  hdr = ft_read_header(meglist(imeg).name);  info = hdr.orig.infods;  datestring{imeg} = [info(find(strcmp('_DATASET_COLLECTIONDATETIME', {info.name}))).data]; %_DATASET_LASTMODIFIEDDATETIME  _DATASET_COLLECTIONDATETIMEendmegdates = datetime(datestring', 'TimeZone','local', 'InputFormat','yyyyMMddHHmmss');[~,idx] = sort(megdates);meglist = meglist(idx); % sorted meg filesmegdates = megdates(idx); % sorted meg acq datesdisp('put eye files in chronological order')eyedir = fullfile(PREIN, SUBJ, ses );cd(eyedir);eyelist = dir('nk*.edf'); % filename does NOT guarantee chrono orderif length(eyelist) ~= length(meglist)  error('N meg and eye files does not match!')enddatestring = {};for ieye = 1:length(eyelist)  datestring{ieye} = eyelist(ieye).name(end-22:end-4);endeyedates = datetime(datestring', 'TimeZone','local', 'InputFormat','yyyy-MM-dd''_''HH.mm.ss');[~,idx] = sort(eyedates);eyelist = eyelist(idx); % sorted eye fileseyedates = eyedates(idx); % sorted eye acq datesdisp('time between meg and eye recordings')datediff = megdates - eyedates;disp(datediff)if any(datediff(2:end) > duration(0,5,0)) % > 5 minutes between meg and eye start is suspicious  error('Eye and meg data do not match, apart maybe from run 1')else  disp('Eye and meg data within 5 mins from each other from run 2 onwards')endeyefile = eyelist(irun);megfile = meglist(irun); % dir(sprintf('*%d.ds', irun));disp('Selected files to process:')disp(megfile)disp(eyefile)cd(megdir);%%disp('load raw meg data')cfg=[]; % HP filter% cfg.hpfilter = 'yes';% cfg.hpfreq =  0.5;  %% cfg.hpfiltord = 4;cfg.dataset = megfile.name;cfg.channel = {'EEG057', 'EEG058', 'EEG059', 'MEG'}; %  58 59 keep %cfg1.channel;cfg.continuous = 'yes';% cfg.bsfilter = 'yes'; % try out% cfg.bsfreq = [49 51];data = ft_preprocessing(cfg);% linenoise_rem = 'bandstop';switch linenoise_rem  case 'zapline-plus'    disp 'Run zapline-plus'    cfg=[];    cfg.resample = 'yes';    cfg.resamplefs = 350;    cfg.detrend = 'no';    data = ft_resampledata(cfg, data);    dat = data.trial{1}(1:end-3,:);        %   clean_data_with_zapline(dat, data.fsample, 'nkeep', 100, 'minsigma', 2, 'noisecompdetectsigma', 2.5);    %     [cleanData, resNremoveFinal, resScores, cfg, plothandles] = ...    %       clean_data_with_zapline(dat, data.fsample, 'noisefreqs', 50,'chunklength',30, 'adaptivesigma',1); % 'searchIndividualNoise',0,        [cleanData, zaplineConfig, analyticsResults, plothandles] = clean_data_with_zapline_plus(dat, data.fsample, 'noisefreqs', 50); %    saveas(gcf, fullfile(PREOUT, 'figures',  sprintf('%s_%s_run%d_zapline-plus_50Hz.png', SUBJ, cond, irun)) )%     [cleanData, zaplineConfig, analyticsResults, plothandles] = clean_data_with_zapline_plus(dat, data.fsample, 'noisefreqs', 50, 'nkeep', 90); %%     saveas(gcf, fullfile(PREOUT, 'figures',  sprintf('%s_%s_run%d_zapline-plus_50/3Hz.png', SUBJ, cond, irun)) )    %     [cleanData, resNremoveFinal, resScores, cfg, plothandles] = clean_data_with_zapline(dat, data.fsample, 50, 'sigmaIncrease', 0, 'nkeep', 271, 'initialSigma', 2.5);    data.trial{1}(1:end-3,:) = cleanData;    %     dat = data.trial{1}(1:end-3,:);%     [cleanData, zaplineConfig, analyticsResults, plothandles] = clean_data_with_zapline_plus(dat, data.fsample, 'noisefreqs', 60, 'maxsigma', 5); %%     data.trial{1}(1:end-3,:) = cleanData;    ft_postamble previous   data   % this copies the data.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"    ft_postamble history    data  % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg  case 'zapline'    x = data.trial{1}(1:end-3,:)'; % TODO make sure it only uses MEG chans!    fline = 50/data.fsample;    nremove = 1;    p = [];    plotflag = 1;    y = nt_zapline(x, fline, nremove, p , plotflag);% %     nsmp = length(data.time{1});%     trl = [];%     trl = [1:12000:nsmp; 12000:12000:nsmp+11999]'; % 10 sec segments%     trl = trl(1:end-1,:); % remove last segment ...%     trl(end,2) = nsmp; % ... and add it to the previous one%     trl(:,3) = 0;%     cfg = [];%     cfg.trl = trl;%     data = ft_redefinetrial(cfg, data);%     %     disp 'Zapline noisetools'%     cfg = [];%     cfg.fline = 50/1200;%     cfg.nremove = 10;%     cfg.p_nkeep = 100;    % p.nfft = 256;%     cfg.plotflag = 1;%     data = ft_nt_zapline(cfg, data);  case 'DFT'    nsmp = length(data.time{1});    trl = [];    trl = [1:12000:nsmp; 12000:12000:nsmp+11999]'; % 10 sec segments    trl = trl(1:end-1,:); % remove last segment ...    trl(end,2) = nsmp; % ... and add it to the previous one    trl(:,3) = 0;    cfg = [];    cfg.trl = trl;    data = ft_redefinetrial(cfg, data);        disp 'DFT filter 50 Hz on segments'    cfg=[];    cfg.dftfilter     = 'yes';    cfg.dftfreq = [49.5:0.1:50.5, 99.5:0.1:100.5,  149.5:0.1:150.5] %[50 100 150]; %[49:0.1:51, 99:0.1:101,  149:0.1:151, 199:0.1:201];  % [50 100 150]; line noise removal using discrete fourier transform    cfg.padding = 0;    data = ft_preprocessing(cfg, data);  case 'bandstop'    disp 'bandstop filter 50 Hz'    cfg=[];    cfg.bsfilter = 'yes'; % try out    cfg.bsfilttype    = 'firws'; % per MNE notch filter    cfg.bsfiltdf = 0.25;    bsfreq = transpose(50:50:100); % per MNE notch filter    width = bsfreq/200;        cfg.bsfreq = [bsfreq-width/2 bsfreq+width/2];%     cfg.bpinstabilityfix = 'yes';    data = ft_preprocessing(cfg, data);endusecleanline = 0;if usecleanline  disp 'tryout cleanline, TODO make work for chunks? Best use zapline with either chunking OR cleanline'  addpath(genpath('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/EEG-Clean-Tools/PrepPipeline'))  signal=[];  signal.data = data.trial{1};  signal.srate = data.fsample;  lineNoiseIn = struct('Fs', 1200, 'lineNoiseChannels', 1:length(data.label)-3, 'lineFrequencies', 50:50:550, 'taperBandWidth', 2, 'taperWindowSize', 4, ...    'taperWindowStep', 1, 'fScanBandWidth', 2, 'tau', 100, 'pad', 0, 'fPassBand', [0 600], 'maximumIterations', 10 ); % ,  tic;[signalclean, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn); toc  dataclean = data;  dataclean.trial = {signalclean.data};end   plotit=1;if plotit   disp 'check if zapline did well'%   cfg=[];%   cfg.trl = [data.sampleinfo(1) data.sampleinfo(end) 0];%   data = ft_redefinetrial(cfg, data)  cfgfreq              = [];  cfgfreq.output       = 'pow';  %   cfgfreq.channel      = 'all';  cfgfreq.method       = 'mtmfft';  cfgfreq.taper        = 'hanning';  cfgfreq.keeptrials   = 'no';  cfgfreq.foilim       = [0 min(256/2, 200)];  cfgfreq.pad='nextpow2';  tempfreq = ft_freqanalysis(cfgfreq, data);% %   figure; plot(tempfreq.freq, mean(tempfreq.powspctrm(1:end-3,:)));%   figure; semilogy(tempfreq.freq, mean(tempfreq.powspctrm(1:end-3,:)))  figure; semilogy(tempfreq.freq, mean(tempfreq.powspctrm(1:end-3,:)))%   saveas(gcf, fullfile(PREOUT, 'figures',  sprintf('%s_%s_run%d_zapline-plus_ftfreq.png', SUBJ, cond, irun)) )  clear tempfreq%   tempfreqcl = ft_freqanalysis(cfgfreq, dataclean)%   hold on; semilogy(tempfreq.freq, mean(tempfreqcl.powspctrm(1:end-3,:)))enddisp('high pass filter the data')cfg=[]; % HP filtercfg.hpfilter = 'yes';cfg.hpfreq =  0.5;  %cfg.hpfiltord = 4;cfg.continuous = 'yes';data = ft_preprocessing(cfg, data);% disp('make 1 s segments for dft line filter, remove line noise')% cfg = [];% cfg.length = 1;% cfg.overlap = 0;% data = ft_redefinetrial(cfg, data);% for itrial = 1:numel(data.trial) % apply notch filters%   for ifreq = 50:50:250%     data.trial{itrial} = ft_preproc_dftfilter(data.trial{itrial}, data.fsample, ifreq);%   end% end% cfg = [];% cfg.trl = [data.sampleinfo(1) data.sampleinfo(end) 0];% data = ft_redefinetrial(cfg, data);%%disp 'Detecting heartbeats . . .' % MOVED TO MEG2afc_readbehavior!!!cfg=[];cfg.trl = [(data.fsample*3) length(data.trial{1})-(data.fsample*3) 0];cfg.continuous = 'yes';cfg.artfctdef.ecg.channel = {'EEG059'};if ismac  cfg.artfctdef.ecg.feedback = 'yes';else  cfg.artfctdef.ecg.feedback = 'no';end[cfg_heartbeats, heartbeats] = ft_artifact_ecg(cfg, data);cfg_heartbeats.heartbeats = heartbeats;inter_beat_durs = diff(heartbeats(:,1)) / data.fsample; % in secinter_beat_durs = inter_beat_durs(zscore(inter_beat_durs) < 3); % remove outlier durs, indicates ECG is loose bpm = length(inter_beat_durs) / (sum(inter_beat_durs)/60);cfg_heartbeats.bpm = bpm;if ~exist(fullfile(PREOUT, 'heartbeats'))  mkdir(fullfile(PREOUT, 'heartbeats'))endoutfile = fullfile(PREOUT, 'heartbeats', sprintf('%s_%s_run%d.mat', SUBJ, cond, irun)); % runno appended belowfprintf('Saving %s\n', outfile)save(outfile, 'cfg_heartbeats');if do_only_heartbeats  disp('Only doing heartbeats')  returnenddisp('define trials')cfg = [];cfg.dataset = megfile.name;cfg.fsample = 1200; % data.fsamplecfg.trialfun = 'sortTrials_MEGhh_2afc';cfg.trialdef.trg = 'stim'; %baseline, stim or resp% cfg.trialdef.begtim = -1.0;  % before stim onset% cfg.trialdef.endtim = 1; % after reportcfg.trialdef.begtim = -0.75;  % before stim onsetcfg.trialdef.endtim = 0.5; % after reportcfg.datatype = 'MEG';cfg.irun = irun;cfg.ses = ses;cfg = ft_definetrial(cfg); % define trials and plot trlevent = cfg.event;trl = cfg.trl;trl(:,1:3) = round(trl(:,1:3)/ (1200/data.fsample)); % resample trl if data was already downsampledcfg = [];cfg.trl = trl;data = ft_redefinetrial(cfg, data);% preproc the eye data to add to meg datadisp('preproc the eye data to add to meg data')[~,eyename] = fileparts( eyefile.name);filename_eye = sprintf('%s.asc', eyename);cd(eyedir)if ~exist(filename_eye)  system(sprintf('%s %s', edf2asc, eyefile.name )); %% convert edf to ascenddisp('preprocess eye data')cfg = [];cfg.dataset          = filename_eye;cfg.montage.tra      = eye(4);cfg.montage.labelorg = {'1', '2', '3', '4'};cfg.montage.labelnew = {'EYE_TIMESTAMP', 'EYE_HORIZONTAL', 'EYE_VERTICAL', 'EYE_DIAMETER'};data_eye = ft_preprocessing(cfg);%disp('interpolate blinks')hdr = ft_read_header(filename_eye); %, 'headerformat', 'eyelink_asc');% eyedate = hdr.orig.header{2}; % TODO check again if files matchdata_eye = interpolate_blinks(hdr, data_eye);disp('low pass filter eye data')cfg = [];cfg.lpfilter = 'yes';cfg.lpfreq = 6; % cf de gee 2019 biorxivcfg.lpfiltord = 3;cfg.channel = 'EYE_DIAMETER';data_pupil = ft_preprocessing(cfg, data_eye);data_eye.trial{1}(4,:) = data_pupil.trial{1}; % put filtered data back in with other chansdisp('make trials in eye using meg triggers')%   1) find block start trg in both meg and eye  % 2) subtract offset between meg and eye  % 3) resample trl: / 1.2runstartmsg = hdr.orig.msg( cellfun(@(x) not(isempty(x)), (strfind(hdr.orig.msg, 'trigger 127'))) );runstartmsg = tokenize(runstartmsg{1});runstart.eye = str2double(runstartmsg{2}) - hdr.orig.dat(1,1) + 1; % sample of run start wrt rec start, only 26?trgval1 = [event(find(strcmp('UPPT001',{event.type}))).value];trgsmp1 = [event(find(strcmp('UPPT001',{event.type}))).sample];runstart.meg = trgsmp1 (find(trgval1 == 127) );runstart.meg = round(runstart.meg / 1.2); % resample meg sample to eye sampling rate: 1200 to 1000 Hzrunstart.diff = runstart.meg - runstart.eye; % delay between runstart in the 2 recstrleye = trl;% trleye(:,1:3) = round( trleye(:,1:3) / 1.2); % trl already at 350 Hz! has to become 1000 Hz (1200/data.fsample)trleye(:,1:3) = round( trleye(:,1:3) / (data.fsample/1000)); % trl already at 350 Hz! has to become 1000 Hz (1200/data.fsample)trleye(:,1:2) = trleye(:,1:2) - runstart.diff;cfg=[];cfg.trl = trleye;data_eye = ft_redefinetrial(cfg, data_eye);if ismac  cfg=[];  %   cfg.event = event;  cfg.preproc.demean = 'yes';  cfg.viewmode = 'vertical';  cfg.channel = 'EYE_DIAMETER';  ft_databrowser(cfg, data_eye)endcd('meg')disp 'MEG artifact rejection'disp 'drop trials with range > 5e-12 Tesla'cfg = [];cfg.channel = 'MEG';tempdata = ft_selectdata(cfg, data);maxpertrl = cell2mat(cellfun(@(x) max(abs(x(:))), tempdata.trial, 'uni', false));cfg = [];cfg.trials = find(maxpertrl < 5e-12);data = ft_selectdata(cfg, data);clear tempdatatrl = trl(cfg.trials,:); % remove from trl too to remove from jump and muscle rejectiondisp 'Looking for JUMP artifacts . . .'% cfg     = [];% cfg.dataset = megfile.name;% cfg.trl = trl;% cfg.continuous = 'yes';% if ismac%   cfg.artfctdef.jump.interactive = 'no';% end% [cfg, artifact_jump] = ft_artifact_jump(cfg);cfg     = [];cfg.trl = trl;cfg.continuous = 'no';if ismac  cfg.artfctdef.jump.interactive = 'yes';end[cfg, artifact_jump] = ft_artifact_jump(cfg, data);disp 'Looking for MUSCLE artifacts . . .'% cfg.artfctdef.muscle.channel = {'MEG'}; % if ismac%   cfg.artfctdef.muscle.interactive = 'no';% end% cfg.artfctdef.muscle.cutoff = -1; % threshold op de abs(min(Z)+cutoff), ft_artifact_zvalue edit% [~, artifact_muscle] = ft_artifact_muscle(cfg);if ismac  cfg.artfctdef.muscle.interactive = 'yes';endcfg.artfctdef.muscle.cutoff = -1;[~, artifact_muscle] = ft_artifact_muscle(cfg, data);disp 'reject artifact trials'cfg = [];cfg.artfctdef.jump.artifact = artifact_jump;cfg.artfctdef.muscle.artifact = artifact_muscle;data  = ft_rejectartifact(cfg, data);disp 'reject meg artefact trials from eye data too'cfg=[];cfg.trials = data.trialinfo(:,7);data_eye = ft_selectdata(cfg, data_eye);% disp 'downsample meg'% cfg=[];% cfg.resample = 'yes';% cfg.resamplefs = 400;% cfg.detrend = 'no';% data = ft_resampledata(cfg, data);disp 'downsample eye to meg'cfg=[];cfg.time = data.time;data_eye = ft_resampledata(cfg, data_eye);% % data and eye could now be appended% datacomb = ft_appenddata([], data, data_eye);disp 'append ECG to eye_data struct'cfg = [];% cfg.channel =  'EEG059';cfg.channel =  {'EEG057' 'EEG058' 'EEG059'};ECGEOGdata = ft_selectdata(cfg, data);data_eye = ft_appenddata([], data_eye, ECGEOGdata);clear ECGEOGdataif ismac  cfg=[]  cfg.viewmode = 'vertical';  cfg.demean = 'yes';  cfg.channel =  {'EEG057' 'EEG058' 'EEG059' 'EYE_HORIZONTAL' 'EYE_VERTICAL'};  ft_databrowser(cfg, data_eye)end  disp 'run ICA'cfg = [];cfg.channel = 'MEG';cfg.method = 'runica';% cfg.method = 'fastica';cfg.runica.stop = 0.00000014;%   cfg.numcomponent = 100;%  cfg.trials = 1; % testingcomp = ft_componentanalysis(cfg, data);if ismac  cfg = [];  cfg.viewmode = 'vertical'; % component  %         cfg.channel = 1:50;       % specify the component(s) that should be plotted  ft_databrowser(cfg, comp)  f = gcf;  f.Position = [69 58 774 1045];    cfg = [];  cfg.component = 1:50;       % specify the component(s) that should be plotted  cfg.layout = 'CTF275';  cfg.comment   = 'no';  cfg.marker = 'off';  figure('units','normalized','outerposition', [0.9995 0.0367 1 0.8775] )  ft_topoplotIC(cfg, comp)enddisp('corr all IC''a with ECG per trial and average across trials')cfg = [];cfg.channel =  'EEG059';ECGdata = ft_selectdata(cfg, data);rho = nan(length(ECGdata.trial), length(comp.label));for itrial = 1:length(ECGdata.trial)  rho(itrial,:) = corr(comp.trial{itrial}', ECGdata.trial{itrial}');endrho = mean(rho);% figure; plot(rho)% [~,ecgcomp] = max(rho); % 1 component with max correlating component% [~,ecgcomp] = find(rho > 0.3); % components with correlating > 0.3[~,ecgcomp] = find(abs(zscore(rho)) > 3); % use zscoredisp('corr all IC''a with eyelink blinks and saccades')eye_events = {'EYE_BLINKS', 'EYE_SACCADES'};eyecomp = {};  cfg = [];if ismac; f = figure; endfor iev = 1:2  cfg.channel = eye_events{iev};  blinkdata = ft_selectdata(cfg, data_eye);    rho = nan(length(blinkdata.trial), length(comp.label));  for itrial = 1:length(blinkdata.trial)    rho(itrial,:) = corr(comp.trial{itrial}', blinkdata.trial{itrial}');  end  rho = nanmean(rho);  if ismac; subplot(1,2,iev); plot(rho); end % figure; plot(zscore(rho))  [~,eyecomp{iev}] = find(abs(zscore(rho)) > 5); % use zscoreenddisp('reject ekg, blink and sacc components')cfg = [];cfg.component = unique([ecgcomp [eyecomp{:}]]); % to be removed component(s)data = ft_rejectcomponent(cfg, comp);%% oostenveld 2019 lmcv beamforming variance based artifact rejection% run after blink and heart artifact removal with ICA to reject trials that escapeddisp 'Find trial variance outliers and index them'cfg = [];cfg.channel = 'MEG';tempdata = ft_selectdata(cfg, data);cfg = [];cfg.cov_cut    = [2, 98]; % not used with zscorecutcfg.badtrs     = [];cfg.bad_trials = [];cfg.method = 'zscorecut'; % zscorecut (abs(min)+1 threshold) or maxmin_perct (original)[selecttrials, cfg] = NM_ft_varcut3(tempdata, cfg, ismac); %https://dx.doi.org/10.1101/795799if ismac  cfg3=[];  cfg3.layout = 'CTF275.lay';  cfg3.metric = 'var';  cfg.channel = 'MEG';  ft_rejectvisual(cfg3, tempdata)enddisp 'Reject bad trials'old_trs= size(data.trial,2);cfg2 = [];cfg2.trials = selecttrials;data = ft_selectdata(cfg2, data);data_eye = ft_selectdata(cfg2, data_eye);fprintf('\nRemaining #trials = %d - %d = %d trials .........\nRemoved trials: ',...  old_trs, length(cfg.bad_trials), size(data.trial,2)); disp(cfg.bad_trials)clear tempdata%%outfile = fullfile(PREOUT, sprintf('%s_%s_run%d_zapline.mat', SUBJ, cond, irun)); % runno appended belowfprintf('Saving %s\n', outfile)save(outfile, 'data', 'data_eye');