function freq = MEG2afc_freq(cfg)
% check resplocked freq, better do fixed effects with appenddata?, no good for
% removing ERP per run
PREIN = cfg.PREIN;
SUBJ = cfg.SUBJ;
sesname = cfg.sesname;
outfile = cfg.outfile;

cd(PREIN)

freq = {};
runlist = dir(sprintf('%s_%s_run*.mat', SUBJ, sesname)); % put simon together for now
% rep_probs = nan(length(runlist),1);
trls = [];
for irun = 1:length(runlist)
  data = {};
  disp('Loading');    disp(runlist(irun).name)
  load(runlist(irun).name);

  trl = ft_findcfg(data_eye.cfg, 'trl');
  disp 'drop trials w RT <0.3 s'
  cfg=[];
  cfg.trials = data.trialinfo(:,5) > (0.3*1200);
  data = ft_selectdata(cfg, data);
%   disp 'compute behavior per run'
%   disp 'rep prob'
%   %   trl = data.cfg.previous.previous{1}.previous.previous.previous.previous.trl;
%   trl = ft_findcfg(data_eye.cfg, 'trl');
%   button = trl(:,6);  % the button pressed by the subject
%   rep_probs(irun,1) = sum(diff(button) == 0) / (numel(button)-1); % trials with same button as previous trial / ntrials
 
  disp('realigning MEG')
  cfg=[];
  if ismac
    cfg.template       = {'/Users/kloosterman/Dropbox/tardis_code/MATLAB/MEG_HH_analysis/megrealign/ctf275.mat'};
    cfg.headmodel   = fullfile('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/MRI', SUBJ, [SUBJ '_hdm.mat']);
  else
    cfg.template       = {'/home/mpib/kloosterman/MATLAB/MEG_HH_analysis/megrealign/ctf275.mat'};
    cfg.headmodel   = fullfile('/home/mpib/kloosterman/projectdata/MEG2afc/MRI', SUBJ, [SUBJ '_hdm.mat']);
  end 
  cfg.inwardshift = 1;
  cfg.feedback = 'no';
  data = ft_megrealign(cfg, data);
  
  disp(' synthetic planar computation')
  cfg              = [];
  cfg.feedback     = 'no';
  cfg.method       = 'template';
  if ismac
    cfg.template     = 'ctf275_neighb.mat';
  else
    cfg.template     = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
  end
  cfg.planarmethod = 'sincos';
  cfg.channel      = {'MEG'};
  cfg.trials       = 'all';
  cfg.neighbours   = ft_prepare_neighbours(cfg, data);
  data =  ft_megplanar(cfg, data);
  
  disp('compute resplocked')
  temp = data;    data = {};
  data{1} = temp;   clear temp
  cfg=[];
  cfg.offset = -round((data{1}.trialinfo(:,5) / ft_findcfg(data{1}.cfg, 'origfs')) * ft_findcfg(data{1}.cfg, 'resamplefs'));
  data{2} = ft_redefinetrial(cfg,data{1});
  
  for ifreq = 1:2 % low high
    for idiff = 1:2
      for itrig = 1:2
        cfg = [];
        cfg.output = 'pow';
        cfg.channel = 'MEG';
        cfg.keeptapers  = 'no';
        cfg.pad = 7;
        cfg.method = 'mtmconvol';
        if ifreq == 1
          cfg.taper = 'hanning'; % low frequency-optimized analysis
          cfg.foi = 2:35;
          cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.5;
          cfg.tapsmofrq = ones(length(cfg.foi),1) .* 2;
        elseif ifreq == 2
          cfg.taper = 'dpss'; % high frequency-optimized analysis (smooth)
          cfg.foi = 36:2:100;
          cfg.t_ftimwin = ones(length(cfg.foi),1) .* 0.4;
          cfg.tapsmofrq = ones(length(cfg.foi),1) .* 8;
        end
        if itrig == 1 % stim
          cfg.toi = -0.3:0.05:0.8; %           cfg.toi = -0.3:0.05:1;
        elseif itrig == 2 % resp
          cfg.toi = -0.8:0.05:0.3; %         cfg.toi = -0.75:0.05:0.3;
        end
        cfg.trials = data{itrig}.trialinfo(:,1) == idiff;
        %         cfg.trials = 'all';
        
        %         freq{1, itrig, ifreq, idiff}{irun} = ft_combineplanar([], ft_freqanalysis(cfg, data{itrig}));
        cfg.keeptrials = 'yes';
        freqtrials = ft_combineplanar([], ft_freqanalysis(cfg, data{itrig}));
        
        freq{1, itrig, ifreq, idiff}{irun} = ft_freqdescriptives([], freqtrials);
        
        % do lateralization
        for iltr = 2:3 %  2=stim, 3=button
          % TODO normalize all trials together before latr -> starting
          % point
          cfg=[];
          cfg.avgoverrpt = 'yes';
          cfg.trials = freqtrials.trialinfo(:,1) == idiff & freqtrials.trialinfo(:,iltr) == 1; % 2=stim, 3=button
          avg_L = ft_selectdata(cfg, freqtrials);
          
          cfg.trials = freqtrials.trialinfo(:,1) == idiff & freqtrials.trialinfo(:,iltr) == 2; % 2=stim, 3=button
          avg_R = ft_selectdata(cfg, freqtrials);
          
          if itrig == 1
            disp 'get baseline from stim for resp'
            cfg=[];
            cfg.latency = [-0.25 0];
            cfg.avgovertime = 'yes';
            freq_stimbl_avg_L = ft_selectdata(cfg, avg_L);
            freq_stimbl_avg_R = ft_selectdata(cfg, avg_R);
            stimbl_powL = freq_stimbl_avg_L.powspctrm;
            stimbl_powR = freq_stimbl_avg_R.powspctrm;
          end
          disp 'repmat baseline matrix over time'
          ntim = length(avg_L.time);
          freq_stimbl_avg_L.powspctrm = repmat(stimbl_powL, 1, 1, ntim);
          freq_stimbl_avg_L.time = avg_L.time;
          freq_stimbl_avg_R.powspctrm = repmat(stimbl_powR, 1, 1, ntim);
          freq_stimbl_avg_R.time = avg_L.time;
          
          disp 'baseline correction resp with stimbaseline'
          cfg=[];
          cfg.parameter = 'powspctrm';
          cfg.operation = '((x1-x2)/x2)*100'; % relchange
          avg_L = ft_math(cfg, avg_L, freq_stimbl_avg_L); % raw
          avg_R = ft_math(cfg, avg_R, freq_stimbl_avg_R); % raw
          
          disp('compute lateralization')
          avg_L.time = avg_R.time; % otherwise lateralizedfreq might complain
          cfg=[];
          cfg.channelcmb = {'MLC11','MRC11';'MLC12','MRC12';'MLC13','MRC13';'MLC14','MRC14';'MLC15','MRC15';'MLC16','MRC16';'MLC17','MRC17';'MLC21','MRC21';'MLC22','MRC22';'MLC23','MRC23';'MLC24','MRC24';'MLC25','MRC25';'MLC31','MRC31';'MLC32','MRC32';'MLC41','MRC41';'MLC42','MRC42';'MLC51','MRC51';'MLC52','MRC52';'MLC53','MRC53';'MLC54','MRC54';'MLC55','MRC55';'MLC61','MRC61';'MLC62','MRC62';'MLC63','MRC63';'MLF11','MRF11';'MLF12','MRF12';'MLF13','MRF13';'MLF14','MRF14';'MLF21','MRF21';'MLF22','MRF22';'MLF23','MRF23';'MLF24','MRF24';'MLF25','MRF25';'MLF31','MRF31';'MLF32','MRF32';'MLF33','MRF33';'MLF34','MRF34';'MLF35','MRF35';'MLF41','MRF41';'MLF42','MRF42';'MLF43','MRF43';'MLF44','MRF44';'MLF45','MRF45';'MLF46','MRF46';'MLF51','MRF51';'MLF52','MRF52';'MLF53','MRF53';'MLF54','MRF54';'MLF55','MRF55';'MLF56','MRF56';'MLF61','MRF61';'MLF62','MRF62';'MLF63','MRF63';'MLF64','MRF64';'MLF65','MRF65';'MLF66','MRF66';'MLF67','MRF67';'MLO11','MRO11';'MLO12','MRO12';'MLO13','MRO13';'MLO14','MRO14';'MLO21','MRO21';'MLO22','MRO22';'MLO23','MRO23';'MLO24','MRO24';'MLO31','MRO31';'MLO32','MRO32';'MLO33','MRO33';'MLO34','MRO34';'MLO41','MRO41';'MLO42','MRO42';'MLO43','MRO43';'MLO44','MRO44';'MLO51','MRO51';'MLO52','MRO52';'MLO53','MRO53';'MLP11','MRP11';'MLP12','MRP12';'MLP21','MRP21';'MLP22','MRP22';'MLP23','MRP23';'MLP31','MRP31';'MLP32','MRP32';'MLP33','MRP33';'MLP34','MRP34';'MLP35','MRP35';'MLP41','MRP41';'MLP42','MRP42';'MLP43','MRP43';'MLP44','MRP44';'MLP45','MRP45';'MLP51','MRP51';'MLP52','MRP52';'MLP53','MRP53';'MLP54','MRP54';'MLP55','MRP55';'MLP56','MRP56';'MLP57','MRP57';'MLT11','MRT11';'MLT12','MRT12';'MLT13','MRT13';'MLT14','MRT14';'MLT15','MRT15';'MLT16','MRT16';'MLT21','MRT21';'MLT22','MRT22';'MLT23','MRT23';'MLT24','MRT24';'MLT25','MRT25';'MLT26','MRT26';'MLT27','MRT27';'MLT31','MRT31';'MLT32','MRT32';'MLT33','MRT33';'MLT34','MRT34';'MLT35','MRT35';'MLT36','MRT36';'MLT37','MRT37';'MLT41','MRT41';'MLT42','MRT42';'MLT43','MRT43';'MLT44','MRT44';'MLT45','MRT45';'MLT46','MRT46';'MLT47','MRT47';'MLT51','MRT51';'MLT52','MRT52';'MLT53','MRT53';'MLT54','MRT54';'MLT55','MRT55';'MLT56','MRT56';'MLT57','MRT57'};
          lrp = ft_lateralizedfreq(cfg, avg_L, avg_R);
          lrp.label = lrp.plotlabel'; % overwrite cmb label
          lrp = rmfield(lrp, 'plotlabel');
          %           lrp.plotlabel = lrp.plotlabel'; % overwrite cmb label
          
          lrp = orderfields(lrp);
          freq{iltr, itrig, ifreq, idiff}{irun} = lrp;
          freq{iltr, itrig, ifreq, idiff}{irun}.cfg = [];
        end
        clear freqtrials lrp
      end
    end
  end
  trls = [trls; trl]; %todo all trials
end

% % disp('average over runs')
% freq = cellfun(@(x) ft_freqgrandaverage([], x{:}), freq); % makes struct array
% [freq(2:end).cfg] = deal([]);  % cfg gets huge after appending data etc. only keep first

disp('keep runs as rpt')
cfg = [];
cfg.keepindividual = 'yes'
freq = cellfun(@(x) ft_freqgrandaverage(cfg, x{:}), freq); % makes struct array
[freq.dimord] = deal('rpt_chan_freq_time');  % cfg gets huge after appending data etc
[freq(2:end).cfg] = deal([]);  % cfg gets huge after appending data etc. only keep first
% [freq.trialinfo] = deal(rep_probs);  
[freq.trialinfo] = deal(trls(:,4:end));  

% lots of hassle with appendfreq
% disp 'append runs into rpt'
% cfg = [];
% cfg.appenddim = 'rpt';
% cfg.parameter = 'powspctrm';
% freq = cellfun(@(x) ft_appendfreq(cfg, x{:}), freq, 'uni', 1); 

% outfile = fullfile(PREOUT, sprintf('%s_%s_freq.mat', SUBJ, sesname));
fprintf('Saving %s\n', outfile)
save(outfile, 'freq');

%% plot lateralization
if ismac
  close all
  cfg=[];
  cfg.layout = 'CTF275.lay';
  cfg.colorbar = 'yes';
  cfg.zlim = 'maxabs';  % zeromax
  figure
  %             ft_multiplotER(cfg, timelockout{1, 1, 3}, timelockout{1, 2, 3})
  %             legend(drugleg);
  
  ft_multiplotTFR(cfg, freq{2, itrig, ifreq, idiff}{irun})
  legend({'drug - plac'});
end

%% tryout normalization and plotting
if ismac
  disp 'stimlocked baseline correction'
  cfg=[];
  cfg.baseline = [-0.25 0];
  cfg.baselinetype = 'relchange';
  %   freq_blc = arrayfun(@(x) ft_freqbaseline(cfg, x), freq(1,:,:,:));
  freq_blc_dft = arrayfun(@(x) ft_freqbaseline(cfg, x), freqdft(1,:,:,:));
  freq_blc_zap = arrayfun(@(x) ft_freqbaseline(cfg, x), freqzap(1,:,:,:));
  % contrast raw zap vs dft
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract';
  freqzapvsdft=ft_math(cfg, freqzap, freqdft)
  
  disp 'get baseline from stim for resp'
  cfg=[];
  cfg.latency = [-0.25 0];
  cfg.avgovertime = 'yes';
  freq_stimbl = arrayfun(@(x) ft_selectdata(cfg, x), freq(1,:,:,:));
  
  disp 'repmat baseline matrix over time'
  ntim = length(freq(2,1,1,1).time);
  powspctrm = arrayfun(@(x) repmat(x.powspctrm, 1, 1, ntim), freq_stimbl, 'uni', false);
  [freq_stimbl.powspctrm] = powspctrm{:};
  [freq_stimbl.time] = freq(2,:,:,:).time;
  
  disp 'baseline correction resp with stimbaseline'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = '((x1-x2)/x2)*100'; % relchange
  freq_blc(2,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), freq(2,:,:,:), freq_stimbl); % raw
  
  disp 'drug?placebo contrast'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract';
  freq_blc(:,:,3,:) = arrayfun(@(x,y) ft_math(cfg, x,y), freq_blc(:,:,1,:), freq_blc(:,:,2,:));
  
  load colormap_jetlightgray.mat
  cfg=[];
  cfg.layout = 'CTF275_helmet.mat';
  %   cfg.baseline = [-0.25 0];
  %   cfg.baselinetype = 'relchange';
  cfg.colorbar = 'yes';
  cfg.colormap = cmap;
  cfg.zlim = 'maxabs';
  cfg.hotkeys = 'yes';
  cfg.zlim = [-0.2 0.2];
  figure
  %   ft_multiplotTFR(cfg, freq_blc(1)) % freq{itrig, ifreq, idrug, idiff}
  ft_multiplotTFR(cfg, freq_blc_dft(1)) % freq{itrig, ifreq, idrug, idiff}
  title('dft')
  figure
  ft_multiplotTFR(cfg, freq_blc_zap(1)) % freq{itrig, ifreq, idrug, idiff}
  title('zap')
  cfg.zlim = 'maxabs';
  figure; ft_multiplotTFR(cfg, freqzapvsdft(1)) % freq{itrig, ifreq, idrug, idiff}
  title('zap-dft raw')
  %   figure; yyaxis left; plot(freqdft.freq, squeeze(mean(mean(freqdft.powspctrm,3))))
  %   yyaxis right; plot(freqzap.freq, squeeze(mean(mean(freqzap.powspctrm,3))))
  figure; plot(freqdft.freq, squeeze(mean(mean(freqdft.powspctrm,3)))); hold on; plot(freqzap.freq, squeeze(mean(mean(freqzap.powspctrm,3))))
  cfg=[];
  cfg.channel = 'MLO32' %'MLO32'; % {'MLO*' 'MRO*'}; %{'MLO11', 'MLO12', 'MLO13', 'MLO14', 'MLO21', 'MLO22, 'MLO23, 'MLO24, 'MLO31, 'MLO32, 'MLO33, 'MLO34, 'MLO43, 'MLO44, 'MLP31, 'MLP41, 'MLP51, 'MLP52, 'MLP53, 'MLP54, 'MLT16, 'MLT27, MRO11, MRO12, MRO13, MRO21, MRO22, MRO23, MRO31, MRO32, MRO33, MRO43, MRP31, MRP41, MRP42, MRP51, MRP52, MRP53, MRP54, MZO01, MZO02, MZP01}
  cfg.avgoverchan = 'no';
  cfg.avgovertime = 'yes';
  freqdftocc = ft_selectdata(cfg, freqdft);
  freqzapocc = ft_selectdata(cfg, freqzap);
  %   figure; yyaxis left; plot(freqdftocc.freq, freqdftocc.powspctrm); yyaxis right; plot(freqzapocc.freq, freqzapocc.powspctrm)
  figure; plot(freqdftocc.freq, freqdftocc.powspctrm); hold on; plot(freqzapocc.freq, freqzapocc.powspctrm)
  %   % same trials in datadft and zap
  %   trls = intersect(data{1}.trialinfo(:,7), datazap{1}.trialinfo(:,7));
  %   cfg=[];
  %   cfg.trials = ismember(datazap{1}.trialinfo(:,7), trls)
  %   datazap = ft_selectdata(cfg, datazap{1});
end