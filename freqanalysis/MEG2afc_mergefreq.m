function [megdat] = MEG2afc_mergefreq( baseline )
%UNTITLED6 Summary of this function goes here
%   BLC = baseline correction: raw or BLC

if nargin == 0
  baseline = 'alreadynormalized'; % BLC
end

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/'; %yesno or 2afc kloosterman
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
end

% PREIN = fullfile(basepath, 'freqzap'); % todo freq instead of freq_avg?
% PREIN = fullfile(basepath, 'freqzapruns');
% PREIN = fullfile(basepath, 'freqzap_log_BLperrun_');
% PREIN = fullfile(basepath, 'freqzap_BLperrun');

% PREIN = fullfile(basepath, 'freqzap_BLperrun_incBL');
% PREIN = fullfile(basepath, 'freqzap-plus');
% PREIN = fullfile(basepath, 'freqDFT');
% PREIN = fullfile(basepath, 'freqbandstop');
PREIN = fullfile(basepath, 'freqzapline-plus'); % sigmaincrease set to 0
cd(PREIN)

drugleg = {'drug', 'plac'};
motorleg = {'contra' 'ipsi'};

SUBJ= [1:5, 7:9, 11:21]; % all 19 subj
% SUBJ=2
% SUBJ= [1, 3:5, 7:9, 11:21]; % TODO not NK2 incomplete, 18 total
% SUBJ= [1, 3, 5, 7:9, 11:13, 15:21]; % NK2 incomplete, NK4 artifact drug-plac?, NK14 magnitude bigger responses
% SUBJ= [1, 3, 5, 7:9, 11:13, 15:21]; % NK2 incomplete, NK4 artifact drug-plac?, NK14 magnitude bigger responses

% SUBJ= [1, 3:5, 7:9, 11:21]; % NK2 incomplete
% SUBJ= [1, 4:5, 7:9, 11:21]; % NK2 incomplete, NK3 has nans, FIX

% behavpath = '/Users/kloosterman/gridmaster2012/MATLAB/MEG_HH_analysis/behavior/historybias';
% load(fullfile(behavpath, 'corrstat_betalat.mat')); %corrstat comes out: for mask!
% mask = corrstat.negclusterslabelmat == 1;
% rp_lat_runs = nan(length(SUBJ), 8, 2,2,2,2); %dimord subj runs drug ses diff repprobvslat

megdat = []; % output
megdat.SUBJ = SUBJ;
megdat.PREOUT = fullfile(PREIN, 'plots');  mkdir(megdat.PREOUT);
megdat.dimord = 'drug_mod_trig_freq_diff';

nsub = length(SUBJ);
freqall = {};
for isub = 1:nsub
  for idrug = 1:2
    sesfiles = dir(sprintf('NK%d_%s*.mat', SUBJ(isub), drugleg{idrug}));
    for is = 1:length(sesfiles)
      fprintf(sesfiles(is).name)
      %       if contains(PREIN, 'freqzap_BLperrun') || contains(PREIN, 'freqzap_log_BLperrun_') || contains(PREIN, 'freqDFT') || contains(PREIN, 'plus')
      try
        load(sesfiles(is).name, 'freq_avg')
        if ~isfield(freq_avg(1), 'powspctrm')
          fprintf('not found\n')
          %         continue
          load(sesfiles(1).name, 'freq_avg'); %NK2 misses 1 ses
        end
        freq=freq_avg; clear freq_avg
      catch
        %       else % TODO fix this hack
        load(sesfiles(is).name, 'freq')
        if ~isfield(freq(1), 'powspctrm')
          fprintf('not found\n')
          %         continue
          load(sesfiles(1).name); %NK2 misses 1 ses
        end
      end
      
      fprintf(' %d runs\n', size(freq(1).powspctrm,1))
      
      if contains(freq(1).dimord, 'rpt')
        %         nrpt = size(freq(1).powspctrm,1);
        %         for idiff = 1:2
        %           iltr = 3; itrig = 2; ifreq = 1; % lat repprob
        %           freqoi = squeeze(freq(iltr, itrig, ifreq, idiff));
        %           rp_lat_runs(isub, 1:nrpt, idrug, is, idiff, 1) = ...
        %             mean(freqoi.powspctrm(:,mask),2); % meg latr
        %           rp_lat_runs(isub, 1:nrpt, idrug, is, idiff, 2) = ...
        %             freqoi.trialinfo; % repprob
        %         end
        freq = arrayfun(@(x) ft_freqdescriptives([], x), freq); % avg over runs
      end
      %       figure; ft_singleplotTFR([], freq(10))
      freqall(isub,idrug,is,:,:,:,:) = num2cell(freq); % freq{iltr, itrig, ifreq, idiff}
      
    end
  end
end

disp 'collect subjects'
% cfg=[]; % TODO use no loops
% cfg.parameter = 'powspctrm';
% cfg.keepindividual = 'yes';
% freqtmp2 = cellfun(@(x) ft_freqgrandaverage(cfg, x), freqall(:,:));
temp = {};
% dimord:        freq{2, itrig, ifreq, idiff}{irun} = lrp;

nmod = size(freqall, 4); % latr wrt prevresp not there for freqzap
cfg=[];
cfg.keepindividual = 'yes';
for idrug = 1:2
  for imod = 1:nmod % modulation or lateralization wrt stim, resp
    for is = 1:2
      for itrig = 1:2
        for ifreq = 1:2
          for idiff = 1:2
            temp{idrug, is, imod, itrig, ifreq, idiff} =  ft_freqgrandaverage(cfg, freqall{:, idrug, is, imod, itrig, ifreq, idiff } );
          end
        end
      end
    end
  end
end
disp 'avg motor ses'
% freqall = squeeze(cellfun(@(x,y) ft_freqgrandaverage([], x,y), temp(:,1,:,:,:,:), temp(:,2,:,:,:,:)));

temp = cellfun(@(x,y) ft_freqgrandaverage([], x,y), temp(:,1,:,:,:,:), temp(:,2,:,:,:,:));
z=size(temp);
freqall = reshape(temp, z([1,3:end])); 
clear temp

% megdat.rp_lat_runs = rp_lat_runs;
% megdat.rp_lat_runs_dimord = 'subj_runs_drug_ses_diff_repprobvslat';

disp 'normalization' % TODO complete this
if strcmp(baseline, 'raw') | strcmp(baseline, 'alreadynormalized')
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract';
  disp 'easy-hard avg'
  freqall(:,:,:,:,3) = squeeze(arrayfun(@(x,y) ft_freqgrandaverage([], x,y), freqall(:,:,:,:,1), freqall(:,:,:,:,2)));
  disp 'easy-hard contrast'
  
  %   freqall(:,:,:,:,3) = arrayfun(@(x,y) ft_math(cfg, x,y), freqall(:,:,:,:,1), freqall(:,:,:,:,2)); % WTF???!!
  
  freqall(:,:,:,:,4) = arrayfun(@(x,y) ft_math(cfg, x,y), freqall(:,:,:,:,1), freqall(:,:,:,:,2));
  disp 'drug-plac contrast'
  freqall(4,:,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x,y), freqall(1,:,:,:,:), freqall(2,:,:,:,:));
  
  megdat.freq = freqall;
  
else
  disp 'stimlocked baseline correction'
  cfg=[];
  cfg.baseline = [-0.25 0];
  cfg.baselinetype = 'relchange';% only BLC mod, not latr
  freq_blc = arrayfun(@(x) ft_freqbaseline(cfg, x), freqall(:,1,1,:,:)); % mod, stimlocked
  
  disp 'get baseline from stim for resp'
  cfg=[];
  cfg.latency = [-0.25 0];
  cfg.avgovertime = 'yes';
  freq_stimbl = arrayfun(@(x) ft_selectdata(cfg, x), freqall(:,1,1,:,:));
  
  disp 'repmat baseline matrix over time'
  ntim = length(freqall(1,1,2,1,1).time);
  powspctrm = arrayfun(@(x) repmat(x.powspctrm, 1, 1, 1, ntim), freq_stimbl, 'uni', false);
  [freq_stimbl.powspctrm] = powspctrm{:};
  clear powspctrm
  [freq_stimbl.time] = freqall(:,1,2,:,:).time;
  
  disp 'baseline correction resp with stimbaseline'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract'; % relchange in 2 steps, faster than  cfg.operation = '((x1-x2)/x2)*100';
  freq_blc(:,1,2,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), freqall(:,1,2,:,:), freq_stimbl); % raw
  cfg.operation = 'divide'; % relchange
  freq_blc(:,1,2,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), freq_blc(:,1,2,:,:), freq_stimbl); % raw
  clear freq_stimbl
  
  % add latr which is already BL corrected
  freq_blc(:,2:3,:,:,:) = freqall(:,2:3,:,:,:);
  
  disp 'drug avg'
  freq_blc(3,:,:,:,:) = arrayfun(@(x,y) ft_freqgrandaverage([], x,y), freq_blc(1,:,:,:,:), freq_blc(2,:,:,:,:));
  disp 'diff avg'
  freq_blc(:,:,:,:,3) = arrayfun(@(x,y) ft_freqgrandaverage([], x,y), freq_blc(:,:,:,:,1), freq_blc(:,:,:,:,2));
  
  disp 'drug-plac contrast'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract';
  freq_blc(4,:,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x,y), freq_blc(1,:,:,:,:), freq_blc(2,:,:,:,:));
  disp 'easy-hard contrast'
  freq_blc(:,:,:,:,4) = arrayfun(@(x,y) ft_math(cfg, x,y), freq_blc(:,:,:,:,1), freq_blc(:,:,:,:,2));
  
  megdat.freq = freq_blc;
end

disp 'add behavior'
load '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/behav/behavstruct.mat';
megdat.behavior = behavior; % check dropped subjects in behav and meg

megdat.megleg = {'mod' 'latrstim' 'latrresp' 'latrprevresp' 'baseline'};

%%
testplot=1
if testplot
  % multiplot for testing if all looks ok
%   close all
  idrug = 4; idiff = 3; imod = 1; itrig = 1; ifreq = 1;
  plotsinglesubj = 1;
  
  load colormap_jetlightgray.mat
  cfg=[];
  if imod == 1
    cfg.layout = 'CTF275_helmet.mat';
  else
    cfg.layout = 'CTF275_helmet_latr.mat';
  end
  %   cfg.layout = 'CTF275.lay';
  %   cfg.baseline = [-0.25 0];
  %   cfg.baselinetype = 'relchange';
  cfg.colorbar = 'yes';
  cfg.colormap = cmap;
  if ifreq == 2
    cfg.ylim = [35 100];
  end
  cfg.zlim = 'maxabs';
  % cfg.zlim = [-20 20];
  cfg.hotkeys = 'yes';
  if itrig == 1; cfg.xlim = [-0.2 0.6]; else cfg.xlim = [-0.8 0.2]; end
  
  % cfg.maskparameter = 'mask';
  % cfg.maskalpha = 0.25;
  % f = figure;
  % f.Position = [   680   444   781   654];
  freq = megdat.freq(idrug, imod, itrig, ifreq, idiff);
  subjfreq=freq;
  if plotsinglesubj
    for isub = 1%:size(freq.powspctrm,1)
      subjfreq.powspctrm = freq.powspctrm(isub,:,:,:);
      figure; title(isub)
      try
        ft_multiplotTFR(cfg, subjfreq) % idrug, itrig, ifreq, idiff
      catch; end
    end
  else
    ft_multiplotTFR(cfg, freq) % idrug, itrig, ifreq, idiff
  end
  %   ft_multiplotTFR(cfg, freqall(3,1,2,1)) % idrug, itrig, ifreq, idiff
  % title(stat{idrug, imod, itrig, ifreq, idiff}.posclusters(1).prob)
  title(sprintf('Drug%d mod%d diff%d', idrug, imod, idiff))
  
end
