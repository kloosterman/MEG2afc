function [megdat] = MEG2afc_mergedva( baseline )
%UNTITLED6 Summary of this function goes here
%   BLC = baseline correction: raw or BLC

if ismac
  basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
end

PREIN = fullfile(basepath, 'dvalock');
cd(PREIN)

drugleg = {'drug', 'plac'};
motorleg = {'contra' 'ipsi'};

% dvalockall = struct([]);
SUBJ= [1:5, 7:9, 11:21]; % all 19 subj
SUBJ= [1, 3:5, 7:9, 11:21]; % TODO not NK2 incomplete, 18 total
SUBJ= [1, 3, 5, 7:9, 11:13, 15:21]; % NK2 incomplete, NK4 artifact drug-plac?, NK14 magnitude bigger responses
SUBJ= [1, 3, 5, 7:9, 11:13, 15:21]; % NK2 incomplete, NK4 artifact drug-plac?, NK14 magnitude bigger responses

megdat = []; % output
megdat.SUBJ = SUBJ;
megdat.PREOUT = fullfile(PREIN, 'plots');

% nsub = 3 
nsub = length(SUBJ);
dvalockall = cell(nsub,2,2,2,3);
for isub = 1:nsub
  for idva = 1:2
    for idrug = 1:2
      sesfiles = dir(sprintf('NK%d_%s*.mat', SUBJ(isub), drugleg{idrug}));
      
      for is = 1:length(sesfiles)
        disp(sesfiles(is).name)
        load(sesfiles(is).name) %dvalock{itrig, idiff}
        if isempty(dvalock)
          dvalock = dvalockall{isub,idrug,1,:,:,:};
        end
        dvalockall(isub,idrug,is,:,:) = dvalock;  % dvalock{idva, idiff} 
        %       dvalockall(isub,idrug,is,:,:,:).SUBJ = SUBJ;
      end
    end
  end
end

disp 'collect subjects'
% cfg=[]; % TODO use no loops
% cfg.parameter = 'powspctrm';
% cfg.keepindividual = 'yes';
% freqtmp2 = cellfun(@(x) ft_timelockgrandaverage(cfg, x), dvalockall(:,:));
temp = {};
cfg=[];
cfg.keepindividual = 'yes';
for idrug = 1:2
  for is = 1:2
    for idva = 1:2
      for idiff = 1:3
        temp{idrug, is, idva,idiff} =  ft_timelockgrandaverage(cfg, dvalockall{:, idrug, is, idva, idiff } );
      end
    end
  end
end

disp 'avg motor ses'
cfg=[];
cfg.parameter = 'individual';
dvalockall = squeeze(cellfun(@(x,y) ft_timelockgrandaverage(cfg, x,y), temp(:,1,:,:), temp(:,2,:,:))); %

dvalockall = rmfield(dvalockall, {'var', 'dof'});

if strcmp(baseline, 'raw') %% normalization
  disp 'easy-hard contrast'
  cfg=[];
  cfg.parameter = 'avg';
  cfg.operation = 'subtract';
  dvalockall(:,:,4) = arrayfun(@(x,y) ft_math(cfg, x,y), dvalockall(:,:,1), dvalockall(:,:,2));
  disp 'drug-plac contrast'
  dvalockall(3,:,:) = arrayfun(@(x,y) ft_math(cfg, x,y), dvalockall(1,:,:), dvalockall(2,:,:));
%   disp 'easy-hard avg'
%   dvalockall(:,:,:,4) = squeeze(arrayfun(@(x,y) ft_timelockgrandaverage([], x,y), dvalockall(:,:,:,1), dvalockall(:,:,:,2)));
  
  megdat.dvalock = dvalockall;
  megdat.dimord = 'drug_dva_diff';

else
  disp 'stimlocked baseline correction'
  cfg=[];
  cfg.baseline = [-0.25 0];
  cfg.baselinetype = 'relchange';
  freq_blc = arrayfun(@(x) ft_freqbaseline(cfg, x), dvalockall(:,1,:,:));
  
  disp 'get baseline from stim for resp'
  cfg=[];
  cfg.latency = [-0.25 0];
  cfg.avgovertime = 'yes';
  freq_stimbl = arrayfun(@(x) ft_selectdata(cfg, x), dvalockall(:,1,:,:));
  
  disp 'repmat baseline matrix over time'
  ntim = length(dvalockall(1,2,1,1).time);
  powspctrm = arrayfun(@(x) repmat(x.powspctrm, 1, 1, 1, ntim), freq_stimbl, 'uni', false);
  [freq_stimbl.powspctrm] = powspctrm{:};
  [freq_stimbl.time] = dvalockall(:,2,:,:).time;
  
  disp 'baseline correction resp with stimbaseline'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract'; % relchange in 2 steps, faster than  cfg.operation = '((x1-x2)/x2)*100';
  freq_blc(:,2,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), dvalockall(:,2,:,:), freq_stimbl); % raw
  cfg.operation = 'divide'; % relchange
  freq_blc(:,2,:,:) = arrayfun(@(x,y) ft_math(cfg, x, y), freq_blc(:,2,:,:), freq_stimbl); % raw
  
  % Contrasts
  disp 'easy-hard contrast'
  cfg=[];
  cfg.parameter = 'powspctrm';
  cfg.operation = 'subtract';
  freq_blc(:,:,:,3) = arrayfun(@(x,y) ft_math(cfg, x,y), freq_blc(:,:,:,1), freq_blc(:,:,:,2));
  disp 'drug-plac contrast'
  freq_blc(3,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x,y), freq_blc(1,:,:,:), freq_blc(2,:,:,:));
  disp 'easy-hard avg'
  freq_blc(:,:,:,4) = squeeze(arrayfun(@(x,y) ft_timelockgrandaverage([], x,y), freq_blc(:,:,:,1), freq_blc(:,:,:,2)));
  
  megdat.dvalock = freq_blc;
  
end
