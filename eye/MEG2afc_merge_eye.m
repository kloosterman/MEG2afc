function [timelock_blc, eyeall] = MEG2afc_merge_eye(  )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if ismac
  basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
else
  basepath = '/home/beegfs/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
end

PREIN = fullfile(basepath, 'eye', 'timelock');
cd(PREIN)

drugleg = {'drug', 'plac'};
motorleg = {'contra' 'ipsi'};

% eyeall = struct([]);
SUBJ= [1:5, 7:9, 11:21]; % all
SUBJ= [1, 3:5, 7:9, 11:21]; % TODO not NK2 incomplete
% nsub = 3 
nsub = length(SUBJ);
eyeall = cell(nsub,2,2,2,2);
for isub = 1:nsub
  for idrug = 1:2
    sesfiles = dir(sprintf('NK%d_%s*.mat', SUBJ(isub), drugleg{idrug}));
%     timelocktmp = struct([]);
    for is = 1:length(sesfiles)
      disp(sesfiles(is).name)
      load(sesfiles(is).name)
      timelock = removefields(timelock, {'var', 'dof'});
      if isempty(timelock)
        timelock = eyeall{isub,idrug,1,:,:,:};
      end
      eyeall(isub,idrug,is,:,:) = num2cell(timelock);
%       eyeall(isub,idrug,is,:,:,:).SUBJ = SUBJ;
    end
  end
end

disp 'collect subjects'
% cfg=[]; % TODO use no loops
% cfg.parameter = 'avg';
% cfg.keepindividual = 'yes';
% timelocktmp2 = cellfun(@(x) ft_timelockgrandaverage(cfg, x), eyeall(:,:));
temp = {};
cfg=[];
cfg.keepindividual = 'yes';
for idrug = 1:2
  for is = 1:2
    for itrig = 1:2
      for idiff = 1:2
        temp{idrug, is, itrig, idiff} =  ft_timelockgrandaverage(cfg, eyeall{:, idrug, is, itrig, idiff } );
      end
    end
  end
end

disp 'avg motor ses'
cfg=[];
cfg.parameter = 'individual';
eyeall = squeeze(cellfun(@(x,y) ft_timelockgrandaverage(cfg, x,y), temp(:,1,:,:), temp(:,2,:,:)));

%% normalization and plotting
% warning('off','MATLAB:lang:cannotClearExecutingFunction');
disp 'stimlocked baseline correction'
cfg=[];
cfg.baseline = [-0.25 0];
timelock_blc = arrayfun(@(x) ft_timelockbaseline(cfg, x), eyeall(:,1,:));

disp 'get baseline from stim for resp'
cfg=[];
cfg.latency = [-0.25 0];
cfg.avgovertime = 'yes';
timelock_stimbl = arrayfun(@(x) ft_selectdata(cfg, x), eyeall(:,1,:));

disp 'repmat baseline matrix over time'
ntim = length(eyeall(1,2,1).time);
avg = arrayfun(@(x) repmat(x.avg, 1, 1, ntim), timelock_stimbl, 'uni', false);
[timelock_stimbl.avg] = avg{:};
[timelock_stimbl.time] = eyeall(:,2,:).time;

      eyeall = removefields(eyeall, {'var', 'dof'});
      timelock_stimbl = removefields(timelock_stimbl, {'var', 'dof'});
      timelock_blc = removefields(timelock_blc, {'var', 'dof'});

disp 'baseline correction resp with stimbaseline'
cfg=[];
cfg.parameter = 'avg';
cfg.operation = 'subtract'; % relchange in 2 steps, faster than  cfg.operation = '((x1-x2)/x2)*100';
timelock_blc(:,2,:) = arrayfun(@(x,y) ft_math(cfg, x, y), eyeall(:,2,:), timelock_stimbl); % raw
% cfg.operation = 'divide'; % relchange
% timelock_blc(:,2,:) = arrayfun(@(x,y) ft_math(cfg, x, y), timelock_blc(:,2,:), timelock_stimbl); % raw


%%
disp 'easy-hard contrast'
cfg=[];
cfg.parameter = 'avg';
cfg.operation = 'subtract';
timelock_blc(:,:,3) = arrayfun(@(x,y) ft_math(cfg, x,y), timelock_blc(:,:,1), timelock_blc(:,:,2));
eyeall(:,:,3) = arrayfun(@(x,y) ft_math(cfg, x,y), eyeall(:,:,1), eyeall(:,:,2));

disp 'drug-plac contrast'
timelock_blc(3,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x,y), timelock_blc(1,:,:,:), timelock_blc(2,:,:,:));
eyeall(3,:,:,:) = arrayfun(@(x,y) ft_math(cfg, x,y), eyeall(1,:,:,:), eyeall(2,:,:,:));

disp 'easy-hard avg'
timelock_blc(:,:,4) = squeeze(arrayfun(@(x,y) ft_timelockgrandaverage([], x,y), timelock_blc(:,:,1), timelock_blc(:,:,2)));
eyeall(:,:,4) = squeeze(arrayfun(@(x,y) ft_timelockgrandaverage([], x,y), eyeall(:,:,1), eyeall(:,:,2)));

