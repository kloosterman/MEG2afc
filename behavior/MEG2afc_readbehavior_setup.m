function [behavior] = MEG2afc_readbehavior_setup()
% read in behav data from eye data

if ismac
  %   basepath = '/Users/kloosterman/beegfs/projectdata/MEG2afc/'; %yesno or 2afc
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/'; %yesno or 2afc
  backend = 'local';
  compile = 'no';
else
  basepath = '/home/mpib/kloosterman/projectdata/MEG2afc/'; %yesno or 2afc
      backend = 'slurm';
%   backend = 'torque';
  %       backend = 'local';
  %   compile = 'yes';
  compile = 'no';
end
timreq = 5; %in minutes per run
memreq = 3000; % in MB

PREIN = fullfile(basepath, 'preproczap');
PREOUT = fullfile(basepath, 'behav');  
mkdir(PREOUT)

saveddm_mat = 0; % save csv for hddm 
runontardis = 1; % run it or load behav per session from file

% subject issues: 
% NK1: high d' (pilot)  KEEP
% NK2: session missing (fixable??) KEEP
% NK10: bad DDM fits (weird RT distributions, outlier driftbias) DROP

SUBJ = [1:5, 7:9, 11:21]; % all
% SUBJ = [2]; % 
SUBJ_idx= [0:18]; % corresponding Python counting
% SUBJ_idx=1;
% SUBJ= [1, 3:5, 7:9, 11:21]; % drop NK2, misses session
% SUBJ_idx= [0, 2:18]; % corresponding Python counting, drop NK2, misses session

nsub = length(SUBJ)

% sesdirs = {'A' 'B' 'C' 'D'};
% sesnames = {'plac_ipsi', 'drug_ipsi', 'plac_contra', 'drug_contra'};
drugs = {'drug', 'plac'};
ses = {'contra', 'ipsi'};

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
cfglist = {};
for isub = 1:nsub
  for idrug = 1:2
    for ises = 1:2      
      cfg.PREIN = PREIN;
      cfg.SUBJ = sprintf('NK%d', SUBJ(isub));
      cfg.PREOUT = PREOUT;
      cfg.sesname = sprintf('%s_%s', drugs{idrug}, ses{ises}); %sesnames{ises};
      cfg.outfile = fullfile(PREOUT, sprintf('%s_%s_%s_behav.mat', cfg.SUBJ, drugs{idrug}, ses{ises}));
      disp(cfg.outfile)
      if runontardis
        cfglist{isub, idrug, ises} = cfg;
      else
        try
          temp = load(cfg.outfile);
        catch 
          fprintf('%s not found\n', cfg.outfile)
          continue
        end
        behav(isub, idrug, ises) = temp.behav;
      end
    end
  end
end

if runontardis
  fprintf('Running MEG2afc_readbehavior for %d cfgs\n', numel(cfglist))
  
  disp 'todo debug running this on tardis'
  
  if strcmp(backend, 'slurm')
    options = '-D. -c2'; % --gres=gpu:1
  else
    options =  '-l nodes=1:ppn=1'; % torque %-q testing or gpu
  end
  
  setenv('TORQUEHOME', 'yes')
  mkdir('~/qsub'); cd('~/qsub');
  if strcmp(compile, 'yes')
    fun2run = qsubcompile(@MEG2afc_readbehavior, 'toolbox', {'signal', 'stats'}); %
  else
    fun2run = @MEG2afc_readbehavior;
  end
  
  if strcmp(backend, 'local')
    behav = cellfun(fun2run, cfglist);
  else
    behav = qsubcellfun(fun2run, cfglist, 'memreq', memreq, 'timreq', timreq*60, 'stack', 1, ...
      'StopOnError', false, 'backend', backend, 'options', options);
  end
end
disp 'save raw behav struct array'
save(fullfile(PREOUT, 'rawbehav.mat'), 'behav')

%%
if saveddm_mat
  disp 'save csv for HDDM all subj together'
  ddmmat = vertcat(behav.ddmmat_runs);
  outfile = 'ddmdat_histbias.csv';
  outdir = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM';
  outpath = fullfile(outdir, outfile);
  fid = fopen(outpath, 'w') ;
  fprintf(fid, 'subj_idx,stimulus,response,prevresp,correct,rt,drug,simon,run_nr,difficulty,trlctr,megused\n');
  fclose(fid);
  dlmwrite(outpath, ddmmat, '-append', 'precision', 10)
  disp 'now run again without saving ddmmat'
  behavior = []; % struct with all behavior
  
  disp 'save csv for HDDM per subj'
  for isub=1:19
    ddmmat = vertcat(behav(isub,:,:).ddmmat_runs);
    ddmmat = ddmmat(:,2:end); % drop subj index col
    ddmmat(:,8) = ddmmat(:,8) - 1; % let it range from 0-X
    outfile = sprintf('subj%d_ddmdat.csv', isub-1);
    outdir = fullfile(PREOUT, 'HDDM');
    mkdir(outdir)
    outpath = fullfile(outdir, outfile);
    fid = fopen(outpath, 'w') ;
    % note run_nr has become subj_idx
    fprintf(fid, 'stimulus,response,prevresp,correct,rt,drug,simon,subj_idx,difficulty,trlctr,megused\n');
    fclose(fid);
    dlmwrite(outpath, ddmmat, '-append', 'precision', 10)
    disp 'now run again without saving ddmmat'
    behavior = []; % struct with all behavior
  end
  
  return
end

%% collect all behavior in nice struct
behavior = []; % struct with all behavior
behavior.SUBJ = SUBJ;
behavior.SUBJ_idx = SUBJ_idx;

disp 'get bmeas dprime and criterion'; disp 'get RT'
disp 'get p_repeat separately for L and R repeats'
bmeas = {'dprime' 'criterion' 'button_bias' 'p_repeatbalanced' 'RT' 'RTsd' 'ntrials'}; % p_repeatbalanced dim5 is LR
for im = 1:length(bmeas)
  behavior.(bmeas{im}) = reshape([behav.(bmeas{im})], 9, 2, nsub,2,2); % dims: runs diff subj drug motor
  behavior.(bmeas{im}) = permute(behavior.(bmeas{im}), [3 1 4 5 2]); % dims: subj runs drug motor diff
  behavior.([bmeas{im} 'dimord']) = 'subj_runs_drug_motor_diff';
  behavior.(bmeas{im})(:,:,4,:,:) = behavior.(bmeas{im})(:,:,1,:,:) - behavior.(bmeas{im})(:,:,2,:,:);
  behavior.(bmeas{im})(:,:,:,:,3) = mean(behavior.(bmeas{im}), 5); % avg over diff
  behavior.(bmeas{im})(:,:,:,3,:) = nanmean(behavior.(bmeas{im}), 4); % avg over motor ses
  behavior.(bmeas{im})(:,10,:,:,:) = nanmean(behavior.(bmeas{im})(:,1:8,:,:,:), 2); % avg over runs, 9 has  collapsed
end

%%
disp 'make RT distributions'
edges = 0:0.05:2.5;
behavior.RThist = nan(nsub, 2,2, length(edges)-1);
for isub = 1:nsub
  for idrug = 1:2
    for ises = 1:2
      if ~isempty(behav(isub, idrug, ises).ddmmat_runs)
        behavior.RThist(isub, idrug, ises,:) = histcounts(behav(isub, idrug, ises).ddmmat_runs(:,6), edges, 'Normalization', 'probability');
      end
    end
  end
end
behavior.RThist(:,:,3,:) = nanmean(behavior.RThist(:,:,1:2,:),3);
behavior.RThist(:,4,:,:) = behavior.RThist(:,1,:,:) - behavior.RThist(:,2,:,:);
behavior.RThistedges = edges;

%%
disp 'get P(repeat)'
bmeas = {'p_repeatunbalanced' 'basepupil' 'bpm' };
for im = 1:length(bmeas)
  behavior.(bmeas{im}) = reshape([behav.(bmeas{im})], 9, nsub,2,2); % dims: runs subj drug motor
  behavior.(bmeas{im}) = permute(behavior.(bmeas{im}), [2 1 3 4]); % dimord subj runs drug motor
  behavior.(bmeas{im})(:,:,4,:) = behavior.(bmeas{im})(:,:,1,:) - behavior.(bmeas{im})(:,:,2,:);
  behavior.(bmeas{im})(:,:,:,3) = nanmean(behavior.(bmeas{im}), 4); % avg over motor ses
  behavior.(bmeas{im})(:,10,:,:) = nanmean(behavior.(bmeas{im}), 2); % avg over runs,  9 has  collapsed
  behavior.([bmeas{im} 'dimord']) = 'subj_runs_drug_motor';
end

%% load perrun DDM fits
ddmpath = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/';
pars = {'a' 't' 'v' 'z' 'dc'};
drugs = {'atx', 'plac'};
ses = {'contra', 'ipsi'};
diffs = {'easy' 'hard'};
prevresps = {'Lprev' 'Rprev'};

model_names = { 'ddm_histbias_perrun' 'ddm_histbias_perrun_ol' 'ddm_acc_perrun' 'ddm_acc_perrun_ol'}; % 'ddm_histbias_perses'
for imodel = 1:length(model_names)
  model_name = model_names{imodel};
  behavior.(model_name) = [];
  behavior.(model_name).dimord = 'subj_runs_drug_motor_difforprevresp';
  isub=0;
  for subind = SUBJ_idx
    isub=isub+1;
    for idrug = 1:2
      for ises = 1:2
        warning off
        try          
          ddmfit = readtable(fullfile(ddmpath, sprintf('params_%s_subj%d_%s_%s.csv', model_name, subind, drugs{idrug}, ses{ises}) ));
        end
        warning on
        varnames = ddmfit.Properties.VariableNames';
        for ipar = 1:length(pars)
%           behavior.(model_name).(pars{ipar})(isub,1:9,1:4,1:4) = NaN(1,9,4,4);
          cols = startsWith(varnames, pars{ipar}) & ~contains(varnames, 'trans'); % & contains(varnames, drugs{idrug}) & contains(varnames, ses{ises})
          if ~any(cols); continue;  end
          if any(contains(varnames(cols), 'easy')) %  easy vs hard (v)
            allcols = cols;
            for idiff=1:2 % put diff dim last
              cols = allcols & contains(varnames, diffs{idiff});
              dat = table2array(ddmfit(:,cols));
              temp = NaN(1,8);
              temp(1,1:length(dat)) = dat;
              behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises,idiff) = temp;
            end
          elseif any(contains(varnames(cols), 'prev'))
            allcols = cols;
            for iprev=1:2 % put prevresp dim last    
              cols = allcols & contains(varnames, prevresps{iprev});
              dat = table2array(ddmfit(:,cols));
              temp = NaN(1,8);
              temp(1,1:length(dat)) = dat;
              behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises,iprev) = temp;
            end
          else % no diff or prevResp
            dat = table2array(ddmfit(:,cols));
            temp = NaN(1,8);
            temp(1,1:length(dat)) = dat;
            behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises) = temp;
          end
        end
      end
    end
  end
  for ipar = 1:length(pars)
    if isfield(behavior.(model_name), pars{ipar})
      if length(size(behavior.(model_name).(pars{ipar}))) == 4
        behavior.(model_name).(pars{ipar})(:,9,:,:) = nanmean(behavior.(model_name).(pars{ipar})(:,:,:,:),2);
        behavior.(model_name).(pars{ipar})(:,:,3,:) = mean(behavior.(model_name).(pars{ipar})(:,:,:,:),3);
        behavior.(model_name).(pars{ipar})(:,:,:,3) = mean(behavior.(model_name).(pars{ipar})(:,:,:,:),4);
        behavior.(model_name).(pars{ipar})(:,:,4,:) = behavior.(model_name).(pars{ipar})(:,:,1,:,:) - behavior.(model_name).(pars{ipar})(:,:,2,:,:);
        behavior.(model_name).(pars{ipar})(:,:,:,4) = behavior.(model_name).(pars{ipar})(:,:,:,1,:) - behavior.(model_name).(pars{ipar})(:,:,:,2,:);
      elseif length(size(behavior.(model_name).(pars{ipar}))) == 5
        behavior.(model_name).(pars{ipar})(:,9,:,:,:) = nanmean(behavior.(model_name).(pars{ipar})(:,:,:,:,:),2);
        behavior.(model_name).(pars{ipar})(:,:,3,:,:) = mean(behavior.(model_name).(pars{ipar})(:,:,:,:,:),3);
        behavior.(model_name).(pars{ipar})(:,:,:,3,:) = mean(behavior.(model_name).(pars{ipar})(:,:,:,:,:),4);
        behavior.(model_name).(pars{ipar})(:,:,:,:,3) = mean(behavior.(model_name).(pars{ipar})(:,:,:,:,:),5);
        behavior.(model_name).(pars{ipar})(:,:,4,:,:) = behavior.(model_name).(pars{ipar})(:,:,1,:,:) - behavior.(model_name).(pars{ipar})(:,:,2,:,:);
        behavior.(model_name).(pars{ipar})(:,:,:,4,:) = behavior.(model_name).(pars{ipar})(:,:,:,1,:) - behavior.(model_name).(pars{ipar})(:,:,:,2,:);
        behavior.(model_name).(pars{ipar})(:,:,:,:,4) = behavior.(model_name).(pars{ipar})(:,:,:,:,1) - behavior.(model_name).(pars{ipar})(:,:,:,:,2);
      end
    end
  end
  if contains(model_name, 'histbias')
    behavior.(model_name).histshift_z = behavior.(model_name).z(:,:,:,:,2) - behavior.(model_name).z(:,:,:,:,1); % histshift = prevR-prevL
    behavior.(model_name).histshift_dc = behavior.(model_name).dc(:,:,:,:,2) - behavior.(model_name).dc(:,:,:,:,1); % histshift = prevR-prevL
  end
end

%% within-subject HDDM
disp 'load ddm fits basic within-subject HDDM model runs'
model_names = {'withinsubj_HDDM' }; % chi_accuracy_basic_runs
pars = {'a' 't' 'v' }; % 'z' 'dc'
drugs = {'atx', 'plac'};
ses = {'contra', 'ipsi'};
diffs = {'easy' 'hard'};
prevresps = {'Lprev' 'Rprev'};
ddmpath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/behav/HDDM';
for imodel = 1:length(model_names)
  model_name = model_names{imodel};
  behavior.(model_name) = [];
  behavior.(model_name).dimord = 'subj_runs_drug_motor_difforprevresp';
  isub=0;
  for subind = SUBJ_idx
    isub=isub+1;
    ddmfit = readtable(fullfile(ddmpath, sprintf('subj%d_ddmparams.csv', subind) ));
    varnames = ddmfit.Properties.VariableNames';
    for ipar = 1:length(pars)
      for idrug = 1:2
        for ises = 1:2
          cols = startsWith(varnames, pars{ipar}) & contains(varnames, drugs{idrug}) ...
            & contains(varnames, ses{ises}) & ~contains(varnames, 'trans');
          if ~any(cols); continue;  end
          if any(contains(varnames(cols), 'easy')) %  easy vs hard (v)
            allcols = cols;
            for idiff=1:2 % put diff dim last
              temp = NaN(1,8);
              cols = allcols & contains(varnames, diffs{idiff});
              temp(1,1:length(find(cols))) = table2array(ddmfit(1,cols));
              behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises,idiff) = temp;
            end
          elseif any(contains(varnames(cols), 'prev')) %  easy vs hard (v)
            allcols = cols;
            for iprev=1:2 % put prevresp dim last
              temp = NaN(1,8);
              cols = allcols & contains(varnames, prevresps{iprev});
              temp(1,1:length(find(cols))) = table2array(ddmfit(1,cols));
              behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises,iprev) = temp;
            end
          else % no diff or prevResp
            temp = NaN(1,8);
            temp(1,1:length(find(cols))) = table2array(ddmfit(1,cols));
            behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises) = temp;
          end
        end
      end
      if isfield(behavior.(model_name), pars{ipar})
        behavior.(model_name).(pars{ipar})(isub,9,:,:,:) = nanmean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),2);
        behavior.(model_name).(pars{ipar})(isub,:,3,:,:) = mean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),3);
        behavior.(model_name).(pars{ipar})(isub,:,:,3,:) = mean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),4);
        behavior.(model_name).(pars{ipar})(isub,:,:,:,3) = mean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),5);
        behavior.(model_name).(pars{ipar})(isub,:,4,:,:) = behavior.(model_name).(pars{ipar})(isub,:,1,:,:) - behavior.(model_name).(pars{ipar})(isub,:,2,:,:);
        behavior.(model_name).(pars{ipar})(isub,:,:,4,:) = behavior.(model_name).(pars{ipar})(isub,:,:,1,:) - behavior.(model_name).(pars{ipar})(isub,:,:,2,:);
        behavior.(model_name).(pars{ipar})(isub,:,:,:,4) = behavior.(model_name).(pars{ipar})(isub,:,:,:,1) - behavior.(model_name).(pars{ipar})(isub,:,:,:,2);
      end
    end
  end
  if strcmp(model_name, 'chi_prevresp_z_dc_runs')
    behavior.(model_name).histshift_z = behavior.(model_name).z(:,:,:,:,2) - behavior.(model_name).z(:,:,:,:,1); % histshift = prevR-prevL
    behavior.(model_name).histshift_dc = behavior.(model_name).dc(:,:,:,:,2) - behavior.(model_name).dc(:,:,:,:,1); % histshift = prevR-prevL
  end
end



%%
% disp 'load ddm fits basic and prevresp histbias model runs'
% model_names = {'chi_accuracy_basic_runs' 'chi_prevresp_z_dc_runs'};
% pars = {'a' 't' 'v' 'z' 'dc'};
% drugs = {'atx', 'plac'};
% ses = {'contra', 'ipsi'};
% diffs = {'easy' 'hard'};
% prevresps = {'Lprev' 'Rprev'};
% 
% for imodel = 1:length(model_names)
%   model_name = model_names{imodel};
%   behavior.(model_name) = [];
%   behavior.(model_name).dimord = 'subj_runs_drug_motor_diffodatselrevresp';
%   isub=0;
%   for subind = SUBJ_idx
%     isub=isub+1;
%     ddmfit = readtable(fullfile(ddmpath, sprintf('%s_subj%d.csv', model_name, subind) ));
%     varnames = ddmfit.Properties.VariableNames';
%     for ipar = 1:length(pars)
%       for idrug = 1:2
%         for ises = 1:2
%           cols = startsWith(varnames, pars{ipar}) & contains(varnames, drugs{idrug}) ...
%             & contains(varnames, ses{ises}) & ~contains(varnames, 'trans');
%           if ~any(cols); continue;  end
%           if any(contains(varnames(cols), 'easy')) %  easy vs hard (v)
%             allcols = cols;
%             for idiff=1:2 % put diff dim last
%               temp = NaN(1,8);
%               cols = allcols & contains(varnames, diffs{idiff});
%               temp(1,1:length(find(cols))) = table2array(ddmfit(1,cols));
%               behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises,idiff) = temp;
%             end
%           elseif any(contains(varnames(cols), 'prev')) %  easy vs hard (v)
%             allcols = cols;
%             for iprev=1:2 % put prevresp dim last
%               temp = NaN(1,8);
%               cols = allcols & contains(varnames, prevresps{iprev});
%               temp(1,1:length(find(cols))) = table2array(ddmfit(1,cols));
%               behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises,iprev) = temp;
%             end
%           else % no diff or prevResp
%             temp = NaN(1,8);
%             temp(1,1:length(find(cols))) = table2array(ddmfit(1,cols));
%             behavior.(model_name).(pars{ipar})(isub,1:8,idrug,ises) = temp;
%           end
%         end
%       end
%       if isfield(behavior.(model_name), pars{ipar})
%         behavior.(model_name).(pars{ipar})(isub,9,:,:,:) = nanmean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),2);
%         behavior.(model_name).(pars{ipar})(isub,:,3,:,:) = mean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),3);
%         behavior.(model_name).(pars{ipar})(isub,:,:,3,:) = mean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),4);
%         behavior.(model_name).(pars{ipar})(isub,:,:,:,3) = mean(behavior.(model_name).(pars{ipar})(isub,:,:,:,:),5);
%         behavior.(model_name).(pars{ipar})(isub,:,4,:,:) = behavior.(model_name).(pars{ipar})(isub,:,1,:,:) - behavior.(model_name).(pars{ipar})(isub,:,2,:,:);
%         behavior.(model_name).(pars{ipar})(isub,:,:,4,:) = behavior.(model_name).(pars{ipar})(isub,:,:,1,:) - behavior.(model_name).(pars{ipar})(isub,:,:,2,:);
%         behavior.(model_name).(pars{ipar})(isub,:,:,:,4) = behavior.(model_name).(pars{ipar})(isub,:,:,:,1) - behavior.(model_name).(pars{ipar})(isub,:,:,:,2);
%       end
%     end
%   end
%   if strcmp(model_name, 'chi_prevresp_z_dc_runs')
%     behavior.(model_name).histshift_z = behavior.(model_name).z(:,:,:,:,2) - behavior.(model_name).z(:,:,:,:,1); % histshift = prevR-prevL
%     behavior.(model_name).histshift_dc = behavior.(model_name).dc(:,:,:,:,2) - behavior.(model_name).dc(:,:,:,:,1); % histshift = prevR-prevL
%   end
% end
% 
% %%
% disp 'load ddm fits computed per session'
% model_names = {'chi_accuracy_basic' 'chi_prevresp_z_dc'};
% for imodel = 1:length(model_names)
%   model_name = model_names{imodel};
%   behavior.(model_name) = [];
%   behavior.(model_name).dimord = 'subj_drug_motor_difforprevresp';
%   ddmpath = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/';
%   ddmfit = readtable(fullfile(ddmpath, sprintf('%s.csv', model_name) ));
%   disp 'drop NK2'
%   ddmfit = ddmfit([SUBJ_idx+1],:);
%   varnames = ddmfit.Properties.VariableNames';
%   for ipar = 1:length(pars)
%     for idrug = 1:2
%       for ises = 1:2
%         cols = startsWith(varnames, pars{ipar}) & contains(varnames, drugs{idrug}) ...
%           & contains(varnames, ses{ises}) & ~contains(varnames, 'trans');
%         if ~any(cols); continue;  end
%         if any(contains(varnames(cols), 'easy')) %  easy vs hard (v)
%           allcols = cols;
%           for idiff=1:2 % put diff dim last
%             cols = allcols & contains(varnames, diffs{idiff});
%             behavior.(model_name).(pars{ipar})(:,1,idrug,ises,idiff) = table2array(ddmfit(:,cols));
%           end
%         elseif any(contains(varnames(cols), 'prev')) %  easy vs hard (v)
%           allcols = cols;
%           for iprev=1:2 % put prevresp dim last
%             cols = allcols & contains(varnames, prevresps{iprev});
%             behavior.(model_name).(pars{ipar})(:,1,idrug,ises,iprev) = table2array(ddmfit(:,cols));
%           end
%         else % no diff or prevResp         
%           behavior.(model_name).(pars{ipar})(:,1,idrug,ises) = table2array(ddmfit(:,cols));
%         end
%       end
%     end
%     disp 'contrasts and averages'
%     if isfield(behavior.(model_name), pars{ipar})
%       behavior.(model_name).(pars{ipar})(:,:,3,:,:) = mean(behavior.(model_name).(pars{ipar}),3);
%       behavior.(model_name).(pars{ipar})(:,:,:,3,:) = mean(behavior.(model_name).(pars{ipar}),4);
%       behavior.(model_name).(pars{ipar})(:,:,:,:,3) = mean(behavior.(model_name).(pars{ipar}),5);
%       behavior.(model_name).(pars{ipar})(:,:,4,:,:) = behavior.(model_name).(pars{ipar})(:,:,1,:,:) - behavior.(model_name).(pars{ipar})(:,:,2,:,:);
%       behavior.(model_name).(pars{ipar})(:,:,:,4,:) = behavior.(model_name).(pars{ipar})(:,:,:,1,:) - behavior.(model_name).(pars{ipar})(:,:,:,2,:);
%       behavior.(model_name).(pars{ipar})(:,:,:,:,4) = behavior.(model_name).(pars{ipar})(:,:,:,:,1) - behavior.(model_name).(pars{ipar})(:,:,:,:,2);
%     end
%   end
%   if strcmp(model_name, 'chi_prevresp_z_dc')
%     behavior.(model_name).histshift_z = behavior.(model_name).z(:,:,:,:,2) - behavior.(model_name).z(:,:,:,:,1); % histshift = prevR-prevL
%     behavior.(model_name).histshift_dc = behavior.(model_name).dc(:,:,:,:,2) - behavior.(model_name).dc(:,:,:,:,1); % histshift = prevR-prevL
%   end
% end
% 
% %%
% disp 'load ddm fits only 1 ses: collapsed over motor'
% model_names = { 'chi_prevresp_z_dc_nomotor' 'chi_accuracy_basic_nomotor'}; 
% for imodel = 1:length(model_names)
%   model_name = model_names{imodel};
%   behavior.(model_name) = [];
%   behavior.(model_name).dimord = 'subj_drug_motor_difforprevresp';
%   ddmpath = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/HDDM/';
%   ddmfit = readtable(fullfile(ddmpath, sprintf('%s.csv', model_name) ));
%   disp 'drop NK2'
%   ddmfit = ddmfit([SUBJ_idx+1],:);
%   varnames = ddmfit.Properties.VariableNames';
%   for ipar = 1:length(pars)
%     for idrug = 1:2
%       for ises = 1%:2 only 1 ses:collapsed over motor
%         cols = startsWith(varnames, pars{ipar}) & contains(varnames, drugs{idrug}) ...
%           & ~contains(varnames, 'trans'); % & contains(varnames, ses{ises})
%         if ~any(cols); continue;  end
%         if any(contains(varnames(cols), 'easy')) %  easy vs hard (v)
%           allcols = cols;
%           for idiff=1:2 % put diff dim last
%             cols = allcols & contains(varnames, diffs{idiff});
%             behavior.(model_name).(pars{ipar})(:,1,idrug,ises,idiff) = table2array(ddmfit(:,cols));
%           end
%         elseif any(contains(varnames(cols), 'prev')) %  easy vs hard (v)
%           allcols = cols;
%           for iprev=1:2 % put prevresp dim last
%             cols = allcols & contains(varnames, prevresps{iprev});
%             behavior.(model_name).(pars{ipar})(:,1,idrug,ises,iprev) = table2array(ddmfit(:,cols));
%           end
%         else % no diff or prevResp         
%           behavior.(model_name).(pars{ipar})(:,1,idrug,ises) = table2array(ddmfit(:,cols));
%         end
%       end
%     end
%     disp 'contrasts and averages'
%     if isfield(behavior.(model_name), pars{ipar})
%       behavior.(model_name).(pars{ipar})(:,:,3,:,:) = mean(behavior.(model_name).(pars{ipar}),3);
%       behavior.(model_name).(pars{ipar})(:,:,:,3,:) = mean(behavior.(model_name).(pars{ipar}),4);
%       behavior.(model_name).(pars{ipar})(:,:,:,:,3) = mean(behavior.(model_name).(pars{ipar}),5);
%       behavior.(model_name).(pars{ipar})(:,:,4,:,:) = behavior.(model_name).(pars{ipar})(:,:,1,:,:) - behavior.(model_name).(pars{ipar})(:,:,2,:,:);
%       behavior.(model_name).(pars{ipar})(:,:,:,4,:) = behavior.(model_name).(pars{ipar})(:,:,:,1,:) - behavior.(model_name).(pars{ipar})(:,:,:,2,:);
%       behavior.(model_name).(pars{ipar})(:,:,:,:,4) = behavior.(model_name).(pars{ipar})(:,:,:,:,1) - behavior.(model_name).(pars{ipar})(:,:,:,:,2);
%     end
%   end
%   if strcmp(model_name, 'chi_prevresp_z_dc_nomotor')
%     behavior.(model_name).histshift_z = behavior.(model_name).z(:,:,:,:,2) - behavior.(model_name).z(:,:,:,:,1); % histshift = prevR-prevL
%     behavior.(model_name).histshift_dc = behavior.(model_name).dc(:,:,:,:,2) - behavior.(model_name).dc(:,:,:,:,1); % histshift = prevR-prevL
%   end
% end

%%
behavior.PREOUT = '/Users/kloosterman/Dropbox/PROJECTS/MEG2afc/plots';
disp 'save behavior struct array'
save(fullfile(PREOUT, 'behavstruct.mat'), 'behavior')





%% tryout plot 
% % close all
% if ismac
% %   datsel = squeeze(behavior.p_repeat.balanced(:,9,[1:2,4],3,3))
% %   datsel = squeeze(behavior.p_repeat.unbalanced(:,9,[1:2,4],3))
%   datsel = squeeze(behavior.button_bias(:,9,[1:2,4],3,1))-0.5
%   figure;  bar(nanmean(datsel));
%   
% %   ylim([0.46 0.62])
%   [~,p] = ttest(datsel(:,1)-0.5)
%   [~,p] = ttest(datsel(:,2)-0.5)
%   [~,p] = ttest(datsel(:,1), datsel(:,2)); title(p)
%   
%   hold on
%   plotdata = cell(1,2);
%   plotdata{1,1} = datsel(:,1);
%   plotdata{1,2} = datsel(:,2);
%   h = plotSpread(plotdata);
% end


% %% some plotting
% v = behavior.(model_name).v([1, 3:end],:,:,:,:);
% % v = mean(v, 1); %subj
% v = nanmean(v, 2); %runs
% v = mean(v, 4); % motor
% % v=(v); % drug diff
% figure; b = bar(squeeze(mean(v))); legend(diffs);
% ax=gca; ax.XTickLabel = drugs;
% [~,p] = ttest(squeeze(v(:,:,1,:,1)), squeeze(v(:,:,2,:,1)));
% [h,p] = ttest(squeeze(v(:,:,1,:,2)), squeeze(v(:,:,2,:,2)))

% %% corr runs modeled vs sessions modeled ddms
% dc=[]
% dc.one = behavior.chi_prevresp_z_dc.dc;
% dc.runs = behavior.chi_prevresp_z_dc_runs.dc;
% dc.runs = squeeze(nanmean(dc.runs,2));
% corrdat = [dc.one(:), dc.runs(:)];
% corrdat = corrdat(~isnan(corrdat(:,1)),:);
% corrdat = corrdat(~isnan(corrdat(:,2)),:);
% r = corr(corrdat(:,1), corrdat(:,2));
% figure; scatter(corrdat(:,1), corrdat(:,2)); title(r)
% axis tight; box on

%% TODO corr p_repeat with dc
% % histshift = behavior.chi_prevresp_z_dc_runs.dc(:,:,:,:,1) - behavior.chi_prevresp_z_dc_runs.dc(:,:,:,:,2);
% histshift = behavior.chi_prevresp_z_dc.dc(:,:,:,2) - behavior.chi_prevresp_z_dc.dc(:,:,:,1); % histshift = prevR-prevL
% % histshift = behavior.chi_prevresp_z_dc.z(:,:,:,2) - behavior.chi_prevresp_z_dc.z(:,:,:,1); % histshift = prevR-prevL
% histshift = squeeze(mean(histshift,3));
% histshift = histshift(:,1) - histshift(:,2);
% p_rep = squeeze(nanmean(behavior.p_repeat,2));
% p_rep = squeeze(mean(p_rep,3));
% p_rep = p_rep(:,1) - p_rep(:,2);
% 
% % corrdat = [behavior.p_repeat(:) histshift(:)];
% corrdat = [p_rep(:) histshift(:)];
% corrdat = corrdat(~isnan(corrdat(:,1)),:)
% corrdat = corrdat(~isnan(corrdat(:,2)),:)
% r = corr(corrdat(:,1), corrdat(:,2));
% figure; scatter(corrdat(:,1), corrdat(:,2)); title(r)
% axis tight; box on
% 


