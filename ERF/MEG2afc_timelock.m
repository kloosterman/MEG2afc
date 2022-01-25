function megdescriptives = MEG2afc_timelock(cfg)
% check resplocked freq, better do fixed effects with appenddata?, no good for
% removing ERP per run
PREIN = cfg.PREIN;
SUBJ = cfg.SUBJ;
sesname = cfg.sesname;
outfile = cfg.outfile;
cd(PREIN)

runlist = dir(sprintf('%s_run*.mat', sesname)); % put simon together for now

trls = [];
timelock_all = {};
clear megdescriptives
for irun = 1:length(runlist)
  data = {};
  disp('Loading');    disp(runlist(irun).name)
  load(runlist(irun).name);
  
  trl = ft_findcfg(data_eye.cfg, 'trl');
  disp 'drop trials w RT <0.3 s'
  cfg=[];
  cfg.trials = data.trialinfo(:,5) > (0.3*1200);
  data = ft_selectdata(cfg, data);
  
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
  
%   disp(' synthetic planar computation')
%   cfg              = [];
%   cfg.feedback     = 'no';
%   cfg.method       = 'template';
%   if ismac
%     cfg.template     = 'ctf275_neighb.mat';
%   else
%     cfg.template     = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
%   end
%   cfg.planarmethod = 'sincos';
%   cfg.channel      = {'MEG'};
%   cfg.trials       = 'all';
%   cfg.neighbours   = ft_prepare_neighbours(cfg, data);
%   data =  ft_megplanar(cfg, data);
  
%   disp('compute resplocked')
%   temp = data;    data = {};
%   data{1} = temp;   clear temp
%   cfg=[];
%   cfg.offset = -round((data{1}.trialinfo(:,5) / ft_findcfg(data{1}.cfg, 'origfs')) * ft_findcfg(data{1}.cfg, 'resamplefs'));
%   data{2} = ft_redefinetrial(cfg,data{1});
  
  disp('timelockanalysis')
  
  disp(' downsample')
  cfg=[];
  cfg.resample = 'yes';
  cfg.resamplefs = 100;
  data = ft_resampledata(cfg, data);
  
  disp('discard trials with RT < 0.2 and > 2.5 and crop a bit')
  cfg=[];
  cfg.trials = data.trialinfo(:,5)/1200 > 0.2 & data.trialinfo(:,5)/1200 < 2.5;
  data = ft_selectdata(cfg, data);
  
  cfg=[];
  cfg.toilim = [-0.5 0];
  data = ft_redefinetrial(cfg,data);

  disp(' make timelock')
  cfg=[];
  cfg.keeptrials = 'yes';
  timelock = ft_timelockanalysis(cfg,data);  
  megdescriptives(irun).mean =  mean(timelock.trial(:));
  megdescriptives(irun).std =  std(timelock.trial(:));

%   timelock_all{1}{irun} = timelock;
  
end
% disp('average over runs')
% timelock_avg = cellfun(@(x) ft_timelockgrandaverage([], x{:}), timelock_all); % makes struct array
% [timelock_avg(2:end).cfg] = deal([]);  % cfg gets huge after appending data etc. only keep first
% % [timelock_avg.trialinfo] = deal(trls(:,4:end));
% 
% disp('keep runs as rpt')
% cfg = [];
% cfg.keepindividual = 'yes';
% timelock = cellfun(@(x) ft_timelockgrandaverage(cfg, x{:}), timelock_all); % makes struct array
% [timelock.dimord] = deal('rpt_chan_freq_time');  % cfg gets huge after appending data etc
% [timelock(2:end).cfg] = deal([]);  % cfg gets huge after appending data etc. only keep first
% % % [freq.trialinfo] = deal(rep_probs);
% % [freq.trialinfo] = deal(trls(:,4:end));
% 
% fprintf('Saving %s\n', outfile)
% save(outfile, 'freq', 'freq_avg');
% 
% 
% 
% %% %OLD:
% % function MEG2afc_timelockanalysis(subjdir, PREOUT)
% % concatenate preproc data per session
% % compute ERF with timelockanalysis
% % save result
% %%%set paths
% 
% %' sesdirs = { 'drug_ipsi' 'drug_contra' 'plac_ipsi' 'plac_contra'}; % Drug_ipsi, Drug_contra, plac_ipsi, plac_contra;  todo keep track of order per subj
% if nargin == 0
%   %   subjdir = '/Users/kloosterman/beegfs/projectdata/MEG2afc/preproc/NK3';
%   subjdir = '/Users/kloosterman/projectdata/MEG2afc/preproc/NK3';
%   PREOUT = '/';
% end
% 
% disp(subjdir)
% [~,SUBJ] = fileparts(subjdir);
% 
% drugleg = {'drug' 'plac'};
% motorleg = {'ipsi' 'contra'};
% 
% %%
% disp(' Concatenate preproc data')
% timelockout = {};
% for idrug = 1:2
%   concatdata = {};
%   for imotor = 1:2
%     curdir = fullfile(subjdir, sprintf('%s_%s', drugleg{idrug}, motorleg{imotor}));
%     try cd(curdir)
%     catch; fprintf('%s not found\n', curdir)
%       continue
%     end
%     
%     runlist = dir('*data.mat');
%     for irun = 1%:length(runlist)
%       fprintf('Loading %s run %d . . .\n',  curdir, irun)
%       load(runlist(irun).name)
%       
%       data.trialinfo(:,9) = idrug;
%       data.trialinfo(:,10) = imotor;
%       
%       disp(' downsample')
%       cfg=[];
%       cfg.resample = 'yes';
%       cfg.resamplefs = 100;
%       data = ft_resampledata(cfg, data);
%       
%       disp('discard trials with RT < 0.2 and > 2.5 and crop a bit')
%       cfg=[];
%       cfg.trials = data.trialinfo(:,5)/1200 > 0.2 & data.trialinfo(:,5)/1200 < 2.5;
%       data = ft_selectdata(cfg, data);
%       
%       cfg=[];
%       cfg.toilim = [-0.5 3];
%       data = ft_redefinetrial(cfg,data);
%       
%       close all
%       %%
%       disp('realigning MEG')
%       cfg=[];
%       if ismac
%         cfg.template       = {'/Users/kloosterman/Dropbox/tardis_code/MATLAB/MEG_HH_analysis/megrealign/ctf275.mat'};
%         cfg.headmodel   = fullfile('/Volumes/FB-LIP/user/Niels/MEG2afc/MRI/NKdet', SUBJ, [SUBJ '_hdm.mat']);
%       else
%         cfg.template       = {'/home/mpib/kloosterman/MATLAB/MEG_HH_analysis/megrealign/ctf275.mat'};
%         cfg.headmodel   = fullfile('/home/beegfs/kloosterman/projectdata/MEG2afc/MRI', SUBJ, [SUBJ '_hdm.mat']);
%       end
%       cfg.inwardshift = 1;
%       cfg.feedback = 'no';
%       data = ft_megrealign(cfg, data);
%       
%       concatdata{end+1}     = data;
%       clear data
%     end
%   end
%   
%   cfg=[];
%   cfg.appenddim = 'rpt';
%   data = ft_appenddata(cfg, concatdata{:});
%   if iscell(data.cfg.previous)
%     data.cfg.previous = data.cfg.previous{1}; % avoid large cfg struct
%   end
%   clear concatdata
%   
%   baseline = {};
%   baseline2 = {};
%   for itrig = 1:2
%     if itrig == 2
%       cfg=[];
%       cfg.offset = -round((data.trialinfo(:,5) / 1200) * ft_findcfg(data.cfg, 'resamplefs'));
%       data = ft_redefinetrial(cfg,data);
%       cfg=[];
%       cfg.toilim = [-1 0.5];
%       data = ft_redefinetrial(cfg, data);
%     end
%     
%     for idiff = 1:3
%       %%
%       disp(' make timelock')
%       cfg=[];
%       cfg.vartrllength = 2;
%       cfg.keeptrials = 'yes';
%       cfg.trials = data.trialinfo(:,1) == idiff;
%       if idiff == 3
%         cfg.trials = 'all';
%       end
%       timelock = ft_timelockanalysis(cfg,data);
%       
%       if itrig == 1
%         cfg=[];
%         cfg.latency      = [-0.2 0];
%         cfg.avgovertime  = 'yes';
%         baseline{idiff} = ft_selectdata(cfg, timelock);
%       end
%       disp(' baseline correct single trials')
%       timelock.trial = timelock.trial - baseline{idiff}.trial;
%       
%       disp('compute modulation')
%       disp(' average trials')
%       avg = ft_timelockanalysis([], timelock);
%       %%
%       disp(' synthetic planar computation')
%       cfg              = [];
%       cfg.feedback     = 'no';
%       cfg.method       = 'template';
%       if ismac
%         cfg.template     = 'ctf275_neighb.mat';
%       else
%         cfg.template     = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
%       end
%       cfg.planarmethod = 'sincos';
%       cfg.channel      = {'MEG'};
%       cfg.trials       = 'all';
%       cfg.neighbours   = ft_prepare_neighbours(cfg, avg);
%       avg =  ft_megplanar(cfg, avg);
%       avg =  ft_combineplanar([], avg);
%       
%       %%
%       disp('  blc avg again after combineplanar')
%       if itrig == 1
%         cfg=[];
%         cfg.latency      = [-0.2 0];
%         cfg.avgovertime  = 'yes';
%         cfg.parameter    = 'avg';
%         baseline2{idiff}  = ft_selectdata(cfg, avg);
%       end
%       avg.avg = avg.avg - baseline2{idiff}.avg;
%       
%       
%       %%
%       disp('crop stimlocked data once more')
%       if itrig == 1
%         cfg = [];
%         cfg.latency = [-0.5 1];
%         avg = ft_selectdata(cfg, avg);
%       end
%       
%       %%
%       timelockout{1, itrig, idrug, idiff} = avg; % modulation in 1, latr in 2-3
%       
%       %%
%       disp('compute lateralization')
%       for iltr = 2:3
%         disp( [ 'lateralization no.' num2str(iltr)] )
%         
%         cfg=[];
%         cfg.trials = timelock.trialinfo(:,iltr) == 1; % 2=stim, 3=button
%         avg_L = ft_timelockanalysis(cfg, timelock);
%         cfg.trials = timelock.trialinfo(:,iltr) == 2;
%         avg_R = ft_timelockanalysis(cfg, timelock);
%         
%         disp(' synthetic planar computation')
%         cfg              = [];
%         cfg.feedback     = 'no';
%         cfg.method       = 'template';
%         if ismac
%           cfg.template     = 'ctf275_neighb.mat';
%         else
%           cfg.template     = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat';
%         end
%         cfg.planarmethod = 'sincos';
%         cfg.channel      = {'MEG'};
%         cfg.trials       = 'all';
%         cfg.neighbours   = ft_prepare_neighbours(cfg, avg_L);
%         avg_L =  ft_megplanar(cfg, avg_L);
%         avg_R =  ft_megplanar(cfg, avg_R);
%         avg_L =  ft_combineplanar([], avg_L);
%         avg_R =  ft_combineplanar([], avg_R);
%         
%         %%
%         disp('  blc avg again after combineplanar')
%         if itrig == 1
%           cfg=[];
%           cfg.latency      = [-0.2 0];
%           cfg.avgovertime  = 'yes';
%           cfg.parameter    = 'avg';
%           baseline2_L{idiff}  = ft_selectdata(cfg, avg_L);
%           baseline2_R{idiff}  = ft_selectdata(cfg, avg_R);
%         end
%         avg_L.avg = avg_L.avg - baseline2_L{idiff}.avg;
%         avg_R.avg = avg_R.avg - baseline2_R{idiff}.avg;
%         
%         %%
%         disp('compute lateralization')
%         cfg=[];
%         cfg.channelcmb = {'MLC11','MRC11';'MLC12','MRC12';'MLC13','MRC13';'MLC14','MRC14';'MLC15','MRC15';'MLC16','MRC16';'MLC17','MRC17';'MLC21','MRC21';'MLC22','MRC22';'MLC23','MRC23';'MLC24','MRC24';'MLC25','MRC25';'MLC31','MRC31';'MLC32','MRC32';'MLC41','MRC41';'MLC42','MRC42';'MLC51','MRC51';'MLC52','MRC52';'MLC53','MRC53';'MLC54','MRC54';'MLC55','MRC55';'MLC61','MRC61';'MLC62','MRC62';'MLC63','MRC63';'MLF11','MRF11';'MLF12','MRF12';'MLF13','MRF13';'MLF14','MRF14';'MLF21','MRF21';'MLF22','MRF22';'MLF23','MRF23';'MLF24','MRF24';'MLF25','MRF25';'MLF31','MRF31';'MLF32','MRF32';'MLF33','MRF33';'MLF34','MRF34';'MLF35','MRF35';'MLF41','MRF41';'MLF42','MRF42';'MLF43','MRF43';'MLF44','MRF44';'MLF45','MRF45';'MLF46','MRF46';'MLF51','MRF51';'MLF52','MRF52';'MLF53','MRF53';'MLF54','MRF54';'MLF55','MRF55';'MLF56','MRF56';'MLF61','MRF61';'MLF62','MRF62';'MLF63','MRF63';'MLF64','MRF64';'MLF65','MRF65';'MLF66','MRF66';'MLF67','MRF67';'MLO11','MRO11';'MLO12','MRO12';'MLO13','MRO13';'MLO14','MRO14';'MLO21','MRO21';'MLO22','MRO22';'MLO23','MRO23';'MLO24','MRO24';'MLO31','MRO31';'MLO32','MRO32';'MLO33','MRO33';'MLO34','MRO34';'MLO41','MRO41';'MLO42','MRO42';'MLO43','MRO43';'MLO44','MRO44';'MLO51','MRO51';'MLO52','MRO52';'MLO53','MRO53';'MLP11','MRP11';'MLP12','MRP12';'MLP21','MRP21';'MLP22','MRP22';'MLP23','MRP23';'MLP31','MRP31';'MLP32','MRP32';'MLP33','MRP33';'MLP34','MRP34';'MLP35','MRP35';'MLP41','MRP41';'MLP42','MRP42';'MLP43','MRP43';'MLP44','MRP44';'MLP45','MRP45';'MLP51','MRP51';'MLP52','MRP52';'MLP53','MRP53';'MLP54','MRP54';'MLP55','MRP55';'MLP56','MRP56';'MLP57','MRP57';'MLT11','MRT11';'MLT12','MRT12';'MLT13','MRT13';'MLT14','MRT14';'MLT15','MRT15';'MLT16','MRT16';'MLT21','MRT21';'MLT22','MRT22';'MLT23','MRT23';'MLT24','MRT24';'MLT25','MRT25';'MLT26','MRT26';'MLT27','MRT27';'MLT31','MRT31';'MLT32','MRT32';'MLT33','MRT33';'MLT34','MRT34';'MLT35','MRT35';'MLT36','MRT36';'MLT37','MRT37';'MLT41','MRT41';'MLT42','MRT42';'MLT43','MRT43';'MLT44','MRT44';'MLT45','MRT45';'MLT46','MRT46';'MLT47','MRT47';'MLT51','MRT51';'MLT52','MRT52';'MLT53','MRT53';'MLT54','MRT54';'MLT55','MRT55';'MLT56','MRT56';'MLT57','MRT57'};
%         avg = ft_lateralizedpotential(cfg, avg_L, avg_R);
%         avg.label = avg.plotlabel'; % overwrite cmb label
%         
%         disp('crop stimlocked data once more')
%         if itrig == 1
%           cfg = [];
%           cfg.latency = [-0.5 1];
%           avg = ft_selectdata(cfg, avg);
%         end
%         
%         %%
%         timelockout{iltr, itrig, idrug, idiff} = avg; % modulation in 1, latr in 2-3
%         clear avg
%       end
%       
%     end
%   end
%   clear data
% end
% 
% disp(' drug avg and contrast for each trig and difficulty')
% for iltr = 1:3
%   for itrig = 1:2
%     for ictr = 1:3
%       cfg=[];
%       cfg.parameter = 'avg';
%       cfg.operation = '1/2*(x1+x2)';
%       timelockout{iltr, itrig, 3, ictr} = ft_math(cfg, timelockout{iltr, itrig, 1, ictr}, timelockout{iltr, itrig, 2, ictr} );
%       
%       cfg=[];
%       cfg.parameter = 'avg';
%       cfg.operation = 'subtract'; %x1-x2
%       timelockout{iltr, itrig, 4, ictr} = ft_math(cfg, timelockout{iltr, itrig, 1, ictr}, timelockout{iltr, itrig, 2, ictr} );
%       disp('difficulty contrast')
%       timelockout{iltr, itrig, ictr, 4} = ft_math(cfg, timelockout{iltr, itrig, ictr, 1}, timelockout{iltr, itrig, ictr, 2} );
%     end
%   end
% end
% 
% outpath = fullfile(PREOUT, [SUBJ '_timelock.mat']);
% disp('Saving to:')
% disp(outpath)
% save(outpath, 'timelockout');
% 
% if ismac
%   close all
%   cfg=[];
%   cfg.layout = 'CTF275.lay';
%   cfg.colorbar = 'yes';
%   cfg.zlim = 'maxabs';  % zeromax
%   figure
%   %             ft_multiplotER(cfg, timelockout{1, 1, 3}, timelockout{1, 2, 3})
%   %             legend(drugleg);
%   
%   ft_multiplotER(cfg, timelockout{1, 4, 3})
%   legend({'drug - plac'});
%   
%   ft_multiplotER(cfg, timelockout{1, 2, 1}, timelockout{1, 2, 2})
%   legend({'easy' 'hard'})
%   %             figure
%   %             ft_multiplotER(cfg, timelockout{1, 4, 3})
% end
% 
% if ismac
%   close all
%   cfg=[];
%   cfg.layout = 'CTF275.lay';
%   cfg.colorbar = 'yes';
%   cfg.zlim = 'zeromax'; %cfg.xlim = [-0.5 1.5]; cfg.xlim = [-1 0.5];
%   figure
%   %           ft_multiplotER(cfg, avg)
%   iltr = 3, itrig=2, idrug=3, idiff=3
%   ft_multiplotER(cfg, timelockout{iltr, itrig, idrug, idiff});
% end
% 
% 
% 
% %%
% %         disp('Interpolate missing channel')
% %         cfg=[];
% %         cfg.method = 'spline';
% %         if ismac
% %           load ctf275_neighb.mat
% %           cfg.layout = 'CTF275.lay';
% %         else
% %           load('/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/neighbours/ctf275_neighb.mat');
% %           cfg.layout = '/home/mpib/kloosterman/MATLAB/tools/fieldtrip-20170611/template/layout/CTF275.lay';
% %         end
% %         cfg.neighbours = neighbours;
% %         labelcomplete = {neighbours.label}';
% %         cfg.missingchannel = labelcomplete(~ismember(labelcomplete, avg_L.label));
% %         avg_L = ft_channelrepair(cfg, avg_L);
% %         avg_R = ft_channelrepair(cfg, avg_R);
