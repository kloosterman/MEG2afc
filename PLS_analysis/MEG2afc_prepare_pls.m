
addpath(genpath('/Users/kloosterman/Dropbox/MATLAB/toolbox/MEGPLS_PIPELINE_v2.02b/MEGPLS_TOOLBOX/PLS_Tools'))

%resplocked:
load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/stats_5D/correlation/powermod/full/Pearson/freq/freq_idrug4_idiff3_itrig2_itest2_idat1.mat')

% stimlocked:
% load('/Users/kloosterman/gridmaster2012/kloosterman/projectdata/MEG2afc/freq/stats_5D/correlation/powermod/full/Pearson/freq/freq_idrug4_idiff3_itrig1_itest2_idat1.mat')

%  Usage: result = pls_analysis(datamat_lst, num_subj_lst, num_cond, ...
%			[option])

% % all the data
% cfg=[];
% cfg.latency = [-0.25 1];
% % cfg.frequency = [2 35];
% cfg.frequency = [2 99];
% data = ft_selectdata(cfg,data)
% 

% make 3 time bins
latencies = [-1.5 -1; -1 -0.5; -0.5 0];
datamat_lst = [];
for itim = 1:3
    cfg=[];
    cfg.latency = latencies(itim,:);
    data_sel = ft_selectdata(cfg,data);
    datamat_lst = [datamat_lst; reshape(data_sel.powspctrm, 18, [])];
end
datamat_lst = {datamat_lst};


% datamat_lst = {reshape(data.powspctrm, 18, [])};


num_subj_lst = 18;
% num_cond = 1;
num_cond = 3;

% behavdata = zeros(size(data.powspctrm));
% for isub=1:18
%     behavdata(isub,:,:,:) = data.drifts(isub);
% end
% stacked_behavdata = reshape(behavdata, [], 1);
stacked_behavdata = data.drifts';

%%
option = [];
option.method = 3; % [1] | 2 | 3 | 4 | 5 | 6
option.num_perm = 500; %( single non-negative integer )
% option.is_struct = [0] | 1
% option.num_split = 100 %( single non-negative integer )
option.num_boot = 100 %500 % ( single non-negative integer )
% option.clim = ( [95] single number between 0 and 100 )
% option.bscan = ( subset of  1:num_cond )
% option.stacked_designdata = ( 2-D numerical matrix )

% option.stacked_behavdata = stacked_behavdata; %( 2-D numerical matrix )
option.stacked_behavdata = repmat(stacked_behavdata, num_cond, 1); %( 2-D numerical matrix )

% option.meancentering_type = [0] | 1 | 2 | 3
option.cormode = 0; % [0] | 2 | 4 | 6
% option.boot_type = ['strat'] | 'nonstrat'

%
result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option)


%% TODO plot LV's in TSF
lvdat = reshape(result.boot_result.compare_u(:,1), 268, 48, 11);
corrdat = reshape(result.datamatcorrs_lst{1}, 268, 48, 36);
udat = reshape(result.u, 268, 48, 36);

% stat = 
%                    prob: [268×48×36 double]       lvdat
%             posclusters: [1×233 struct]           put lv pval in stat.posclusters.prob
%     posclusterslabelmat: [268×48×36 double]       put 1 for lvdat>3 and <3
%         posdistribution: [1×1000 double]
%             negclusters: [1×219 struct]
%     negclusterslabelmat: [268×48×36 double]
%         negdistribution: [1×1000 double]
%                 cirange: [268×48×36 double]
%                    mask: [268×48×36 logical]
%                    stat: [268×48×36 single]
%                     ref: [268×48×36 double]
%                     rho: [268×48×36 single]
%                  dimord: 'chan_freq_time'
%                    freq: [1×48 double]
%                   label: {268×1 cell}
%                    time: [1×36 double]
%                     cfg: [1×1 struct]

stat.prob = lvdat;
stat.rho = lvdat;
stat.posclusters = [];
stat.posclusters.prob = result.perm_result.sprob;
stat.posclusterslabelmat = lvdat > 3 | lvdat < -3;
stat.negclusterslabelmat = zeros(size(stat.posclusterslabelmat));

% allstat{end} = stat;

%
temp = squeeze(sum(stat.posclusterslabelmat));
% temp = squeeze(sum(udat));
% temp = squeeze(sum(corrdat));
close all
figure; 
imagesc(data.time, data.freq, temp)
colorbar
set(gca, 'Ydir', 'Normal')
% caxis([60 140])
caxis([10 60])
