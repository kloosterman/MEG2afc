function MEG2afc_behaviorplot(b)

%% plot modelfree
SAV=1;
condlabels = {'ATX' 'plac'}; % TODO put in b struct
% behavnames = {'dprime' 'criterion' 'button_bias' 'p_repeatbalanced' 'RT' 'RTsd'}; % 5D matrices
behavnames = {
  {'dprime'};  { 'ntrials' }; %{'criterion'};  { 'RT'} ;{ 'RTsd' };
%   {'basepupil' } ;{'bpm'};  { 'p_repeatbalanced'};  {'p_repeatunbalanced' }; {'button_bias'};
  
%   %   {'chi_accuracy_basic_nomotor' 'v'}; {'chi_accuracy_basic_nomotor' 'a'}; {'chi_accuracy_basic_nomotor' 't'};  {};  {};
%   %   {'chi_accuracy_basic' 'v'}; {'chi_accuracy_basic' 'a'}; {'chi_accuracy_basic' 't'}; {}; {};
%   {'ddm_acc_perrun' 'v'}; {'ddm_acc_perrun' 'a'}; {'ddm_acc_perrun' 't'}; {}; {};
%   {'ddm_acc_perrun_ol' 'v'}; {'ddm_acc_perrun_ol' 'a'}; {'ddm_acc_perrun_ol' 't'}; {}; {};
%   %   {'chi_accuracy_basic_runs' 'v'}; {'chi_accuracy_basic_runs' 'a'}; {'chi_accuracy_basic_runs' 't'}; {}; {};
%   
%   %   {'chi_prevresp_z_dc_nomotor' 'v'};{'chi_prevresp_z_dc_nomotor' 'a'}; {'chi_prevresp_z_dc_nomotor' 't'}; {'chi_prevresp_z_dc_nomotor' 'histshift_dc'}; {'chi_prevresp_z_dc_nomotor' 'histshift_z'};
%   %   {'chi_prevresp_z_dc' 'v'};{'chi_prevresp_z_dc' 'a'}; {'chi_prevresp_z_dc' 't'}; {'chi_prevresp_z_dc' 'histshift_dc'}; {'chi_prevresp_z_dc' 'histshift_z'};
%   {'ddm_histbias_perrun' 'v'};{'ddm_histbias_perrun' 'a'}; {'ddm_histbias_perrun' 't'}; {'ddm_histbias_perrun' 'histshift_dc'}; {'ddm_histbias_perrun' 'histshift_z'};
%   {'ddm_histbias_perrun_ol' 'v'};{'ddm_histbias_perrun_ol' 'a'}; {'ddm_histbias_perrun_ol' 't'}; {'ddm_histbias_perrun_ol' 'histshift_dc'}; {'ddm_histbias_perrun_ol' 'histshift_z'};
%   %   {'chi_prevresp_z_dc_runs' 'v'};{'chi_prevresp_z_dc_runs' 'a'}; {'chi_prevresp_z_dc_runs' 't'}; {'chi_prevresp_z_dc_runs' 'histshift_dc'}; {'chi_prevresp_z_dc_runs' 'histshift_z'};
%   %   {'ddm_histbias_perses' 'v'};{'ddm_histbias_perses' 'a'}; {'ddm_histbias_perses' 't'}; {'ddm_histbias_perses' 'histshift_dc'}; {'ddm_histbias_perses' 'histshift_z'};
  }; % all matrices
close all
nrow=8; ncol=6;
diffleg = {'strong', 'weak', ''};

avgtypestr = {'averaged' 'together'}; % runs averaged or taken together
f = figure; iplot=0;
for ia = 1:2
  if strcmp(avgtypestr{ia}, 'averaged') % computed once collapsed over runs (9), or averaged over runs (10)
    avgtype = 10;
  elseif strcmp(avgtypestr{ia}, 'together') % computed once collapsed over runs (9), or averaged over runs (10)
    avgtype = 9;
  end
  
  for idiff=1:3
    Fontsize = 6;
    f.Position =[   680   467   85*ncol   100*nrow];
    for im = 1:length(behavnames)
      iplot=iplot+1;
      if isempty(behavnames{im})
        continue;
      end
      curb = getfield(b, behavnames{im}{:});
      if length(size(curb)) == 5
        data = squeeze(curb(:,avgtype, 1:2, 3, idiff));  % dims: subj runs drug motor diff
      elseif length(size(curb)) == 4
        data = squeeze(curb(:,9, 1:2, 3));  % dims: subj runs drug motor
      end
      %   disp 'Drop NK1 high drift'
      %   data = data(2:end,:);
      subplot(nrow,ncol,iplot); hold on; % axis tight
      plotspreadfig(data, Fontsize, condlabels, b.SUBJ);
      title(sprintf('%s %s\nruns %s', behavnames{im}{1}, diffleg{idiff}, avgtypestr{ia}), 'Fontsize', Fontsize-1)
    end
  end
end
if SAV
  %   saveas(gcf, fullfile(b.PREOUT, sprintf('behavior.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavior_%s_runs%s.png', diffleg{idiff}, avgtypestr{ia} )))
  cd(b.PREOUT)
end

%% plot things per session
SAV=1;
condlabels = {'ATX' 'plac'}; % TODO put in b struct
behavnames = {
  {'dprime'}; {'criterion'};  { 'RT'} ;{ 'RTsd' }; { 'ntrials' };
  {'ddm_acc_perrun' 'v'}; {'ddm_acc_perrun' 'a'}; {'ddm_acc_perrun' 't'};
  {'ddm_histbias_perrun_ol' 'v'};{'ddm_histbias_perrun_ol' 'a'}; {'ddm_histbias_perrun_ol' 't'}; {'ddm_histbias_perrun_ol' 'histshift_dc'}; {'ddm_histbias_perrun_ol' 'histshift_z'};
  {'basepupil' } ;{'bpm'};  { 'p_repeatbalanced'};  {'p_repeatunbalanced' }; {'button_bias'};
  
  }; % all matrices
close all
nrow=9; ncol=6;
f = figure; iplot=0;
Fontsize = 6;
f.Position =[   680   467   85*ncol   100*nrow];
for im = 1:length(behavnames)
  if isempty(behavnames{im})
    continue;
  end
  curb = getfield(b, behavnames{im}{:});
  for ises = 1:2
    iplot=iplot+1;
    if length(size(curb)) == 5
      if strcmp(behavnames{im}{1}, 'button_bias')
        data = squeeze(curb(:,end, 1:2, ises, 1));  % pick left bp
      else
        data = squeeze(curb(:,end, 1:2, ises, 3));  % dims: subj runs drug motor diff
      end
    elseif length(size(curb)) == 4
      data = squeeze(curb(:,end, 1:2, ises));  % dims: subj runs drug motor
    end
    %     disp 'Drop NK1 high drift'
    %     data = data(2:end,:);
    subplot(nrow,ncol,iplot); hold on; % axis tight
    plotspreadfig(data, Fontsize, condlabels, b.SUBJ);
    title(behavnames{im}, 'Fontsize', Fontsize-1)
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavior_perses.pdf' )))
  cd(b.PREOUT)
end

%% plot bpm and pupil over runs to see how it progresses
SAV=1;
makezscore = false;
condlabels = {'ATX' 'plac'}; % TODO put in b struct
condcol = {'r' 'b'};
behavnames = {
  {'basepupil' } ;{'bpm'}; {'ntrials'};
      {'dprime'};  
% { 'RT'} ;
  
  }; % all matrices

close all
nrow=4; ncol=2;
f = figure; iplot=0;
Fontsize = 6;
f.Position =[   680   467   85*ncol   100*nrow];
for im = 1:length(behavnames)
  if isempty(behavnames{im})
    continue;
  end
  curb = getfield(b, behavnames{im}{:});
  for ises = 3
    iplot=iplot+1;
    if length(size(curb)) == 5
      if strcmp(behavnames{im}{1}, 'button_bias')
        data = squeeze(curb(:,1:6, 1:2, ises, 1));  % pick left bp
      else
        data = squeeze(curb(:,1:6, 1:2, ises, 3));  % dims: subj runs drug motor diff
      end
    elseif length(size(curb)) == 4
      data = squeeze(curb(:,1:6, 1:2, ises));  % dims: subj runs drug motor
    end
    if makezscore
      data = data(:,:);
      data = (data - nanmean(data,2)) ./ nanstd(data,0,2);
      data = reshape(data, size(data,1), 6, 2);
    end
    subplot(nrow,ncol,iplot); hold on; % axis tight
    
    clear h
    for idrug=1:2
      %       plot(data(:,:,idrug)') %, Fontsize, condlabels, b.SUBJ);
      %       plot(squeeze(nanmean(data))) %, 'k', 'Linewidth', 2) %, Fontsize, condlabels, b.SUBJ);
      h(idrug) = shadedErrorBar([], squeeze(nanmean(data(:,:,idrug))), squeeze(nanstd(data(:,:,idrug))/sqrt(size(data(:,:,idrug),1))), condcol{idrug}, 1 )
    end
    title(behavnames{im}, 'Fontsize', Fontsize-1)
    xlabel('Run no.')
    ax=gca;
    ax.FontSize = Fontsize;
    if im == 2
      legend([h.mainLine], condlabels); legend boxoff
    end
  end
end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavior_overruns.pdf' )))
  cd(b.PREOUT)
end


%% plot RT histograms
SAV=1;
condlabels = {'ATX' 'plac' 'drugavg' 'ATX-plac'}; % TODO put in b struct
condcol = {'r' 'b'};
behavnames = {
  {'RThist' };
  }; % all matrices

close all
nrow=4; ncol=2;
f = figure; iplot=0;
Fontsize = 6;
f.Position =[   680   467   85*ncol   100*nrow];
for im = 1:length(behavnames)
  if isempty(behavnames{im})
    continue;
  end
  curb = getfield(b, behavnames{im}{:});
  for ises = 3
    iplot=iplot+1;
    if length(size(curb)) == 5
      if strcmp(behavnames{im}{1}, 'button_bias')
        data = squeeze(curb(:,1:6, 1:2, ises, 1));  % pick left bp
      else
        data = squeeze(curb(:,1:6, 1:2, ises, 3));  % dims: subj runs drug motor diff
      end
    elseif length(size(curb)) == 4
      data = squeeze(curb(:, 1:2, ises, :));  % dims: subj runs drug motor
    end
    subplot(nrow,ncol,iplot); hold on; % axis tight
    
    clear h m
    for idrug=1:2
      %       plot(data(:,:,idrug)') %, Fontsize, condlabels, b.SUBJ);
      %       plot(squeeze(nanmean(data))) %, 'k', 'Linewidth', 2) %, Fontsize, condlabels, b.SUBJ);
      h(idrug) = shadedErrorBar(b.RThistedges(1:end-1), squeeze(nanmean(data(:,idrug,:))), squeeze(nanstd(data(:,idrug,:))/sqrt(size(data(:,idrug,:),1))), condcol{idrug} )
      [m(idrug) mi(idrug)] = max(squeeze(nanmean(data(:,idrug,:)))); % find mode: 0.65 s
    end
    %     diff(squeeze(nanmean(data(:,idrug,:)))) % at bin 7 starts to
    %     rise: 0.3000 poststim
    title(behavnames{im}, 'Fontsize', Fontsize-1)
    xlabel('RT (s)')
    ylabel('Probability')
    ax=gca;
    ax.FontSize = Fontsize;
    if im == 1
      legend([h.mainLine], condlabels); legend boxoff
    end
  end
end
xlim([0 2.5])
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavior_RThist.pdf' )))
  cd(b.PREOUT)
end



%% correlate prepeat and drift bias
corrnames = {
  { 'p_repeatbalanced'} {'ddm_histbias_perrun_ol' 'histshift_dc'}
  {'p_repeatunbalanced'} {'ddm_histbias_perrun_ol' 'histshift_dc'}
  { 'p_repeatbalanced'} {'ddm_histbias_perrun' 'histshift_dc'}
  {'p_repeatunbalanced'} {'ddm_histbias_perrun' 'histshift_dc'}
  { 'p_repeatbalanced'} {'ddm_histbias_perrun' 'histshift_z'}
  {'p_repeatunbalanced'} {'ddm_histbias_perrun' 'histshift_z'}
  {'ddm_histbias_perrun' 'histshift_dc'} {'ddm_histbias_perrun' 'histshift_z'}
  {'dprime'} {'ddm_acc_perrun' 'v'}
  };

nrow=4; ncol=2;

for idrug = [1,2,4]
  f = figure; iplot=0;
  Fontsize = 6;
  f.Position =[   680   467   125*ncol   125*nrow];
  for im = 1:length(corrnames)
    iplot=iplot+1;
    if isempty(corrnames{im})
      continue;
    end
    acorr = getfield(b, corrnames{im,1}{:});
    if length(size(acorr)) == 5
      acorr = squeeze(acorr(:,end, idrug, 3, 3));  % dims: subj runs drug motor diff
    elseif length(size(acorr)) == 4
      acorr = squeeze(acorr(:,end, idrug, 3));  % dims: subj runs drug motor
    end
    bcorr = getfield(b, corrnames{im,2}{:});
    if length(size(bcorr)) == 5
      bcorr = squeeze(bcorr(:,end, idrug, 3, 3));  % dims: subj runs drug motor diff
    elseif length(size(bcorr)) == 4
      bcorr = squeeze(bcorr(:,end, idrug, 3));  % dims: subj runs drug motor
    end
    %   disp 'Drop NK1 high drift'
    %   acorr = acorr([1:4, 6:end],:);
    %   bcorr = bcorr([1:4, 6:end],:);
    %   bcorr = bcorr(2:end,:);
    subplot(nrow,ncol,iplot); hold on; % axis tight
    
    scatter(acorr, bcorr, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 30);
    box on; axis square;
    ax=gca;
    ax.FontSize = Fontsize;
    %   text(acorr, bcorr, num2cell(b.SUBJ_idx), 'Fontsize', 7);
    [rho, p] = corr(acorr, bcorr);
    if p < 0.05
      lsline
    end
    title(sprintf('r=%1.4f\np=%1.4f', rho, p))
    xlabel(corrnames{im,1})
    ylabel(corrnames{im,2})
    axis tight
    %   title(behavnames{im}, 'Fontsize', Fontsize-1)
  end
  if SAV
    saveas(gcf, fullfile(b.PREOUT, sprintf('behavcorr_%s.pdf', condlabels{idrug} )))
    cd(b.PREOUT)
  end
end
  
  %% correlate drift, bound and NDT among each other for atx, plac, atx-plac
  corrnames = {
    {'ddm_acc_perrun' 'v'} {'ddm_acc_perrun' 'a'}
    {'ddm_acc_perrun' 'v'} {'ddm_acc_perrun' 't'}
    {'ddm_acc_perrun' 'a'} {'ddm_acc_perrun' 't'}
    %   {'ddm_acc_perrun_ol' 'v'} {'ddm_acc_perrun_ol' 'a'}
    %   {'ddm_acc_perrun_ol' 'v'} {'ddm_acc_perrun_ol' 't'}
    %   {'ddm_acc_perrun_ol' 'a'} {'ddm_acc_perrun_ol' 't'}
    
    %   { 'p_repeatbalanced'} {'ddm_histbias_perrun_ol' 'histshift_dc'}
    %   {'p_repeatunbalanced'} {'ddm_histbias_perrun_ol' 'histshift_dc'}
    %   { 'p_repeatbalanced'} {'ddm_histbias_perrun' 'histshift_dc'}
    %   {'p_repeatunbalanced'} {'ddm_histbias_perrun' 'histshift_dc'}
    %   { 'p_repeatbalanced'} {'ddm_histbias_perrun' 'histshift_z'}
    %   {'p_repeatunbalanced'} {'ddm_histbias_perrun' 'histshift_z'}
    %   {'ddm_histbias_perrun' 'histshift_dc'} {'ddm_histbias_perrun' 'histshift_z'}
    %   {'dprime'}
    };
  drugleg = {'ATX', 'plac', 'drugavg', 'ATX-plac'};
  
  nrow=4; ncol=3;
  
  f = figure; iplot=0;
  Fontsize = 6;
  f.Position =[   680   467   125*ncol   125*nrow];
  for idrug = [1,2,4]
    for im = 1:length(corrnames)
      iplot=iplot+1;
      if isempty(corrnames{im})
        continue;
      end
      acorr = getfield(b, corrnames{im,1}{:});
      if length(size(acorr)) == 5
        acorr = squeeze(acorr(:,end, idrug, 3, 3));  % dims: subj runs drug motor diff
      elseif length(size(acorr)) == 4
        acorr = squeeze(acorr(:,end, idrug, 3));  % dims: subj runs drug motor
      end
      bcorr = getfield(b, corrnames{im,2}{:});
      if length(size(bcorr)) == 5
        bcorr = squeeze(bcorr(:,end, idrug, 3, 3));  % dims: subj runs drug motor diff
      elseif length(size(bcorr)) == 4
        bcorr = squeeze(bcorr(:,end, idrug, 3));  % dims: subj runs drug motor
      end
      %   disp 'Drop NK1 high drift'
      %   acorr = acorr([1:4, 6:end],:);
      %   bcorr = bcorr([1:4, 6:end],:);
      %   bcorr = bcorr(2:end,:);
      subplot(nrow,ncol,iplot); hold on; % axis tight
      
      scatter(acorr, bcorr, 'filled', 'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'w', 'sizedata', 30);
      box on; axis square
      ax=gca;
      ax.FontSize = Fontsize;
      %     text(acorr, bcorr, num2cell(b.SUBJ_idx), 'Fontsize', 7);
      [rho, p] = corr(acorr, bcorr);
      if p < 0.05
        lsline
      end
      title(sprintf('%s\nrho=%1.2f, p=%1.3f', drugleg{idrug}, rho, p))
      %     xlabel(corrnames{im,1})
      %     ylabel(corrnames{im,2})
      xlabel(corrnames{im,1}{2})
      ylabel(corrnames{im,2}{2})
      %   title(behavnames{im}, 'Fontsize', Fontsize-1)
      axis tight
    end
  end
if SAV
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavcorr_DDMavt.pdf' )))
  saveas(gcf, fullfile(b.PREOUT, sprintf('behavcorr_DDMavt.png' )))
  cd(b.PREOUT)
end

end % function




%% plot bars + lines
function ax = plotspreadfig(data, Fontsize, condlabels, SUBJ)
% lines between pts
line([ones(1,size(data,1)); ones(1,size(data,1))+1], [data(:,1) data(:,2)]', 'Color',[0.66 0.66 0.66 ], 'Linewidth', 0.5)

handle = plotSpread( {data(:,1) data(:,2)}, 'distributionMarkers', 'o', 'distributionColors', [1 0.5 0.5; 0.5 0.5 1] );% b r circles
set(handle{1}(1:2), 'MarkerSize', 3, 'Linewidth', 0.5)
% text(data(:,1)', data(:,2)', num2cell(SUBJ), 'Fontsize', Fontsize)

ax=gca;
ax.FontSize = Fontsize;
ax.XTickLabel = condlabels;

% mean lines
line([0.75 1.25]', [mean(data(:,1)) mean(data(:,1))]',  'Color', 'r', 'Linewidth', 2)
line([1.75 2.25]', [nanmean(data(:,2)) nanmean(data(:,2))]',  'Color', 'b', 'Linewidth', 2)

% stats
[~,p] = ttest(data(:,1), data(:,2));
text(1, ax.YLim(2), sprintf('p=%1.3f', p), 'Fontsize', Fontsize)
%   if p < 0.05 %   sigstar
xlim([0.5 2.5])
end

