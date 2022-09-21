function MEG2afc_rmcorr_allbehavphys(b)
% rm corr all behav_names, plot mosaic

condlabels = {'ATX' 'plac' '' 'ATX-plac'}; %

behavnames = {  {'dprime'}; {'criterion'};  {'RT'};  {'RTsd'}; { 'p_repeatbalanced'}; {'basepupil' } ;  {'bpm'}; };
behavnames = {  {'dprime'};  {'RT'};  {'basepupil' } ;  {'bpm'}; };

% behavnames = fliplr(behavnames)

% get in matrix, columns
close all
corrdat = NaN(228,length(behavnames),4);
for idrug = [1,2,4]
  for im = 1:length(behavnames)
    if isempty(behavnames{im})
      continue;
    end
    curb = getfield(b, behavnames{im}{:});
    
    data_demean = [];
    for ises = 1:2
      if length(size(curb)) == 5
        if strcmp(behavnames{im}{1}, 'button_bias')
          data = squeeze(curb(:,1:6, idrug, ises, 1));  % pick left bp
        else
          data = squeeze(curb(:,1:6, idrug, ises, 3));  % dims: subj runs drug motor diff
        end
      elseif length(size(curb)) == 4
        data = squeeze(curb(:,1:6, idrug, ises));  % dims: subj runs drug motor
      end
      data_demean(:,:,ises) = data - nanmean(data,2);
    end
    data_demean = data_demean(:);
    disp 'drop outliers'
    outliers = 1;
    while outliers
      outliers = (data_demean - nanmean(data_demean)) / nanstd(data_demean);
      outliers = find(abs(outliers) > 3);
      data_demean(outliers) = NaN;
    end
    
    corrdat(:,im,idrug) = data_demean;
  end
  clean_inds = all(~isnan(corrdat(:,:,idrug)),2);
  [corrmat(:,:,idrug), corrmatp(:,:,idrug)] = corr(corrdat(clean_inds,:,idrug));
  
  f=figure;
  f.Position =[        1000         473         964         865];
  cr = corrplot(corrdat(clean_inds,:,idrug));
  
end

%% plot corr matrix
% close all
f=figure;
f.Position =[        1000         473         964         865];
load colormap_jetlightgray.mat
for idrug = [1,2,4]
  subplot(2,2,idrug)
%   im=imagesc(triu(corrmat(:,:,idrug), 0), [-0.3 0.3]);
  im=imagesc(corrmat(:,:,idrug), [-0.3 0.3]);
  alphadat = double(corrmatp(:,:,idrug) < 0.05);
  alphadat(alphadat==0) = 0.2;
  im.AlphaData = triu(alphadat);
%   im.AlphaData = alphadat;
  %   imagesc(corrmat(:,:,idrug), [-0.3 0.3])
  colormap(cmap)
  %   colormap(jet(256))
  colorbar
  ax=gca;
  ax.XTick = 1:length([behavnames{:}]);
  ax.YTick = 1:length([behavnames{:}]);
  ax.XTickLabel = [behavnames{:}];
  ax.YTickLabel = [behavnames{:}];
  ax.XTickLabelRotation = 45;
  %   ax.YDir = 'normal';
  %   ax.XDir = 'reverse';
  axis square
  hold on
  %     contour(corrmatp(:,:,idrug) < 0.05,1)
  title(condlabels{idrug})
  
end
saveas(f, fullfile(b.PREOUT, 'corrmat_behav.pdf'))