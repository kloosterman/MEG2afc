
cd /Users/kloosterman/gridmaster2012/kloosterman/MATLAB/tools/fieldtrip-20150803
dirlist  = dir('template/layout/*');
% filename = {dirlist(~[dirlist.isdir]).name}'
filename = {'CTF275.lay'}
for i=1:length(filename)
  cfg = [];
  cfg.layout = filename{i};
  layout = ft_prepare_layout(cfg);

  figure
  ft_plot_lay(layout);
  title(filename{i});
end