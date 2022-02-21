function ddmpars = MEG2afc_load_ddm_paras(csvfile)
% return DDM para's from csv files:

fid = fopen(csvfile); % get header
a = textscan(fid,'%s',1);
fclose(fid);
parnames = string(strsplit(a{1}{1}, ','));

temp = csvread(csvfile, 1,0); % get data
% parnames = textread(csvfile, ',', 0,0, [0 0 0 size(temp,2)]); 
parnames = strrep(parnames, '(', '_');
parnames = strrep(parnames, ')', '');
parnames = strrep(parnames, '.', '_');

% parnames = string({'a', 'v', 't', 'z', 'dc', 'bic', 'likelihood', 'penalty'});
for ipar = 1:length(parnames)
  if isempty(parnames{ipar})
    continue
  end
%     colinds = startsWith(hdr, parnames{ipar});
%     siz = sum(colinds);
%     if siz == 20
%         ddmpars.(parnames{ipar}) = reshape( temp(:, colinds), [], 10, 2); 
%         ddmpars.dimord = 'subj_bin_cond';
%     else
        ddmpars.(parnames{ipar}) = temp(:, ipar);
%     end
end

