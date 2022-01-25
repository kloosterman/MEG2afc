function data = interpolate_blinks(hdr, data)
% TODO add fixations

eye_events = {'EYE_BLINKS', 'EYE_SACCADES', 'EYE_FIXATIONS'};

eye_ev{1} = hdr.orig.eblink;
eye_ev{2} = hdr.orig.esacc;
eye_ev{3} = hdr.orig.efix;

% eblink = hdr.orig.eblink;
% esacc = hdr.orig.esacc;

timestamps = data.trial{1}(1,:); %ascdat.dat(1,:);

for iev = 1:length(eye_events)
  nchan = length(data.label);
  data.trial{1}(nchan+1,:) = 0; % blink binary
  data.label = [data.label eye_events{iev} ];
  
  if isempty(eye_ev{iev})
    disp(sprintf('No %s found', eye_events{iev} ))
    continue
  end
    
  % old way, not robust due to separator diffs in blink and sacc
%   arg1 = repmat({'%*s%*s%d%d'}, length(eye_ev{iev}), 1);
%   evtimes = cellfun(@sscanf, eye_ev{iev}, arg1, 'UniformOutput', false); % parse evtimes from ascdat
%   evtimes = cell2mat(cellfun(@transpose, evtimes, 'UniformOutput', false)); %transpose and turn into matrix
%   evtimes = evtimes(:,1:2);

  evtimes = cellfun(@(x) tokenize(x, ' ', true), eye_ev{iev}, 'uni', false);
  evtimes = cellfun(@(x) x{3}, evtimes, 'uni', false);
  evtimes = cellfun(@tokenize, evtimes, 'uni', false);
  evtimes = cell2mat(cellfun(@str2double, evtimes, 'uni', false));
  evtimes = evtimes(:,1:2);
  evtimes = evtimes(all(ismember(evtimes, timestamps),2),:); % check if all evtimes there
  evsmp = arrayfun(@(x) find(timestamps == x, 1,'first'), evtimes, 'uni', true ); %find sample indices of evtimes in timestamps
  if iev == 1 % some extra padding for blinks
    evsmp(:,1) = evsmp(:,1) - round(0.15*hdr.Fs); 
    evsmp(:,2) = evsmp(:,2) + round(0.15*hdr.Fs);
  end
  evsmp = evsmp(find(evsmp(:,1) > 0),:);
  for iblink=1:size(evsmp,1) % blinks but also sacc and fix
    if ~any(evsmp(iblink,:) > length(data.trial{1}))
      for ichan = 1:nchan
        data.trial{1}(ichan,evsmp(iblink,1):evsmp(iblink,2)) = ...
          linspace(data.trial{1}(ichan,evsmp(iblink,1)), ...
          data.trial{1}(ichan,evsmp(iblink,2)), ...
          length(evsmp(iblink,1):evsmp(iblink,2))); % interpolate
        
        data.trial{1}(nchan+1,evsmp(iblink,1):evsmp(iblink,2)) = 1;   % make binary chans with blinks
        
      end
    end
  end

end
  %   evsmp
  
  %                             cfg2=[]; %inspect data
  %                             cfg2.channel = 'p';
  %                             cfg2.method = 'channel';
  %                                 ft_rejectvisual(cfg2,data)
  
