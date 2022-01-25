function respavg = MEG2afc_loadtimelocks()

% ALL Subjects:
SUBJ  = {
  'NK1' ... has funny driftrate
  'NK2' ...
  'NK3'   'NK4'   'NK5'   'NK7'   'NK8'  'NK9'   'NK11'  'NK12'   'NK13'    'NK15'...
  'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
  }; %   'NK14' nomri 'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

%   SUBJ  = {
%    'NK4'   'NK7'
%     }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past
nsub = length(SUBJ);
PREIN = '/Users/kloosterman/beegfs/projectdata/MEG2afc/timelock'; % on the cluster
% PREIN = '/Volumes/FB-LIP/user/Niels/MEG2afc/timelock'; % on the cluster

respavg = cell(nsub,3,2,4,4);
for isub = 1:length(SUBJ)
  curpath = fullfile(PREIN, [SUBJ{isub} '_timelock.mat' ]);
  disp(curpath)
  load(curpath)
    
  respavg(isub,:,:,:,:) = timelockout;  %  itrig, idrug, idiff
end

cfg=[];
cfg.keepindividual = 'yes';

timelockkeep = {};
for iltr = 1:3
  for itrig = 1:2
    for idrug = 1:4
      
      for idiff = 1:3
        dat = {respavg{:, iltr, itrig, idrug, idiff }};
        timelockkeep{iltr, itrig, idrug, idiff} =  ft_timelockgrandaverage(cfg, dat{:} );
      end
      
    end
  end
end
respavg = timelockkeep;

