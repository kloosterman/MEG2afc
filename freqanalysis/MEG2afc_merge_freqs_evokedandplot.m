
if ismac
    basepath = '/Users/kloosterman/gridmaster2012/kloosterman'; % on the cluster
else
    basepath = '/home/mpib/kloosterman'; % on the cluster
end
%         PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq_robdetr', analysistype{iband}, trigger_leg{itrig});
% PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq_robdetr/lowfine/baseline')
PREIN = fullfile(basepath, 'projectdata/MEG2afc/freq_robdetr/low/stim')
cd(PREIN)

SUBJ  = {
    'NK1' ... has funny driftrate
     'NK2' ...
    'NK3' ...
    'NK4'   'NK5'   'NK7'   'NK8'  'NK9'  ...
    'NK11' ...
    'NK12'   'NK13'  ...
     'NK14'  ...
    'NK15'...
    'NK16'     'NK17'      'NK18'      'NK19'     'NK20'       'NK21'
    }; %  'NK6' out,  'NK10' don't exist,   %%   'NK12' 'NK15' had bad ddm fits in the past

drugs = {'drug' 'plac'};
motors = {'contra' 'ipsi'};

allfreq=[];
for isub=1:length(SUBJ)
    cd(SUBJ{isub})
    disp(SUBJ{isub})
    
    for idrug=1:2
        for imo = 1:2
            folder = [drugs{idrug}  '_' motors{imo} ];
            disp(folder)
            if ~exist(folder)
                warning('not found')
                continue
            end
            cd(folder)
            list = dir('*evoked_freq.mat');
            if isempty(list)
                warning('not found')
                cd ..
                continue
            end
                
            load(list.name)
%             warning off
            freq = ft_freqdescriptives([], freq);
%             warning on
            allfreq{isub, idrug, imo} = freq;
            cd ..
            
        end
    end
    cd ..
end

% Grand average
cfg=[];
cfg.keepindividual = 'yes';
cfg.foilim = [1 33];

respavg = [];
for idrug = 1:2
    temp = {allfreq{:,idrug,:}};
    temp = temp(~cellfun(@isempty,temp));
    respavg{idrug} = ft_freqgrandaverage(cfg, temp{:} );
end
respavg{3} = respavg{1}; % put in constrast
respavg{3}.powspctrm = (respavg{1}.powspctrm - respavg{2}.powspctrm) ./ respavg{2}.powspctrm * 100;

%% plot spectra for lowfine
sens = MEG2afc_sensorselection();
cols = {'r', 'b', 'g', 'k'};
drugs = {'drug' 'plac' 'drug wrt plac'};

% close all
figure; hold on
clear h
for idrug = 3 % 1:2
    dat = mean(respavg{idrug}.powspctrm(:,sens.ind{1},:),2)    
    h(idrug) = shadedErrorBar(respavg{idrug}.freq, squeeze(mean(dat)), squeeze(std(dat))/sqrt(32), cols{idrug},1 )
end

legend([h.mainLine], drugs(idrug))
ssvep = 60/7;
plot([ssvep ssvep], get(gca, 'YLim'), 'k')
plot([ssvep*2 ssvep*2], get(gca, 'YLim'), 'k')
title('Evoked power spectra')


%% plot TFR for stimlocked
load colormap_jetlightgray.mat
cfg = [];
cfg.layout = 'CTF275.lay';       %CTF275_RH.lay neuromag306cmb neuromag306all neuromag306mag neuromag306planar CTF275_helmet.mat
cfg.hotkeys = 'yes';
cfg.fontsize = 18;
cfg.colormap = cmap;
cfg.zlim = 'maxabs';
cfg.xlim = [-0.5 1.5]
cfg.colorbar = 'yes';
f = figure;
f.Position = [ 680   678   1200   1200];
ft_multiplotTFR(cfg, respavg{2});

