function [sens] = MEG2afc_sensorselection()

if ismac
    load('/Users/kloosterman/gridmaster2012/kloosterman/MATLAB/MEG_HH_analysis/plotting/chlabel.mat')
else
    load('/home/mpib/kloosterman/MATLAB/MEG_HH_analysis/plotting/chlabel.mat')
end

% % motor based on sensor naming (MLC etc)
% leftmotor =  { 'MLC11';'MLC12'  ;'MLC13' ;'MLC14';'MLC15' ;'MLC16';'MLC22';'MLC23';'MLC24';'MLC25';'MLC31';'MLC32';'MLC41';'MLC42';'MLC51';'MLC52'; ...
%     'MLC53';'MLC54';'MLC55';'MLC61';'MLC62';'MLC63'};
% 
% rightmotor = {   'MRC11';'MRC12';'MRC13';'MRC14';'MRC15';'MRC16';'MRC17';'MRC21';'MRC22';'MRC23';'MRC24';'MRC25';'MRC31';'MRC32';'MRC41';'MRC42'; ...
% 'MRC51';'MRC52';'MRC53';'MRC54';'MRC55';'MRC61';'MRC62';'MRC63'};
% 
% motor = [leftmotor; rightmotor; 'MZC01'; 'MZC02'; 'MZC03';  ];

% motor based on resp1-resp2: left-right button press
leftmotor = {    'MLC24' ;'MLC23'; 'MLC16'; 'MLC22'; 'MLC15'; 'MLC31'; 'MLC25'; 'MLC14'; 'MLP45'; 'MLC32'; 'MLP57'; 'MLF65'; 'MLF66'; 'MLC41'; 'MLF64'; 'MLC42'; 'MLP35'; 'MLF63' ;'MLP56'};
rightmotor = {  'MRC24' ;'MRC23'; 'MRC31'; 'MRC22'; 'MRC41'; 'MRC25'; 'MRC42'; 'MRC32'; 'MRC16'; 'MRC15'; 'MRC14'; 'MRP45'; 'MRC17'; 'MRP35'; 'MRC53'; 'MRP23'; 'MRC21'; 'MRC13'; 'MRC54'; 'MRF66'};
motor = [leftmotor; rightmotor];

leftmotorind = match_str(chlabel,ft_channelselection(leftmotor,chlabel));
rightmotorind = match_str(chlabel,ft_channelselection(rightmotor,chlabel));
motorind =  match_str(chlabel,ft_channelselection(motor,chlabel));

% %occipital left right based on sensor naming
% leftocc = { 'MLO11'; 'MLO12'; 'MLO14';  'MLO21';  'MLO22'; 'MLO23'; 'MLO24'; 'MLO31'; 'MLO32'; 'MLO33'; 'MLO34'; 'MLO41'; 'MLO42'; 'MLO43'; 'MLO44'}  ; %'MLO51'; 'MLO52'; 'MLO53'
% rightocc = { 'MRO12'; 'MRO13'; 'MRO14'; 'MRO21'; 'MRO22'; 'MRO23'; 'MRO24'; 'MRO31'; 'MRO32'; 'MRO33'; 'MRO34'; 'MRO41'; 'MRO43'; 'MRO44'  }; % 'MRO51'; 'MRO52'; 'MRO53'
% occ = [leftocc; rightocc; 'MZO01'; 'MZO02'; 'MZO03' ];

%occipital left right based on 40-80 Hz 0-1 s power across all conditions
% 30 sensors, split based on n sensors
% % occ = { 'MLO11', 'MLP52', 'MLO12', 'MLO21', 'MZO01', 'MLP51', 'MLO22', 'MLO24', 'MRP51', 'MLP41', 'MLP53', 'MLO23', 'MZP01', 'MZO02', 'MZC04', 'MRP52', ...
% %     'MRO12', 'MLP31', 'MRP33', 'MLO31', 'MRP31', 'MRO13', 'MLO32', 'MRP11', 'MLO14', 'MRO23', 'MRO21', 'MLT27', 'MLP54', 'MRO22'};
% leftocc = { 'MLT27';  'MLO14'; 'MLO24'; 'MLP54'; 'MLO23'; 'MLP53'; 'MLP41'; 'MLO11'; 'MLP52'; 'MLO12'; 'MLO22'; 'MLP31'; 'MLO31';'MLO32';};
% rightocc ={ 'MZO01';  'MRP51'; 'MZP01'; 'MZO02'; 'MRP52'; 'MRO12'; 'MRP33'; 'MRP31'; 'MRO13'; 'MRP11'; 'MRO23'; 'MRO21'; 'MRO22'};
% occ = [leftocc; rightocc; 'MLO21'; 'MLP51';  'MZC04'];

% occipital based on 40 max gamma sensors 40-80 Hz 0.25-0.5 s wrt stim  across all conditions
% split based on Z, L and R names. most objective
leftocc = {    'MLP52'; 'MLO11'; 'MLO12'; 'MLP51'; 'MLO21'; 'MLO22'; 'MLO24'; 'MLP41'; 'MLP53'; 'MLP31'; 'MLO23'; 'MLO14'; 'MLP54'; 'MLT27'; 'MLP21'; 'MLO31'; 'MLO34'; 'MLO32'; 'MLP32'; 'MLP42'; 'MLO33'; 'MLT16'};
rightocc = {   'MRP51'; 'MRP52'; 'MRO12'; 'MRP31'; 'MRO13'; 'MRO23'; 'MRO24'; 'MRP41'; 'MRO21'; 'MRP53'; 'MRP21'; 'MRO33'; 'MRO22'; 'MRO32'; 'MRO14'; };
occ = [leftocc; rightocc; 'MZO01'; 'MZP01'; 'MZO02']; 
% 
% occipital based on leftstim-rightstim
% leftocc = {    'MLO12'; 'MLO11'; 'MLO23'; 'MLP52'; 'MLO24';'MLO21'; 'MLP54'; 'MLP42'; 'MLO42'; 'MLO32'; 'MLO22';  'MLP53'; 'MLO31' };      
% rightocc = {     'MRO23'; 'MRO13'; 'MRP55'; 'MRO33'; 'MRO24'; 'MRO32'; 'MRP52'; 'MRO44'; 'MRO14'; 'MRO12'; 'MRP54'; 'MRO34'; 'MRT16'; 'MRO22'; 'MRP53' };
% occ = [leftocc; rightocc];

% occipital based on leftchoice-rightchoice
% TODO make poolings
%     'MLT52'    'MRO31'
%     'MLP54'    'MRO41'
%     'MLP32'    'MRO22'
%     'MLO14'    'MLT53'
%     'MLT51'    'MRO21'
%     'MLP53'    'MRO23'
%     'MLF56'    'MRO32'
%     'MLP34'    'MZO02'
%     'MLP55'    'MRO24'
%     'MLP43'    'MRO12'

leftoccind = match_str(chlabel,ft_channelselection(leftocc,chlabel));
rightoccind = match_str(chlabel,ft_channelselection(rightocc,chlabel));
occind =  match_str(chlabel,ft_channelselection(occ,chlabel));

%frontal lateralized theta cluster
frontal = {'MLC13', 'MLC14', 'MLC15', 'MLC22', 'MLC23', 'MLF12', 'MLF13', 'MLF14', 'MLF22', 'MLF23', 'MLF24', 'MLF25', 'MLF32', 'MLF33', 'MLF34', 'MLF35', 'MLF43', 'MLF44', 'MLF45', 'MLF46', 'MLF53', 'MLF54', 'MLF55', 'MLF56', 'MLF62', 'MLF63', 'MLF64', 'MLF65', 'MLT11', 'MLT12', 'MLT21', 'MLT22', 'MLT31', 'MLT32', 'MLT41', 'MRF13', 'MRF14', 'MRF24', 'MRF25', 'MRF34', 'MRF35', 'MRF45', 'MRF46', 'MRF55', 'MRF56', 'MRF65', 'MRF66', 'MRT11', 'MRT12', 'MRT13', 'MRT21', 'MRT22', 'MRT23', 'MRT31', 'MRT32', 'MRT33', 'MRT41', 'MRT42', 'MRT51'};

% frontal gamma cluster that seems to differ between easy and hard
% frontal = {'MLC14', 'MLC22', 'MLC23', 'MLC31', 'MLC32', 'MLC41', 'MLC42', 'MLC51', 'MLC52', 'MLC53', 'MLC54', 'MLC55', 'MLC61', 'MLC62', 'MLC63', 'MLP12', 'MLP23', 'MRC21', 'MRC41', 'MRC51', 'MRC52', 'MRC53', 'MRC54', 'MRC55', 'MRC61', 'MRC62', 'MRC63', 'MZC02', 'MZC03', 'BG1', 'BG2', 'BG3', 'BP1', 'BP2', 'BP3', 'BR1', 'BR2', 'BR3', 'G11', 'G12', 'G13', 'G22', 'G23', 'P11', 'P12', 'P13', 'P22', 'P23', 'Q11', 'Q12', 'Q13', 'Q22', 'Q23', 'R11', 'R12', 'R13', 'R22', 'R23'};
frontalind =  match_str(chlabel,ft_channelselection(frontal,chlabel));

% large frontoparietal cluster with strong baseline drug effect, see MEG2afc_prestim_spectrum
frontal_drug = {'MLC11';'MLC12';'MLC13';'MLC14';'MLC15';'MLC16';'MLC22';'MLC23';'MLC51';'MLC52';'MLF11';'MLF12';'MLF13';'MLF14';'MLF22';'MLF23';'MLF24';'MLF25';'MLF31';'MLF32';'MLF33';'MLF34';'MLF35';'MLF41';'MLF42';'MLF43';'MLF44';'MLF45';'MLF46';'MLF51';'MLF52';'MLF53';'MLF54';'MLF55';'MLF56';'MLF61';'MLF62';'MLF63';'MLF64';'MLF65';'MLF66';'MLT11';'MLT12';'MLT13';'MLT21';'MLT22';'MLT23';'MLT24';'MLT25';'MLT31';'MLT32';'MLT33';'MLT34';'MLT35';'MLT36';'MLT41';'MLT42';'MLT43';'MLT51';'MLT52';'MLT53';'MRC11';'MRC12';'MRC13';'MRC14';'MRC15';'MRC16';'MRC21';'MRC22';'MRC23';'MRC24';'MRC31';'MRC41';'MRC42';'MRC51';'MRC52';'MRF11';'MRF12';'MRF13';'MRF14';'MRF21';'MRF22';'MRF23';'MRF24';'MRF25';'MRF31';'MRF32';'MRF33';'MRF34';'MRF35';'MRF41';'MRF42';'MRF43';'MRF44';'MRF45';'MRF46';'MRF51';'MRF52';'MRF53';'MRF54';'MRF55';'MRF56';'MRF61';'MRF63';'MRF64';'MRF65';'MRF66';'MRT11';'MRT12';'MRT13';'MRT14';'MRT21';'MRT22';'MRT23';'MRT24';'MRT25';'MRT26';'MRT31';'MRT32';'MRT33';'MRT34';'MRT35';'MRT36';'MRT37';'MRT41';'MRT42';'MRT43';'MRT44';'MRT45';'MRT46';'MRT51';'MRT52';'MRT53';'MRT55';'MRT56';'MZC01';'MZC02';'MZF01';'MZF02';'MZF03'};
frontal_drug_ind =  match_str(chlabel,ft_channelselection(frontal_drug,chlabel));

% occipital region that shows power decrease
occipital_drug = {'MLC25';'MLC32';'MLC42';'MLC54';'MLC55';'MLC63';'MLO11';'MLO12';'MLP11';'MLP12';'MLP21';'MLP22';'MLP23';'MLP31';'MLP32';'MLP33';'MLP34';'MLP35';'MLP41';'MLP42';'MLP43';'MLP44';'MLP51';'MLP52';'MLP53';'MLP54';'MRC54';'MRC55';'MRC62';'MRC63';'MRO12';'MRO13';'MRP11';'MRP12';'MRP21';'MRP22';'MRP23';'MRP31';'MRP32';'MRP33';'MRP34';'MRP41';'MRP42';'MRP43';'MRP44';'MRP51';'MRP52';'MRP53';'MRP54';'MRP55';'MZC04';'MZO01';'MZP01'};
occipital_drug_ind =  match_str(chlabel,ft_channelselection(occipital_drug,chlabel));

sens=[];
sens.ind = {occind; motorind; (1:length(chlabel))'; frontalind; frontal_drug_ind; occipital_drug_ind};
sens.leg = {'occipital'; 'motor'; 'allsensors'; 'frontal'; 'frontal_drug'; 'occipital_drug'};


% make lateralization per channel vectors
left_chan = ft_channelselection('ML*', chlabel);

% go through left chan, find right counterpart
LR_subtract_mat = [];
for ichan = 1:length(left_chan)
    left_chan_name = left_chan{ichan};
    right_chan_name = left_chan_name;
    right_chan_name(2) = 'R';

    left_chan_ind = find(strcmp(left_chan_name, chlabel)); % get ind in chlabel
    right_chan_ind = find(strcmp(right_chan_name, chlabel)); % get ind in chlabel
    if isempty(right_chan_ind)
        fprintf('Matching counterpart of %s not found in %s\n', left_chan_name, right_chan_name)
        LR_subtract_mat(ichan,1) = NaN;
        LR_subtract_mat(ichan,2) = NaN;
        continue
    end
    % make leftrightmat
    LR_subtract_mat(ichan,1) = left_chan_ind;
    LR_subtract_mat(ichan,2) = right_chan_ind;
end
sens.LR_subtract_mat = LR_subtract_mat(~isnan(LR_subtract_mat(:,1)),:); % remove nan channels
% inspect matching:
% [chlabel(LR_subtract_mat(:,1))'; chlabel(LR_subtract_mat(:,2))' ]

chancmb = {'MLC11','MRC11';'MLC12','MRC12';'MLC13','MRC13';'MLC14','MRC14';'MLC15','MRC15';'MLC16','MRC16';'MLC17','MRC17';'MLC21','MRC21';'MLC22','MRC22';'MLC23','MRC23';'MLC24','MRC24';'MLC25','MRC25';'MLC31','MRC31';'MLC32','MRC32';'MLC41','MRC41';'MLC42','MRC42';'MLC51','MRC51';'MLC52','MRC52';'MLC53','MRC53';'MLC54','MRC54';'MLC55','MRC55';'MLC61','MRC61';'MLC62','MRC62';'MLC63','MRC63';'MLF11','MRF11';'MLF12','MRF12';'MLF13','MRF13';'MLF14','MRF14';'MLF21','MRF21';'MLF22','MRF22';'MLF23','MRF23';'MLF24','MRF24';'MLF25','MRF25';'MLF31','MRF31';'MLF32','MRF32';'MLF33','MRF33';'MLF34','MRF34';'MLF35','MRF35';'MLF41','MRF41';'MLF42','MRF42';'MLF43','MRF43';'MLF44','MRF44';'MLF45','MRF45';'MLF46','MRF46';'MLF51','MRF51';'MLF52','MRF52';'MLF53','MRF53';'MLF54','MRF54';'MLF55','MRF55';'MLF56','MRF56';'MLF61','MRF61';'MLF62','MRF62';'MLF63','MRF63';'MLF64','MRF64';'MLF65','MRF65';'MLF66','MRF66';'MLF67','MRF67';'MLO11','MRO11';'MLO12','MRO12';'MLO13','MRO13';'MLO14','MRO14';'MLO21','MRO21';'MLO22','MRO22';'MLO23','MRO23';'MLO24','MRO24';'MLO31','MRO31';'MLO32','MRO32';'MLO33','MRO33';'MLO34','MRO34';'MLO41','MRO41';'MLO42','MRO42';'MLO43','MRO43';'MLO44','MRO44';'MLO51','MRO51';'MLO52','MRO52';'MLO53','MRO53';'MLP11','MRP11';'MLP12','MRP12';'MLP21','MRP21';'MLP22','MRP22';'MLP23','MRP23';'MLP31','MRP31';'MLP32','MRP32';'MLP33','MRP33';'MLP34','MRP34';'MLP35','MRP35';'MLP41','MRP41';'MLP42','MRP42';'MLP43','MRP43';'MLP44','MRP44';'MLP45','MRP45';'MLP51','MRP51';'MLP52','MRP52';'MLP53','MRP53';'MLP54','MRP54';'MLP55','MRP55';'MLP56','MRP56';'MLP57','MRP57';'MLT11','MRT11';'MLT12','MRT12';'MLT13','MRT13';'MLT14','MRT14';'MLT15','MRT15';'MLT16','MRT16';'MLT21','MRT21';'MLT22','MRT22';'MLT23','MRT23';'MLT24','MRT24';'MLT25','MRT25';'MLT26','MRT26';'MLT27','MRT27';'MLT31','MRT31';'MLT32','MRT32';'MLT33','MRT33';'MLT34','MRT34';'MLT35','MRT35';'MLT36','MRT36';'MLT37','MRT37';'MLT41','MRT41';'MLT42','MRT42';'MLT43','MRT43';'MLT44','MRT44';'MLT45','MRT45';'MLT46','MRT46';'MLT47','MRT47';'MLT51','MRT51';'MLT52','MRT52';'MLT53','MRT53';'MLT54','MRT54';'MLT55','MRT55';'MLT56','MRT56';'MLT57','MRT57'}

% %% Freq bands of interest
% % faxis = 5:2:149;
% beta_band_ind = faxis >= 12 & faxis <= 35;
% faxis(beta_band_ind)
% % gamma_band_ind = faxis >= 40 & faxis <= 80;
% % gamma_band_ind = faxis >= 56 & faxis <= 90;
% gamma_band_ind = faxis >= 40 & faxis <= 65;
% faxis(gamma_band_ind)
% bands = {'beta' 'gamma'};


