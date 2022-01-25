function [blinkdat] = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [BLINKS_NR,PUPIL_RAW_B,PUPIL_RAW_B_PSC,RT,RUN_NR,TRIAL_NR] =
%   IMPORTFILE(FILENAME) Reads data from text file FILENAME for the default
%   selection.
%
%   [BLINKS_NR,PUPIL_RAW_B,PUPIL_RAW_B_PSC,RT,RUN_NR,TRIAL_NR] =
%   IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [blinks_nr,pupil_raw_b,pupil_raw_b_psc,rt,run_nr,trial_nr] = importfile('nk1_B_pupil_data.csv',2, 961);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/09/20 15:15:01

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column2: double (%f)
%	column3: double (%f)
%   column4: double (%f)
%	column5: double (%f)
%   column6: double (%f)
%	column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
blinks_nr = dataArray{:, 1};
pupil_raw_b = dataArray{:, 2};
pupil_raw_b_psc = dataArray{:, 3};
rt = dataArray{:, 4};
run_nr = dataArray{:, 5};
trial_nr = dataArray{:, 6};





%% Allocate imported array to column variable names
nruns = run_nr(end);
blinkdat = cell(1,nruns);

for irun = 1:nruns
    
    run_ind = run_nr == irun;

    blinkdat{irun}.trial_nr = trial_nr(run_ind);
    blinkdat{irun}.blinks_nr = blinks_nr(run_ind);
    blinkdat{irun}.rt = rt(run_ind);
    blinkdat{irun}.baseline_pupil = pupil_raw_b(run_ind);
    blinkdat{irun}.baseline_pupil_psc = pupil_raw_b_psc(run_ind);

end
