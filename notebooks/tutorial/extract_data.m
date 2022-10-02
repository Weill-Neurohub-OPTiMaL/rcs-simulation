%% Function Name: extract_data()
%
% Description: Converts .json files from RC+S recordings into .csv files
% and MATLAB tables that can be used with the rcssim module. 
%
% The Analysis-rcs-data repo must be included in the path to run 
% (https://github.com/openmind-consortium/Analysis-rcs-data/tree/master)
%
% Inputs:
%     user input : select a single folder containing the .json files for
%         your recording
%
% Outputs:
%     data : (X, n) table
%         FFT amplitude data given in internal RC+S units, or converted to 
%         match the scaled mV units that the device outputs in data logs if
%         specified by the `output_in_mv` parameter.
%     settings : (a, b) table
%         Unix timestamps for the corresponding FFT data samples.
%
% Author: Tanner Chas Dixon, tcdixon44@gmail.com
% Date last updated: May 27, 2022
%---------------------------------------------------------
 
function [data, settings] = extract_data()

%% Select a single folder containing .json files and process the session
% data using Analysis-rcs-data
data_dir = uigetdir;
[unifiedDerivedTimes, timeDomainData, ~, ~, ~, ~, ~, PowerData, ~, ~, ~,...
    ~, ~, AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
    fftSettings, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
    stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
    AdaptiveEmbeddedRuns_StimSettings, ~] = ProcessRCS(data_dir, 2);
dataStreams = {timeDomainData, PowerData, AdaptiveData};
[combinedDataTable] = createCombinedTable(dataStreams, ...
                                          unifiedDerivedTimes, metaData);

%% Create a simplified data table

% timestamp
timestamp = combinedDataTable.DerivedTime/1000;

% Time-Domain data
td1 = combinedDataTable.TD_key0;
td2 = combinedDataTable.TD_key1;
td3 = combinedDataTable.TD_key2;
td4 = combinedDataTable.TD_key3;

% Power-Band data
pb1 = combinedDataTable.Power_Band1;
pb2 = combinedDataTable.Power_Band2;
pb3 = combinedDataTable.Power_Band3;
pb4 = combinedDataTable.Power_Band4;
pb5 = combinedDataTable.Power_Band5;
pb6 = combinedDataTable.Power_Band6;
pb7 = combinedDataTable.Power_Band7;
pb8 = combinedDataTable.Power_Band8;

% LD outputs
ld1 = combinedDataTable.Adaptive_Ld0_output;
ld1 = correct_ld(ld1);
ld2 = combinedDataTable.Adaptive_Ld1_output;
ld2 = correct_ld(ld2);

% state
state = isolate_state_vector(combinedDataTable);

% stim
stim = isolate_stim_vector(combinedDataTable);

% combine into a single data table and write to file
data = table(timestamp, ...
             td1, td2, td3, td4, ...
             pb1, pb2, pb3, pb4, pb5, pb6, pb7, pb8, ...
             ld1, ld2, state, stim);
writetable(data, fullfile(data_dir, 'data.csv'))


%% Create a settings file
% If settings were changed in the middle of the selected session, this may
% not work appropriately.
fs_td = timeDomainSettings.samplingRate(1);
fft_size = fftSettings.fftConfig.size;
interval = fftSettings.fftConfig.interval;
bit_shift = str2num(fftSettings.fftConfig.bandFormationConfig(6));

band_edges_hz = '[';
for i = 1:8
    td_ch = num2str(ceil(i/2));
    band = strcat('[', td_ch, ',[', ...
                  num2str(powerSettings.powerBands(1).lowerBound(i)), ...
                  ',', ...
                  num2str(powerSettings.powerBands(1).upperBound(i)), ...
                  ']],');
    band_edges_hz = strcat(band_edges_hz, band);
end
band_edges_hz(end) = ']';

% update_rate = 


end
