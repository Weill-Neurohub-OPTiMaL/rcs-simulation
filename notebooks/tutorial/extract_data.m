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
[unifiedDerivedTimes, timeDomainData, ~, ~, ~, ~, ~, PowerData, ~, ~, ~,...
    ~, ~, AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
    fftSettings, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
    stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
    AdaptiveEmbeddedRuns_StimSettings, ~] = ProcessRCS(uigetdir, 2);
dataStreams = {timeDomainData, PowerData, AdaptiveData};
[combinedDataTable] = createCombinedTable(dataStreams, ...
                                          unifiedDerivedTimes, metaData);

%% Create a simplified data table

timestamp = combinedDataTable.DerivedTime/1000;

td = cell(1,4);
for i = 1:4
    td{i} = combinedDataTable(:,2+i);
end

pb = cell(1,8);
for i = 1:4
    td{i} = combinedDataTable(:,15+i);
end

ld1 = combinedDataTable.Adaptive_Ld0_output;
ld1 = correct_ld(ld1);
ld2 = combinedDataTable.Adaptive_Ld1_output;
ld2 = correct_ld(ld2);

state = isolate_state_vector(combinedDataTable);





end
