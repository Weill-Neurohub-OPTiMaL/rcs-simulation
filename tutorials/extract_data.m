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

% Select a single folder containing .json files and process the session
% data using Analysis-rcs-data
[unifiedDerivedTimes, timeDomainData, ~, ~, ~, ~, ~, PowerData, ~, ~, ~,...
    ~, ~, AdaptiveData, ~, ~, timeDomainSettings, powerSettings,...
    fftSettings, eventLogTable, metaData, stimSettingsOut, stimMetaData,...
    stimLogSettings, DetectorSettings, AdaptiveStimSettings, ...
    AdaptiveEmbeddedRuns_StimSettings, ~] = ProcessRCS(uigetdir, 2);
dataStreams = {timeDomainData, PowerData, AdaptiveData};
[data] = createCombinedTable(dataStreams, unifiedDerivedTimes, metaData);

% Organize and reformat the data table
data = removevars(data, {'localTime', 'TD_samplerate',...
                         'Power_ExternalValuesMask','Power_FftSize',...
                         'Power_IsPowerChannelOverrange',...
                         'Power_ValidDataMask', ...
                         'Adaptive_IsInHoldOffOnStartup',...
                         'Adaptive_PreviousAdaptiveState',...
                         'Adaptive_SensingStatus',...
                         'Adaptive_StateEntryCount',...
                         'Adaptive_StateTime','Adaptive_StimFlags',...
                         'Adaptive_StimRateInHz'});
data.Properties.VariableNames{1} = 'timestamp';
data.Properties.VariableNames{2} = 'TD0';
data.Properties.VariableNames{3} = 'TD1';
data.Properties.VariableNames{4} = 'TD2';
data.Properties.VariableNames{5} = 'TD3';
data.Properties.VariableNames{6} = 'PB1';
data.Properties.VariableNames{7} = 'PB2';
data.Properties.VariableNames{8} = 'PB3';
data.Properties.VariableNames{9} = 'PB4';
data.Properties.VariableNames{10} = 'PB5';
data.Properties.VariableNames{11} = 'PB6';
data.Properties.VariableNames{12} = 'PB7';
data.Properties.VariableNames{13} = 'PB8';
data.Properties.VariableNames{14} = 'LD_state';
data.Properties.VariableNames{15} = 'stim';


end
