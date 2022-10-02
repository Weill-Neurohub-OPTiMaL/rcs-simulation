%% Function Name: isolate_state_vector()
%
% Description: Turns the state field from openmind combinedDataTable into a
% vector of integers.
%
% Inputs:
%     table : (num_samples, 1) table
%         OpenMind output, commonly called combinedDataTable
%
% Outputs:
%     state : (1, num_samples) array
%         Vector of current states at each timepoint
%
% Author: Tanner Chas Dixon, tanner.dixon@ucsf.edu.
% Date last updated: August 24, 2022
%---------------------------------------------------------

function [state] = isolate_state_vector(table)

state = zeros(height(table), 1);
for i = 1:length(state)
    if strcmp(table.Adaptive_CurrentAdaptiveState{i}, 'No State') ...
            | isnan(table.Adaptive_CurrentAdaptiveState{i})
        state(i) = NaN;
    else
        state(i) = str2double(table.Adaptive_CurrentAdaptiveState{i}(7));
    end
end

end