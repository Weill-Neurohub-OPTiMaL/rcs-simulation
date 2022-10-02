%% Function Name: isolate_stim_vector()
%
% Description: Turns the stim field from openmind combinedDataTable into a
% vector
%
% Inputs:
%     table : (num_samples, 1) table
%         OpenMind output, commonly called combinedDataTable
%
% Outputs:
%     stim : (1, num_samples) array
%         Vector of current stim at each timepoint
%
% Author: Tanner Chas Dixon, tanner.dixon@ucsf.edu.
% Date last updated: September 1, 2022
%---------------------------------------------------------

function [stim] = isolate_stim_vector(table)

stim = zeros(height(table), 1);
for i = 1:length(stim)
    if isnan(table.Adaptive_CurrentProgramAmplitudesInMilliamps{i})
        stim(i) = NaN;
    else
        stim(i) = table.Adaptive_CurrentProgramAmplitudesInMilliamps{i}(1);
    end
end

end