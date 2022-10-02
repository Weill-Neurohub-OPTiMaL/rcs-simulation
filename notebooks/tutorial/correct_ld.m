%% Function Name: correct_ld()
%
% Description: Corrects LD outputs by accounting for (1) incorrectly  
% assigned negative bit during integer interpretation and (2) fractional   
% fixed pointvalue.
%
% Inputs:
%     ld : (num_samples, 1) array
%         LD outputs
%
% Outputs:
%     ld_corrected
%
% Author: Tanner Chas Dixon, tanner.dixon@ucsf.edu.
% Date last updated: August 24, 2022
%---------------------------------------------------------

function [ld] = correct_ld(ld)

ld(ld>2^31) = ld(ld>2^31) - 2^32;
ld = ld/(2^10);

end
