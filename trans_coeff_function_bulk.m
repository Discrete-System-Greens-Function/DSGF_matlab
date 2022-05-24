function [ Trans_spectral_bulk ] = trans_coeff_function_bulk( Trans, ind_bulk_vector, ind_1, ind_2 )

% This function calculates the transmission coefficient between bulk
% objects given the transmission coefficient between every pair of dipoles
% for a given frequency.

% INPUTS:  Trans            N x N matrix of transmission coefficient between
%                           every pair of dipoles for a given frequency
%                           [unitless]
%          ind_bulk_vector  Vector containing the index of the start of every bulk object in the r array of subvolume coordinates
%          ind_1            Number of bulk object that's emitting
%          ind_2            Number of bulk object that's absorbing
%
%
%
% OUTPUTS: Trans_spectral_bulk    (N x N) transmission coefficient vector
%                                 Note: Trans_bulk(ind_1 --> ind_2)

if ind_1 == length(ind_bulk_vector) && ind_2 == length(ind_bulk_vector)
    Trans_spectral_bulk = sum(sum(Trans(ind_bulk_vector(ind_1):end, ind_bulk_vector(ind_2):end)));
elseif ind_1 == length(ind_bulk_vector) && ind_2 ~= length(ind_bulk_vector)
    Trans_spectral_bulk = sum(sum(Trans(ind_bulk_vector(ind_1):end, ind_bulk_vector(ind_2):ind_bulk_vector(ind_2 + 1)-1)));
elseif ind_1 ~= length(ind_bulk_vector) && ind_2 == length(ind_bulk_vector)
    Trans_spectral_bulk = sum(sum(Trans(ind_bulk_vector(ind_1):ind_bulk_vector(ind_1 + 1)-1, ind_bulk_vector(ind_2):end)));
else
    Trans_spectral_bulk = sum(sum(Trans(ind_bulk_vector(ind_1):ind_bulk_vector(ind_1 + 1)-1, ind_bulk_vector(ind_2):ind_bulk_vector(ind_2 + 1)-1)));
end