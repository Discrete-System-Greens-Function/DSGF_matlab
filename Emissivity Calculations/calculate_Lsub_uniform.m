function [ L_sub, delta_V ] = calculate_Lsub_uniform( r )

% This function calculates the subvolume side length (i.e., cubic lattice
% discretization length scale) for a uniform cubic lattice discretization
% r.

%
% INPUTS:  r         (N x 3) uniform cubic lattice discretization
%
% OUTPUTS: L_sub     Subvolume side length (i.e., cubic lattice discretization length scale)
%          delta_V   Subvolume volume (= L_sub^3)
%

for ii = 1:3 % Loop through all columns

    % All differences
    D = abs(bsxfun(@minus,r(:,ii), r(:,ii).')); % square form;

    % Find nonzero minimum difference
    L_sub = min(D(D>0));

    if isempty(L_sub)
        continue
    else
        break
    end

end % End loop through columns

% Subvolume volume
delta_V = L_sub^3;

