function [ d_min_center, d_min_edge, r_1_min, r_2_min ] = calculate_surface_separation( r_1, r_2, L_sub_1, L_sub_2)

% This function calculates the minimum surface-to-surface separation
% between two discretized object.

% INPUTS:  r_1          (N_1 x 3) matrix containing points of all cubic lattice points of object #1
%          r_2          (N_2 x 3) matrix containing points of all cubic lattice points of object #2
%          L_sub_1      Length of a side of a subvolume in object #1 (i.e., dipole diameter) [m]
%          L_sub_2      Length of a side of a subvolume in object #1 (i.e., dipole diameter) [m]
%         
%
% OUTPUTS: d_min_center  Minimum separation from center of subvolume in object
%                        #1 to center of subvolume in object #2
%          d_min_edge    Minimum surface-to-surface separation [m]
%          r_1_min       Subvolume coordinate of object #1 for minimum separation [m]
%          r_2_min       Subvolume coordinate of object #2 for minimum separation [m]
%

% Calculate number of subvolumes in each object
[N_1, ~] = size(r_1);
[N_2, ~] = size(r_2);

% Preallocate
d_mag = zeros(N_1, N_2);

for ii = 1:N_1
    for jj = 1:N_2
        d = r_1(ii,:) - r_2(jj,:);        % Distance between subvolume i and subvolume j
        d_mag(ii,jj) = sqrt(sum(d.^2));   % Magnitude of r_ij vector
    end
end

% Calculate minimum separation distance from center-to-center of subvolumes
[d_min_center, ind] = min(d_mag);
[d_min_center, ind_2] = min(d_min_center);
ind_1 = ind(ind_2);

% Calculate minimum separation distance from edge-to-edge of subvolumes
d_min_edge = d_min_center - L_sub_1/2 - L_sub_2/2;

% Locations of minimum separation subvolumes 
r_1_min = r_1(ind_1,:);
r_2_min = r_2(ind_2,:);


        
        
        