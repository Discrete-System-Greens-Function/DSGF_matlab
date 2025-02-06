function [ r_centered ] = center_of_mass( r )

% This function translates a discretized object so that its center of mass
% is located at the origin [0,0,0].

%
% INPUTS:  r              (N x 3) cubic lattice discretization
%
% OUTPUTS: r_centered     Discretization with center-of-mass at the origin
%            

[N,~] = size(r);                                % Number of subvolumes in discretization
center_of_mass = mean(r);                       % Calculate center of mass coordinates
r_centered = r - repmat(center_of_mass, N, 1);  % Put center of mass at origin

end