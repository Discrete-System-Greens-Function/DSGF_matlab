function [ G_0_4D, G_0_2D ] = G_0_free_space( r, omega, delta_V, epsilon_ref, N )

% This function calculates the free-space Green's function for every dipole
% interaction.  The G_0_4D output is an (N x N x 3 x 3) 4-dimensional
% matrix.  The G_0_2D output is the same output as G_0_4D, except in 2D
% matrix format.
 
% INPUTS:  r              (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%          omega          radial frequency [rad/s]
%          delta_V        (N x 1) vector of all subvolume sizes [m^3]
%          epsilon_ref    dielectric function for background reference medium (constant)
%          N              total number of subvolumes (i.e. dipoles)
%
%
% OUTPUTS: G_0_4D         free-space Green's function matrix for all subvolume interactions in 4D matrix format (N x N x 3 x 3)
%          G_0_2D         free-space Green's function matrix for all subvolume interactions in 2D matrix format (3N x 3N)


% Constants
epsilon_0 = 8.8541878e-12;    % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]


% Wave vector in background reference medium
k = omega*sqrt(epsilon_ref*epsilon_0*mu_0);

% 3-by-3 unit matrix
I = eye(3); %creating the identity matrix 

% Preallocate free-space Green's function 4D matrix
G_0_4D = zeros(N,N,3,3);

for jj = 1:N-1                         % Loop through all discretized volume locations
    for ii = jj+1:N                    % Loop through all discretized volume locations (again)
        r_ij = r(ii,:) - r(jj,:);      % Vector from location i to location j
        r_ij_mag = sqrt(sum(r_ij.^2)); % Magnitude of r_ij vector
        r_ij_hat = r_ij./r_ij_mag;     % Unit direction vector from subvolume i to j
        r_ij_outer_r_ij = (r_ij_hat')*(r_ij_hat); % Outer product of unit direction vectors
        
        % Constants
        const_1 = exp(1i*k*r_ij_mag)/(4*pi*r_ij_mag);
        const_2 = 1 - 1/((k*r_ij_mag)^2) + 1i/(k*r_ij_mag);
        const_3 = 1 - 3/((k*r_ij_mag)^2) + 3*1i/(k*r_ij_mag);
        
        % Upper triangular part of Green's function 
        G_0_4D(ii,jj,:,:) = const_1*( (const_2*I) - (const_3*r_ij_outer_r_ij) );      
    end
end




% Self-interaction terms of subvolume with itself (found from principal volume method)
% Original T-DDA principal value

a = (3*delta_V./(4*pi)).^(1/3); % Equivalent radius of subvolume
self_terms = (1./(3*delta_V*(k^2))).*(2*(exp(1i*k*a).*(1 - 1i*a*k) - 1) - 1);
self_terms = diag(repelem(self_terms, 1, 3));
%self_terms = diag(repelem(self_terms, 3, 1));

% Convert from 4D matrix to a 2D matrix
G_0_2D = reshape(permute(G_0_4D,[3,1,4,2]),[3*N,3*N]);
G_0_2D = G_0_2D + G_0_2D.' + self_terms;

