function [ Trans ] = trans_coeff_function( G_sys_2D, omega, epsilon, delta_V, N )

% This function calculates the transmission coefficient for every dipole
% for the system Green's function approach.

% INPUTS:  G_sys_2D      system Green's function matrix for all subvolume interactions (3N x 3N)
%          omega         radial frequency [rad/s]
%          epsilon       (N x 1) vector of dielectric functions for each subvolume
%          delta_V       (N x 1) vector of all subvolume sizes [m^3]
%          N             total number of subvolumes (i.e. dipoles)
%
%
% OUTPUTS: Trans         (N x N) transmission coefficient matrix


% Constants
epsilon_0 = 8.8542e-12;  % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);   % Permeability of free space [H/m]

% Wave vector in free space
k_0 = omega*sqrt(epsilon_0*mu_0);

% Preallocate matrices
Trans = zeros(N, N);     % Preallocate transmission coefficient matrix


for j = 1:N              % Loop through all discretized volume locations
    for i = j:N          % Loop through all discretized volume locations (again)
        
        % Extract one Green's function component
        G_element = G_sys_2D((3*i)-2:3*i, (3*j)-2:3*j);
        
        % Transmission coefficient for dipole i with dipole j (3x3 matrix for every i-j dipole combination)
        Trans(i,j) = 4*(k_0^4)*delta_V(i)*delta_V(j)*imag(epsilon(i))*imag(epsilon(j))*sum(diag((G_element*(G_element'))));
        Trans(j,i) = Trans(i,j);
        
    end
end





