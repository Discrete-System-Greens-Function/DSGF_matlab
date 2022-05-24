function [A_2D, alpha_0 ] = A_matrix_function( G_0_2D, omega, epsilon, epsilon_ref, delta_V, N )

% This function populates the deterministic interaction matrix for the
% system Green's function approach.

% INPUTS:  G_0_2D           free-space Green's function matrix for all subvolume interactions in 2D matrix format (3N x 3N)
%          omega            radial frequency [rad/s]
%          epsilon      	(N x 1) vector of dielectric functions for each subvolume
%          epsilon_ref      dielectric function for background reference medium (constant)
%          delta_V          (N x 1) vector of all subvolume sizes [m^3]
%          N                total number of subvolumes (i.e. dipoles)
%
%
% OUTPUTS: A_2D             (3N x 3N) interactive matrix
%          alpha_0          bare radiative polarizability of each subvolume (i.e. dipole) [m^3]


% Constants
epsilon_0 = 8.8542e-12;  % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);   % Permeability of free space [H/m]

% Wave vector in free space
k_0 = omega*sqrt(epsilon_0*mu_0);

% Bare polarizability [m^3]
alpha_0 = delta_V.*(epsilon - epsilon_ref);

% Calculate full interaction matrix
A_2D = eye(3*N) - (k_0^2)*G_0_2D*diag((repelem(alpha_0,3)));








