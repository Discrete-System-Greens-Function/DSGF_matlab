function [ G_sys_2D, Trans, Q_omega_bulk, Q_omega_subvol ] = direct_function( omega, r, epsilon, epsilon_ref, delta_V, T_vector, ind_bulk )

% This function calculates the system Green's function, the monochromatic 
% thermal power dissipated in each subvolume, and the monochromatic thermal
% power dissipated in each bulk object.

% INPUTS:  omega            Radial frequency [rad/s]
%          r                (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%          epsilon          (N x 1) vector of all subvolume dielectric functions
%          epsilon_ref      Dielectric function for background reference medium (constant)
%          delta_V          (N x 1) vector of all subvolume sizes [m^3]
%          T_vector         (N x 1) vector of all subvolume temperatures [K]
%          ind_bulk         Indices of first subvolume in a given bulk object
%
% OUTPUTS: G_sys_2D         System Green's function matrix for all subvolume interactions in 2D matrix format (3N x 3N)
%          Trans            Transmission coefficient matrix (3N x 3N)
%          Q_omega_bulk     Spectral heat dissipation within bulk objects [W/m]
%          Q_omega_subvol   Spectral heat dissipation within subvolumes [W/m*(m^3)]
%



% Determine total number of subvolumes
[N,~] = size(r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate background medium Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ ~, G_0_2D] = G_0_function(r, omega, delta_V, epsilon_ref, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate deterministic interaction matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ A_2D, ~ ] = A_matrix_function(G_0_2D, omega, epsilon, epsilon_ref, delta_V, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = toc;
G_sys_2D = A_2D\G_0_2D;  % System Green's function matrix
t2 = toc;


disp(['Time for matrix inversion = ' num2str(t2-t1) ' s = ' num2str((t2-t1)/60) ' minutes'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate transmission coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t3 = toc;
[ Trans ] = trans_coeff_function(G_sys_2D, omega, epsilon, delta_V, N);
t4 = toc;

disp(['Time to calculate transmission coefficient = ' num2str(t4-t3) ' s = ' num2str((t4-t3)/60) ' minutes'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate monochromatic power dissipated in bulk objects %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Q_omega_bulk, Q_omega_subvol ] = Q_function(Trans, omega, T_vector, N, ind_bulk);



