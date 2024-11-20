function [G_sys_2D] = G_sys_function(omega, r, epsilon, epsilon_ref, delta_V)

% This function calculates the system Green's function

% INPUTS:  omega            Radial frequency [rad/s]
%          r                (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%          epsilon          (N x 1) vector of all subvolume dielectric functions
%          epsilon_ref      Dielectric function for background reference medium (constant)
%          delta_V          (N x 1) vector of all subvolume sizes [m^3]
%
% OUTPUTS: G_sys_2D         System Green's function matrix for all subvolume interactions in 2D matrix format (3N x 3N)



% Determine total number of subvolumes
[N,~] = size(r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate background medium Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%original 
[ ~, G_0_2D] = G_0_free_space(r, omega, delta_V, epsilon_ref, N);%free space greens 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populate deterministic interaction matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ A_2D, ~ ] = A_matrix_function(G_0_2D, omega, epsilon, epsilon_ref, delta_V, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%original 
t1 = toc;
%G_sys_2D = Invers(A_2D).*G_0_2D;
G_sys_2D = A_2D\G_0_2D;  % System Green's function matrix
t2 = toc;





