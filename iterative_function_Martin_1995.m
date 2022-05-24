function [ G_sys_2D, Trans, Q_omega_bulk, Q_omega_subvol ] = iterative_function_Martin_1995( omega, r, epsilon, epsilon_ref, delta_V, T_vector, ind_bulk )

% This function calculates the system Green's function using the
% iterative technique outlined in the paper O. J. F. Martin, C. Girard,
% and A. Dereux, Generalized Field Propagator for Electromagnetic
% Scattering and Light Confinement, Phys. Rev. Lett. 74, 526 (1995).

% INPUTS:  omega          radial frequency [rad/s]
%          r              (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%          epsilon        (N x 1) vector of all subvolume dielectric functions
%          epsilon_ref    dielectric function for background reference medium (constant)
%          delta_V        (N x 1) vector of all subvolume sizes [m^3]
%          T_vector       (N x 1) vector of all subvolume temperatures [K]
%          ind_bulk       Indices of first subvolume in a given bulk object
%
% OUTPUTS: G_0_4D         free-space Green's function matrix for all subvolume interactions in 4D matrix format (N x N x 3 x 3)
%          G_0_2D         free-space Green's function matrix for all subvolume interactions in 2D matrix format (3N x 3N)



% Determine total number of subvolumes
[N,~] = size(r);

% Constants
epsilon_0 = 8.8542e-12;     % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]

% Wave vector in background reference medium
k = omega*sqrt(epsilon_ref*epsilon_0*mu_0);

% 3-by-3 unit matrix
I = eye(3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate background medium Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ ~, G_0_2D] = G_0_function(r, omega, delta_V, epsilon_ref, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate system Green's function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scattering dielectric function
epsilon_s = epsilon - epsilon_ref;

% Preallocate and set initial values
G_sys_2D_new = zeros(3*N, 3*N);

% Set initial values to free-space Green's function values
G_sys_2D_old = G_0_2D;

for mm = 1:N
    ind_mm = (3*mm - 2):(3*mm); % Set indices

   

    % First, solve ii = mm system of equations.
    A_mm = I - (k^2).*epsilon_s(mm).*delta_V(mm).*G_sys_2D_old(ind_mm, ind_mm);
    for jj = 1:N % Only loop through remaining perturbations
        ind_jj = (3*jj - 2):(3*jj); % Set indices

        G_sys_2D_new(ind_mm, ind_jj) = A_mm\G_sys_2D_old(ind_mm, ind_jj);

    end % End loop through jj

    % Next, solve all systems of equations for ii not equal to mm.
    for ii = 1:N
        ind_ii = (3*ii - 2):(3*ii); % Set indices

        if ii ~= mm % Only loop through ii is not equal to mm
            for jj = 1:N % Only loop through remaining perturbations
                ind_jj = (3*jj - 2):(3*jj); % Set indices

                G_sys_2D_new(ind_ii, ind_jj) = G_sys_2D_old(ind_ii, ind_jj) + (k^2)*epsilon_s(mm)*delta_V(mm)*G_sys_2D_old(ind_ii, ind_mm)*G_sys_2D_new(ind_mm, ind_jj);

            end % End loop through jj
        end % End ii ~= mm check
    end % End loop through ii

    % Update matrices
    G_sys_2D_old = G_sys_2D_new;

end % End loop through mm

G_sys_2D = G_sys_2D_new;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate transmission coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Trans ] = trans_coeff_function(G_sys_2D, omega, epsilon, delta_V, N);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate monochromatic power dissipated in bulk objects %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Q_omega_bulk, Q_omega_subvol ] = Q_function(Trans, omega, T_vector, N, ind_bulk);



