function [ G_omega_bulk, G_bulk ] = conductance_bulk( Trans_omega_bulk, T, omega )

% This function calculates the radiative conductance (both spectral and 
% integrated values) between two bulk objects given the spectral 
% transmission coefficient between the two bulk objects.

% INPUTS:  Trans_omega_bulk      (N_omega x 1) Spectral transmission coefficient vector
%                                 Note: Trans_bulk(ind_1 --> ind_2)
%                                 [unitless]
%          T                      Temperature at which to evaluate conductance [K]
%          omega                 (N_omega x 1) Vector of frequencies [rad/s]
%
%
%
% OUTPUTS: G_omega_bulk    (N_omega x 1) Spectral radiative conductance between bulk objects 1 and 2
%                           Note: G_omega_bulk(ind_1 --> ind_2)
%          G_bulk           Integrated radiative conductance between bulk objects 1 and 2
%                           Note: G_bulk(ind_1 --> ind_2)


% Constants
h_bar = 1.0545718e-34;      % Planck's constant [J*s]
k_b = 1.38064852e-23;       % Boltzmann constant [J/K]

% Derivative of mean energy of an electromagnetic state with respect to temperature
num = ((h_bar.*omega).^2).*exp((h_bar.*omega)./(k_b.*T));
denom = k_b.*(T.^2).*((exp((h_bar.*omega)./(k_b.*T)) - 1).^2);
dtheta_dT = num./denom;

% Spectral radiative conductance between bulk objects
G_omega_bulk = dtheta_dT.*Trans_omega_bulk;

% Integrate for total radiative conductance between bulk objects
if omega(1) > omega(2)
    G_bulk = (1/(2*pi)).*trapz(flipud(omega), flipud(G_omega_bulk));
else
    G_bulk = (1/(2*pi)).*trapz(omega, G_omega_bulk);
end
