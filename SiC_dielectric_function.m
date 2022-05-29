function [ epsilon ] = SiC_dielectric_function( omega , constants)

% This code uses the Drude-Lorentz model to calculate the dielectric
% function of silicon carbide (SiC).
% See Francoeur et al. 2011.

% INPUTS:  omega         (N_omega x 1) vector of frequencies at which to calculate the dielectric function [rad/s]
%	   constants	 struct containing all the constants needed
%
%
% OUTPUTS: epsilon       dielectric function [-]

% Wave energy
E_joules = constants.h_bar*omega; % [J]
E_eV = E_joules/constants.q;      % [eV]

% Dielectric function approximated with a Lorentz model
omega_TO = 1.494e14;              % Transverse optical phonon frequency [rad/s]
omega_LO = 1.825e14;              % Longitudinal optical phonon frequency [rad/s]
Gamma_D = 8.966e11;               % Damping constant [rad/s]
epsilon_inf = 6.7;                % High-frequency limit of permittivity
epsilon = epsilon_inf*( ((omega.^2) - (omega_LO^2) + 1i*Gamma_D*omega)./ ((omega.^2) - (omega_TO^2) + 1i*Gamma_D*omega) );

% Complex refractive indices
%m = sqrt(epsilon);
