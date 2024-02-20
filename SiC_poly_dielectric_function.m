function [ epsilon ] = SiC_poly_dielectric_function( omega )

% This code uses the Drude-Lorentz model to calculate the dielectric
% function of silicon carbide (SiC).

% INPUTS:  omega         vector of frequencies at which to calculate the dielectric function [rad/s]
%
%
% OUTPUTS: epsilon       dielectric function [-]


% Dielectric function approximated with a Lorentz model
omega_TO = 1.486e14; %1.494e14;              % Transverse optical phonon frequency [rad/s]
omega_LO = 1.801e14; %1.825e14;              % Longitudinal optical phonon frequency [rad/s]
Gamma_D = 3.767e12; %8.966e11;               % Damping constant [rad/s]
epsilon_inf = 8; %6.7;                % High-frequency limit of permittivity
epsilon = epsilon_inf*( ((omega.^2) - (omega_LO^2) + 1i*Gamma_D*omega)./ ((omega.^2) - (omega_TO^2) + 1i*Gamma_D*omega) );

% Complex refractive indices
%m = sqrt(epsilon);
