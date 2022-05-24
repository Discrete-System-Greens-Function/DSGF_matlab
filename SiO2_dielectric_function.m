function [ epsilon ] = SiO2_dielectric_function( omega )

% This code uses the Lorentz oscillator model to calculate the dielectric
% function of silicon dioxide (SiO2).
%
% Values taken from B. Czapla and A. Narayanaswamy, “Thermal radiative
% energy exchange between a closely-spaced linear chain of spheres and its
% environment,” J. Quant. Spectrosc. Radiat. Transf., vol. 227, pp. 4–11, 
% 2019.

% INPUTS:  omega         (N_omega x 1) vector of frequencies at which to calculate the dielectric function [rad/s]
%
%
% OUTPUTS: epsilon       dielectric function [-]


% Constants
q = 1.60218e-19;            % Number of joules per eV [J/eV]
h_bar = 1.0545718e-34;      % Planck's constant [J*s]
% k_b = 1.38064852e-23;       % Boltzmann constant [J/K]
% epsilon_0 = 8.8542e-12;     % Permittivity of free space [F/m]
% mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]

% % Wave energy
% E_joules = h_bar*omega; % [J]
% E_eV = E_joules/q;      % [eV]

% Dielectric function approximated with a Lorentz model
epsilon_inf = 2.03843;                          % High-frequency limit of permittivity
h_bar_omega_0 = [0.05624, 0.09952, 0.13355];    % [eV]
omega_0 = (h_bar_omega_0*q)/(h_bar);            % [rad/s]
%lambda_0 = [22.04432, 12.45818, 9.28364]*1e-6;  % [m]
S = [0.93752, 0.05050, 0.60642];                % [-]
Gamma = [0.09906, 0.05511, 0.05246];            % [-]

summation = zeros(length(omega),1); % Preallocate vector
for k = 1:length(S)
    num = S(k);
    term_1 = (omega./omega_0(k)).^2;
    term_2 = 1i*Gamma(k).*(omega./omega_0(k));
    denom = 1 - term_1 - term_2;
    summation = summation + (num./denom);
end
epsilon = epsilon_inf + summation;


% omega_TO = 1.494e14;              % Transverse optical phonon frequency [rad/s]
% omega_LO = 1.825e14;              % Longitudinal optical phonon frequency [rad/s]
% Gamma_D = 8.966e11;               % Damping constant [rad/s]
% epsilon_inf = 6.7;                % High-frequency limit of permittivity
% epsilon = epsilon_inf*( ((omega.^2) - (omega_LO^2) + 1i*Gamma_D*omega)./ ((omega.^2) - (omega_TO^2) + 1i*Gamma_D*omega) );

% Complex refractive indices
%m = sqrt(epsilon);
