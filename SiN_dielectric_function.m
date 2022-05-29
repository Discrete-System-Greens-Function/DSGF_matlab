function [ epsilon_r ] = SiN_dielectric_function( omega, constants)

% This code uses the Lorentz oscillator model to calculate the dielectric
% function of silicon nitride (SiN). The dielectric function Lorentz 
% oscillator terms come from G. Cataldo, J. A. Beall, H.-M. Cho, B. 
% McAndrew, M. D. Niemack, and E. J. Wollack, “Infrared dielectric 
% properties of low-stress silicon nitride,” Opt. Lett., vol. 37, no. 20,
% pp. 4200–4202, 2012.  These are the same values used and reported in the
% Supplementary Materials of Y. Wu et al., “Enhanced thermal conduction by 
% surface phonon-polaritons,” Sci. Adv., vol. 6, p. eabb4461, 2020.

% INPUTS:  omega         (N_omega x 1) vector of frequencies at which to calculate the dielectric function [rad/s]
%	   constants     struct containing all the necessary constants
%
%
% OUTPUTS: epsilon_r     dielectric function [-]


M = 5; % Number of oscillations
epsilon_real = [7.582, 6.754, 6.601, 5.430, 4.601, 4.562];
epsilon_imag = [0, 0.3759, 0.0041, 0.1179, 0.2073, 0.0124];
epsilon = epsilon_real + 1i.*epsilon_imag;
epsilon_inf = epsilon(M+1);
omega_T = (2*pi).*(1e12).*[13.913, 15.053, 24.521, 26.440, 31.724];
Gamma = (2*pi).*(1e12).*[5.810, 6.436, 2.751, 3.482, 5.948];
alpha = [0.0001, 0.3427, 0.0006, 0.0002, 0.0080];

summation = 0;

for jj = 1:M
    
    num_1 = ((omega_T(jj))^2) - omega.^2;
    denom_1 = omega.*Gamma(jj);
    Gamma_prime = Gamma(jj).*exp(-alpha(jj).*((num_1./denom_1).^2));
    
    delta_epsilon = epsilon(jj) - epsilon(jj+1);
    
    num_2 = delta_epsilon*((omega_T(jj))^2);
    denom_2 = ((omega_T(jj))^2) - omega.^2 - 1i.*omega.*Gamma_prime;
    
    summation = summation + num_2./denom_2;
    
end

epsilon_r = epsilon_inf + summation;


% % Wave energy
% E_joules = constants.h_bar*omega; % [J]
% E_eV = E_joules/constants.q;      % [eV]
% 
% % Dielectric function approximated with a Lorentz model
% omega_TO = 93e12;              % Transverse optical phonon frequency [rad/s]
% omega_LO = 96e12;              % Longitudinal optical phonon frequency [rad/s]
% Gamma_TO = 149e12;               % Damping constant [rad/s]
% Gamma_LO = 214e12;               % Damping constant [rad/s]
% 
% epsilon_inf = 3.3;                % High-frequency limit of permittivity (I just guessed this value by looking at the graph.)
% 
% num = omega_LO.^2 - omega.^2 - 1i.*Gamma_LO.*omega;
% denom = omega_TO.^2 - omega.^2 - 1i.*Gamma_TO.*omega;
% epsilon = epsilon_inf + (num./denom);


