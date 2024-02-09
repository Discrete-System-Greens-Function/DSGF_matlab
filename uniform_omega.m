function [omega] = uniform_omega(start_omega, end_omega, N_omega)
% calculates a uniformly discretized vector of angular frequencies for a given range of vacuum wavelengths
% and how many elements it should have
% 
%	Inputs:
%		start_lambda - beginning wavelength [m]
%		end_lambda - end wavelength [m]
%		N_omega - number of evenly spaced angular frequencies
%
%	Outputs:
%		omega - angular frequency vector [rad/s]
%

	omega_pre = linspace(start_omega, end_omega, N_omega); % [m]
    
    omega = omega_pre';
    
end
