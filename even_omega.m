function [omega] = even_omega(start_lambda, end_lambda, N_omega)
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

	lambda = linspace(start_lambda, end_lambda, N_omega); % [m]
	c_0 = 299792458;            % Speed of light in vacuum [m/s]
	omega = (2*pi*c_0./lambda).'; % [rad/s]

end
