function [omega] = even_omega(start_lambda, end_lambda, N_omega)
% calculates an even omega vector based on the beginning and ending wavelength
% and how many elements it should have
% 
%	Inputs:
%		start_lambda - beginning wavelength
%		end_lambda - end wavelength
%		N_omega - number of evenly spaced angular velcoties
%
%	Outputs:
%		omega - angular frequency vector
%

	lambda = linspace(start_lambda, end_lambda, N_omega); % [m]
	c_0 = 299792458;            % Speed of light in vacuum [m/s]
	omega = (2*pi*c_0./lambda).'; % [rad/s]

end
