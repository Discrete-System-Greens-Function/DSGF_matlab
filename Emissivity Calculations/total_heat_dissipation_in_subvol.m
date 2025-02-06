function [Q_total_subvol, Q_density_subvol, Q_total_subvol_matrix, Q_total_bulk] = total_heat_dissipation_in_subvol(N, omega, Q_omega_subvol, delta_V_vector, r, ind_bulk)
	
% This function calculates the total heat dissipated in every subvolume and the power density for every subvolume.
	
% INPUTS:  
%                N                                Total number of subvolumes in all objects.
%                omega                       (N_omega x 1) vector of radial frequencies [rad/s]
%                Q_omega_subvol       (N_omega x N) matrix of the spectral conductance of every subvolume  [W/(rad/s)]
%                delta_V_vector            (N x 1) vector of all subvolume sizes [m^3]
%                r                                  (N x 3) matrix containing points of all cubic lattice points of thermal objects [m]
%
%
% OUTPUTS: 
%                  Q_total_subvol               Total power dissipated in each subvolume [W]
%                  Q_density_subvol           Total power density dissipated in a subvolume [W/(m^3)]
%                  Q_total_subvol_matrix    Restructure heat dissipated in a subvolume from a vector to a matrix with indices of coordinates intact

	% Total heat dissipated in a subvolume [W]
	Q_total_subvol = zeros(1,N); % Preallocate
	for ii = 1:N
	    Q_total_subvol(ii) = trapz(omega, Q_omega_subvol(:,ii));  % Heat dissipated in each subvolume [W]
	end

	% Total heat density dissipated in a subvolume [W/(m^3)]
	Q_density_subvol = Q_total_subvol'./delta_V_vector; % Heat density dissipated in each subvolume [W/(m^3)]

	% Restructure heat dissipated in a subvolume from a vector to a matrix with
	% indices of coordinates intact
	Q_total_subvol_matrix = [r, Q_total_subvol.'];
    
    % Total heat dissipated in each bulk object [W]
	Q_total_bulk = zeros(1,2); % Preallocate
    
    % Calculate total power dissipated in each bulk object [W]
for m = 1:length(ind_bulk)  % Loop through each bulk thermal object
    if m ~= length(ind_bulk)
        Q_total_bulk(m) = sum(Q_total_subvol(ind_bulk(m):(ind_bulk(m+1)-1)));
    else
        Q_total_bulk(m) = sum(Q_total_subvol(ind_bulk(m):end));
    end
end


end
