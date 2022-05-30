function [Q_total_subvol, Q_density_subvol, Q_total_subvol_matrix] = total_heat_dissipation_in_subvol(N, omega, Q_omega_subvol, delta_V_vector, r)
	

	% Total heat dissipated in a subvolume [W]
	Q_total_subvol = zeros(1,N); % Preallocate
	for ii = 1:N
	    Q_total_subvol(ii) = trapz(omega, Q_omega_subvol(:,ii));  % Heat dissipated in each subvolume [W]
	end

	% Total heat density dissipated in a subvolume [W/(m^3)]
	Q_density_subvol = Q_total_subvol./delta_V_vector; % Heat density dissipated in each subvolume [W/(m^3)]

	% Restructure heat dissipated in a subvolume from a vector to a matrix with
	% indices of coordinates intact
	Q_total_subvol_matrix = [r, Q_total_subvol.'];

end
