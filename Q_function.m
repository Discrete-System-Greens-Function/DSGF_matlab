function [ Q_omega_bulk, Q_omega_subvol ] = Q_function( Trans, omega, T_vector, N, ind_bulk )

% This function calculates the monochromatic thermal power dissipated in 
% each subvolume and in each bulk object.

% INPUTS:  Trans            (N x N) transmission coefficient matrix
%          omega            Radial frequency [rad/s] (a single scalar value)
%          T_vector         (N x 1) vector of all subvolume temperatures [K]
%          N                Total number of subvolumes (i.e. dipoles)
%          ind_bulk         Indices of first subvolume in a given bulk object
%
%
% OUTPUTS: Q_omega_bulk     Total power dissipated in each subvolume [W]
%          Q_omega_subvol   Total power dissipated in each bulk object [W]     


% Constants
h_bar = 1.0545718e-34;      % Planck's constant [J*s]
k_b = 1.38064852e-23;       % Boltzmann constant [J/K]


% Mean energy of an electromagnetic state for all subvolumes
theta = (h_bar*omega)./(exp(h_bar*omega./(k_b*T_vector)) - 1);

% Preallocate vectors
Q_omega_subvol = zeros(1,N);
inner_sum = zeros(1,N);
Q_omega_bulk = zeros(1,length(ind_bulk));

% % Calculate total power dissipated in each subvolume [W]
% for i = 1:N       % Loop through all subvolumes
%     for j = 1:N   % Loop through all subvolumes
%         
%         % Using Eric's Eq. 26
%         if j ~= i
%             inner_sum(j) = (theta(j) - theta(i))*Trans(i,j);
%         else
%             inner_sum(j) = 0;
%         end
%     end
%     Q_omega_subvol(i) = (1/(2*pi))*sum(inner_sum);
% end

% Vectorized version
theta_matrix = repmat(theta, 1, N);
theta_ji = theta_matrix - (theta_matrix.');
Q_omega_subvol = (1/(2*pi))*theta_ji.*Trans;
Q_omega_subvol = sum(Q_omega_subvol);

% Calculate total power dissipated in each bulk object [W]
for m = 1:length(ind_bulk)  % Loop through each bulk thermal object
    if m ~= length(ind_bulk)
        Q_omega_bulk(m) = sum(Q_omega_subvol(ind_bulk(m):(ind_bulk(m+1)-1)));
    else
        Q_omega_bulk(m) = sum(Q_omega_subvol(ind_bulk(m):end));
    end
end



