function [ check_1, check_2, check_3, check_4, check_5,...
    omega_failed_2, N_failed_2, omega_failed_3, N_failed_3,...
    omega_failed_4, N_failed_4, omega_failed_5, N_failed_5 ] ...
    = convergence_check_function( delta_V_vector, d, omega, epsilon, epsilon_ref, tol_1, tol_2, tol_3, tol_4, tol_5 )

% This function performs convergence checks for a given discretization and
% and material when used with the SGF-TDDA method.  Criteria are taken from
% S. Edalatpour, M. ?uma, T. Trueax, R. Backman, and M. Francoeur, 
% “Convergence analysis of the thermal discrete dipole approximation,” 
% Phys. Rev. E, vol. 91, no. 6, p. 063307, 2015.

% INPUTS:  delta_V_vector    (N x 1) vector of all subvolume sizes [m^3]
%          d                 Closest surface-to-surface distance between bulk objects
%          omega             (N_omega x 1) vector of all angular frequencies [rad/s]
%          epsilon           (N_omega x 1) vector of dielectric function of thermal objects [-]
%          epsilon_ref       Dielectric function of background reference medium [-]
%          tol_1             Tolerated minimum ratio of d/[delta_V^(1/3)]
%                            (__2_ is usually an acceptable number)
%          tol_2             Tolerated minimum ratio of (1/k_therm)/[delta_V^(1/3)]
%                            (__10_ is usually an acceptable number)
%          tol_3             Tolerated minimum ratio of [lambda_0/Im{m}]/[delta_V^(1/3)]
%                            (__2_ is usually an acceptable number)
%          tol_4             Tolerated minimum ratio of lambda_0/[delta_V^(1/3)]
%                            (__2_ is usually an acceptable number)
%          tol_5             Tolerated minimum ratio of [lambda_0/Im{m}]/[delta_V^(1/3)]
%                            (__2_ is usually an acceptable number)
%                            (Constant free space Green's function for each subvolume.)
%
% 
% OUTPUTS: check_1           "Passed"/"Failed" designation for d/[delta_V^(1/3)] >> 1 convergence criteria.
%          check_2           "Passed"/"Failed" designation for (1/k_therm)/[delta_V^(1/3)] >> 1 convergence criteria.
%          check_3           "Passed"/"Failed" designation for [lambda_0/Im{m}]/[delta_V^(1/3)] >> 1 convergence criteria.
%          check_4           "Passed"/"Failed" designation for [delta_V^(1/3)]/lambda_0 >> 1 convergence criteria.
%          omega_failed_2    Frequencies at which (1/k_therm)/[delta_V^(1/3)] >> 1 convergence criteria failed.
%          N_failed_2        Number of frequencies for which (1/k_therm)/[delta_V^(1/3)] >> 1 convergence criteria failed.
%          omega_failed_3    Frequencies at which [lambda_0/Im{m}]/[delta_V^(1/3)] >> 1 convergence criteria failed.
%          N_failed_3        Number of frequencies for which [lambda_0/Im{m}]/[delta_V^(1/3)] >> 1 convergence criteria failed.
%          omega_failed_4    Frequencies at which lambda_0/[delta_V^(1/3)] >> 1 convergence criteria failed.
%          N_failed_4        Number of frequencies for which lambda_0/[delta_V^(1/3)] >> 1 convergence criteria failed.
%          omega_failed_5    Frequencies at which [delta_V^(1/3)]/lambda_0 >> 1 convergence criteria failed.
%          N_failed_5        Number of frequencies for which [delta_V^(1/3)]/lambda_0 >> 1 convergence criteria failed.
%


% Constants
epsilon_0 = 8.8542e-12;     % Permittivity of free space [F/m]
mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]
c_0 = 299792458;            % Speed of light in vacuum [m/s]

% Number of frequencies
N_omega = length(omega);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Criteria #1: Subvolume length scale must be much smaller than distance
%              between closest surfaces of objects.
%              ( delta_V^(1/3) << d )
char_length = delta_V_vector.^(1/3);
check_1 = find((d./char_length) < tol_1);
if isempty(check_1) == 0
    disp(['Failed convergence criteria #1 (d/[delta_V^(1/3)] > ' num2str(tol_1) ').']);
    check_1 = 'Failed';
else
    check_1 = 'Passed';
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Criteria #2: Subvolume length scale must be much smaller than one over the
% %              wavenumber in the thermal objects.
% %              ( delta_V^(1/3) << 1/k )
% count_2 = 0;
% k_therm = abs(omega.*sqrt(epsilon.*mu_0*epsilon_0));
% omega_failed_2 = zeros(1,N_omega);
% 
% for omega_loop = 1:N_omega % Loop through all frequencies
%     check_2 = find((k_therm(omega_loop).*char_length) > tol_2);
%     if isempty(check_2) == 0
%         disp(['Failed convergence criteria #2 ([1/k_therm]/[delta_V^(1/3)] > ' num2str(tol_2) ') for omega = ' num2str(omega(omega_loop)) ' rad/s.']);
%         count_2 = count_2 + 1;
%         check_2 = 'Failed';
%         omega_failed_2(omega_loop) = omega(omega_loop);
%     end
% end
% 
% if sum(omega_failed_2) == 0
%     check_2 = 'Passed';
% end
% 
% N_failed_2 = count_2; % Number of frequencies for which Criteria #2 failed.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Criteria #2: Subvolume length scale must be much smaller than the material
%              wavelength.
%              ( delta_V^(1/3) << lambda_m )

% Calculate index of refraction
m_real = sqrt(0.5*(sqrt((real(epsilon)).^2 + (imag(epsilon)).^2) + real(epsilon)));
m_imag = sqrt(0.5*(sqrt((real(epsilon)).^2 + (imag(epsilon)).^2) - real(epsilon)));
m = m_real + 1i.*m_imag; % Complex index of refraction

% Free space wavevector
k_0 = omega.*sqrt(epsilon_0*mu_0);

% Free space wavelength
lambda_0 = (2*pi)./k_0;

% Material wavelencth
lambda_m = lambda_0./m_real;

count_2 = 0;
omega_failed_2 = zeros(1,N_omega);

for omega_loop = 1:N_omega % Loop through all frequencies
    check_2 = find(lambda_m(omega_loop)./char_length < tol_2);
    if isempty(check_2) == 0
        disp(['Failed convergence criteria #2 ([lambda_0/Re{m}]/[delta_V^(1/3)] > ' num2str(tol_2) ') for omega = ' num2str(omega(omega_loop)) ' rad/s.']);
        count_2 = count_2 + 1;
        omega_failed_2(omega_loop) = omega(omega_loop);
    end
end

N_failed_2 = count_2; % Number of frequencies for which Criteria #2 failed.

if N_failed_2 == 0
    check_2 = 'Passed';
else
    check_2 = 'Failed';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Criteria #3: Subvolume length scale must be smaller than the decay
%              length of the electric field.
%              ( delta_V^(1/3) << lambda_0/Im{m} )

% Calculate index of refraction
m_real = sqrt(0.5*(sqrt((real(epsilon)).^2 + (imag(epsilon)).^2) + real(epsilon)));
m_imag = sqrt(0.5*(sqrt((real(epsilon)).^2 + (imag(epsilon)).^2) - real(epsilon)));
m = m_real + 1i.*m_imag; % Complex index of refraction

% Free space wavevector
k_0 = omega.*sqrt(epsilon_0*mu_0);

% Free space wavelength
lambda_0 = (2*pi)./k_0;

% Preallocate
count_3 = 0;
omega_failed_3 = zeros(1,N_omega);

for omega_loop = 1:N_omega % Loop through all frequencies
    check_3 = find((lambda_0(omega_loop)/m_imag(omega_loop))./char_length < tol_3);
    if isempty(check_3) == 0
        disp(['Failed convergence criteria #3 ([lambda_0/Im{m}]/[delta_V^(1/3)] > ' num2str(tol_3) ') for omega = ' num2str(omega(omega_loop)) ' rad/s.']);
        count_3 = count_3 + 1;
        omega_failed_3(omega_loop) = omega(omega_loop);
    end
end

N_failed_3 = count_3; % Number of frequencies for which Criteria #3 failed.

if N_failed_3 == 0
    check_3 = 'Passed';
else
    check_3 = 'Failed';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Criteria #4: Subvolume length scale must be much smaller than the free
%              space wavelength.
%              ( lambda_0 >> delta_V^(1/3))

% Preallocate
count_4 = 0;
omega_failed_4 = zeros(1,N_omega);

for omega_loop = 1:N_omega % Loop through all frequencies
    check_4 = find(lambda_0(omega_loop)./char_length < tol_4);
    if isempty(check_4) == 0
        disp(['Failed convergence criteria #4 (lambda_0/[delta_V^(1/3)] > ' num2str(tol_4) ') for omega = ' num2str(omega(omega_loop)) ' rad/s.']);
        count_4 = count_4 + 1;
        omega_failed_4(omega_loop) = omega(omega_loop);
    end
end

N_failed_4 = count_4; % Number of frequencies for which Criteria #4 failed.

if N_failed_4 == 0
    check_4 = 'Passed';
else
    check_4 = 'Failed';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Criteria #5: Constant free-space Green's function inside each subvolume.
%              The ratio of free space wavelegth divided by the distance
%              between cubic lattice points must be small.
%              ( lambda_0/(r-r') << 1)
% THIS CRITERIA SEEMS TO BE IN CONTRADICTION TO CRITERIA #4!

% Preallocate
count_5 = 0;
omega_failed_5 = zeros(1,N_omega);

for omega_loop = 1:N_omega % Loop through all frequencies
    check_5 = find(char_length/lambda_0(omega_loop) > tol_5);
    if isempty(check_5) == 0
        disp(['Failed convergence criteria #4 ([delta_V^(1/3)]/lambda_0 > ' num2str(tol_5) ') for omega = ' num2str(omega(omega_loop)) ' rad/s.']);
        count_5 = count_5 + 1;
        omega_failed_5(omega_loop) = omega(omega_loop);
    end
end

N_failed_5 = count_5; % Number of frequencies for which Criteria #4 failed.

if N_failed_5 == 0
    check_5 = 'Passed';
else
    check_5 = 'Failed';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out results

disp(' ')
disp(['Convergence criteria #1 (d/[delta_V^(1/3)] > ' num2str(tol_1) '): ' check_1]);
% disp(['Convergence criteria #2 ([1/k_therm]/[delta_V^(1/3)] > ' num2str(tol_2) '): ' check_2]);
disp(['Convergence criteria #2 ([lambda_0/Re{m}]/[delta_V^(1/3)] > ' num2str(tol_2) '): ' check_2]);
disp(['Convergence criteria #3 ([lambda_0/Im{m}]/[delta_V^(1/3)] > ' num2str(tol_3) '): ' check_3]);
disp(['Convergence criteria #4 (lambda_0/[delta_V^(1/3)] > ' num2str(tol_4) '): ' check_4]);
disp(['Convergence criteria #5 ([delta_V^(1/3)]/lambda_0 > ' num2str(tol_5) '): ' check_5]);
disp(' ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots

k_0 = omega.*sqrt(epsilon_0*mu_0);
lambda_0 = 2*pi./k_0;

% lambda_ref = lambda_0./sqrt(epsilon_ref);
% epsilon_real = real(epsilon);
% lambda_mat = lambda_0./sqrt(epsilon_real);
% lambda_decay = lambda_0./sqrt(imag(epsilon));

lambda_ref = lambda_0./sqrt(epsilon_ref);
lambda_mat = lambda_0./m_real;
lambda_decay = lambda_0./m_imag;

L_sub = delta_V_vector(1).^(1/3);

ratio_ref = L_sub./lambda_ref;
ratio_mat = L_sub./lambda_mat;
ratio_decay = L_sub./lambda_decay;
%ratio_pen = L_sub./d;

FIG_check_convergence = figure(11);
plot(omega, ratio_ref, omega, ratio_mat, omega, ratio_decay, 'linewidth', 4)
xlabel('Frequency, \omega [rad/s]', 'fontsize', 12)
ylabel('Convergence ratios, L_s_u_b/\lambda', 'fontsize', 12)
title(['Convergence analysis, L_s_u_b = ' num2str(L_sub*(1e9)) ' nm'], 'fontsize', 16)
set(gca, 'fontsize', 22)
legend('Ratio for wavelength in background medium, L_s_u_b/\lambda_r_e_f',...
    'Ratio for wavelength in material, L_s_u_b/\lambda_m',...
    'Ratio for decay length in material, L_s_u_b/\lambda_d_e_c_a_y', 'location', 'best')
axis tight
grid on









