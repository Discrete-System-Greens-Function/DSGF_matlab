%function [] = DSGF_main(description, discretization, material, T, T_conductance, epsilon_ref, omega, wave_type, output, discretization_type, L_char, origin, delta_V)
function [] = DSGF_main(description, discretization, material, T, T_conductance, epsilon_ref, omega, wave_type, output, discretization_type, L_char, delta_V, d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function implements the DSGF method to model near-field radiative
% heat transfer between 3D objects.
%
% INPUTS:   discretization
%           L_char
%           origin
%           material
%           T
%           T_conductance
%           epsilon_ref
%           omega
%           observation_point
%           output
%
%
%
% OUTPUTS: Saved figures and workspace variables.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
disp(['Running MATLAB script ' mfilename])
tic


%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%

constants = struct();

constants.q = 1.60218e-19;            % Number of joules per eV [J/eV]
constants.h_bar = 1.0545718e-34;      % Planck's constant [J*s]
constants.k_b = 1.38064852e-23;       % Boltzmann constant [J/K]
constants.epsilon_0 = 8.8542e-12;     % Permittivity of free space [F/m]
constants.mu_0 = (4*pi)*(10^-7);      % Permeability of free space [H/m]
constants.c_0 = 299792458;            % Speed of light in vacuum [m/s]


%%%%%%%%%%%%%%%%%%%%%%
% Set results export %
%%%%%%%%%%%%%%%%%%%%%%

[filePath_st, file_name_saved] = result_setup(output, description);

%%%%%%%%%%%%%%%%%%
% Figure options %
%%%%%%%%%%%%%%%%%%

% Show figure axes?
show_axes = true; %false


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set calculation options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation approach
% 'direct' = direct matrix inversion using 'mldivide' operator
% 'iterative' =  iterative solver from Martin et al.
calc_approach = CalculationOption.direct;

% Conduct convergence analysis?
convergence_analysis = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************DISCRETIZATION*******************************%



if strcmp('sample',discretization_type)
    
     %[N_each_object, volume, r_each_object, ind_bulk, delta_V_each_object, L_sub_each_object] = read_discretization(discretization, L_char, origin);
     [N_each_object, volume, r_each_object, ind_bulk, delta_V_each_object, L_sub_each_object,origin] = read_discretization(discretization, L_char,d);
     
     % Discretized lattice including subvolumes of all objects in one matrix (N x 3 matrix)
     r = cell2mat(r_each_object);
    
     % Total number of subvolumes
     [N,~] = size(r);
    
     r_1 = r(1:ind_bulk(2)-1,:);
     r_2 = r(ind_bulk(2):N,:);
     
    % Subvolume size for all N subvolumes (N x 1 vectors)
    delta_V_vector = cell2mat(delta_V_each_object); % Volume of subvolumes for all N subvolumes
    L_sub_vector = cell2mat(L_sub_each_object);     % Length of side of a cubic subvolume for all N subvolumes
     
    % UPDATE FOR MORE THAN 2 OBJECTS!
    % Center-of-mass separation distance [m]
    d_center = norm(origin(1,:) - origin(2,:)); 

    % UPDATE FOR MORE THAN 2 OBJECTS!
    % Closest vacuum gap separation distance [m]
    [ d_min_center, d_min_edge, r_1_min, r_2_min ] = calculate_surface_separation( r_each_object{1}, r_each_object{2}, L_sub_each_object{1}(1),  L_sub_each_object{2}(1));

    % UPDATE FOR MORE THAN 2 OBJECTS!
    % Temperature
    T_vector = [T(1).*ones(N_each_object(1),1); T(2).*ones(N_each_object(2),1)];
    
elseif strcmp('user_defined',discretization_type)
        
    discFile = discretization; % File name of discretization
    discDir = "Library/Discretizations/User_defined"; % Directory where discretization is stored

    % Import discretization 
    r = readmatrix(append(append(append(discDir, '/'), discFile), '.txt'));

    % Total number of subvolumes
    [N,~] = size(r);
    N1 = N/2;       % Preallocate
    N2 = N/2;       % Preallocate
    ind_bulk = [1,N1+1]; 

    % Subvolume size for all N subvolumes (N x 1 vectors)
    discFile2 = delta_V;
    delta_V_vector = readmatrix(append(append(append(discDir, '/'), discFile2), '.txt')); % Volume of subvolumes for all N subvolumes
 
    L_sub_vector = delta_V_vector.^(1/3);
    
    % Temperature
    T_vector = [T(1).*ones(N1,1); T(2).*ones(N2,1)];
end

%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%

if strcmp('sample',discretization_type)
    discretization_plotting(r, L_sub_vector, N, show_axes, output, filePath_st.main,ind_bulk,r_1,r_2);
elseif strcmp('user_defined',discretization_type)
    discretization_plotting_user_defined(r, L_sub_vector, N, show_axes, output, filePath_st.main, N1);
end

% % String name describing simulation geometry
% simulation_geometry = [description '_R' num2str((1e9)*radius_1) 'nm_R' num2str((1e9)*radius_2) 'nm_dc' num2str((1e9)*d_center) 'nm_N' num2str(N1 + N2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*********************CALCULATE DIELECTRIC FUNCTION***********************%

% Number of discretized frequencies
N_omega = length(omega);

%%%%%%%%%%%%%%%%%%%%%%%%
% Dielectric Functions %
%%%%%%%%%%%%%%%%%%%%%%%%

switch (material)
	case Material.SiO2

	    epsilon = SiO2_dielectric_function(omega, constants); % (N x 1) vector of all dielectric functions for every frequency

	case Material.SiC

	    %epsilon = SiC_dielectric_function(omega, constants); % (N x 1) vector of all dielectric functions for every frequency

 	    epsilon = SiC_poly_dielectric_function(omega); % (N x 1) vector of all dielectric functions for every frequency

	case Material.SiN

	    epsilon = SiN_dielectric_function(omega, constants); % (N x 1) vector of all dielectric functions for every frequency

	case Material.Si3N4

	    epsilon = SiN_dielectric_function(omega, constants); % (N x 1) vector of all dielectric functions for every frequency
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dielectric function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dielectric_function_plotting(omega, epsilon, material, N_omega, output, filePath_st.main);

%****************END CALCULATION OF DIELECTRIC FUNCTION*******************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***************************CHECK CONVERGENCE*****************************%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check convergence criteria before proceeding %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if convergence_analysis
    % Set tolerances
    tol_1 = 1;
    tol_2 = 6;
    tol_3 = 6;
    tol_4 = 6;
    tol_5 = 6;

    [ check_1, check_2, check_3, check_4, check_5,...
        omega_failed_2, N_failed_2, omega_failed_3, N_failed_3,...
        omega_failed_4, N_failed_4, omega_failed_5, N_failed_5 ] ...
        = convergence_check_function( delta_V_vector, d_min_edge, omega, epsilon, epsilon_ref, tol_1, tol_2, tol_3, tol_4, tol_5 );
end % End convergence check

%*************************END CONVERGENCE CHECK***************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate DSGF and power dissipation %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate vectors
Trans_omega_12 = zeros(N_omega,1);
%Trans_omega_13 = zeros(N_omega,1);
%Trans_lambda_12 = zeros(N_omega,1);
%Trans_lambda_13 = zeros(N_omega,1);

%Q_omega_bulk = zeros(N_omega, length(ind_bulk));
%Q_omega_subvol = zeros(N_omega, N);
Q_w_AB = zeros(N_omega, length(ind_bulk));
Q_w_subvol = zeros(N_omega, N);

for omega_loop = 1:N_omega % Loop through all frequencies
    
    % Add space in Command Window
    disp(' ')
    
    % Record time at start of frequency loop
    t1 = toc;
    
    switch calc_approach
	    
	    case CalculationOption.direct
        
		%[ G_sys_2D, Trans, Q_omega_bulk(omega_loop, :), Q_omega_subvol(omega_loop, :) ] = direct_function(omega(omega_loop), r, epsilon(omega_loop)*ones(N,1), epsilon_ref, delta_V_vector, T_vector, ind_bulk, wave_type);
        	[ G_sys_2D, Trans, Q_w_AB(omega_loop, :), Q_w_subvol(omega_loop, :) ] = direct_function(omega(omega_loop), r, epsilon(omega_loop)*ones(N,1), epsilon_ref, delta_V_vector, T_vector, ind_bulk, wave_type);
        
	     case CalculationOption.iterative
        
			% not available for MATLAB version
        
    end % End direct vs. iterative approach
    

    % Calculate spectral transmission coefficient between object 1 and object 2 [dimensionless]
    Trans_omega_12(omega_loop) = trans_coeff_function_bulk( Trans, ind_bulk, 1, 2 );
    
    % Convert transmission coefficient for bulk objects into units of per
    % wavelength (rather than frequency) [1/(s*m)]
    %Trans_lambda_12(omega_loop) = Trans_omega_12(omega_loop)*(omega(omega_loop)^2)/(((2*pi)^2)*c_0);
    %Trans_lambda_13(omega_loop) = Trans_omega_13(omega_loop)*(omega(omega_loop)^2)/(((2*pi)^2)*c_0);
    
 %{       
    %%%%%%%%%%%%%%%%
    % Save results %
    %%%%%%%%%%%%%%%%
	
    % Save all workspace variables
    if output.save_workspace
        save([filePath_st.main, '/', file_name_saved])
    end
    
        save_DSGF_TRANS_matrix(output, filePath_st, omega_loop, G_sys_2D, Trans);
    
%}    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print status to Command window %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % End time for one frequency loop
    t6 = toc;

    % Print time for one frequency loop
    disp(['t = ' num2str(t6-t1) ' seconds = ' num2str((t6-t1)/60) ' minutes for one frequency loop'])

    % Print frequencies remaining to Command Window
    disp([num2str(length(omega) - omega_loop) ' frequencies remaining'])


end % End loop through all frequencies


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate total heat dissipation in each subvolume %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[Q_total_subvol, Q_density_subvol, Q_total_subvol_matrix, Q_total_bulk] = total_heat_dissipation_in_subvol(N, omega, Q_omega_subvol, delta_V_vector, r, ind_bulk);
[Q_t_subvol, Q_density_subvol, Q_total_subvol_matrix, Q_t] = total_heat_dissipation_in_subvol(N, omega, Q_w_subvol, delta_V_vector, r, ind_bulk);


%% 
%%%%%%%%%%%%%%%%%%
% Plot heat maps %
%%%%%%%%%%%%%%%%%%
    
    % Heat map in [W] and [W/m^3]
    if strcmp('sample',discretization_type)
        %subvol_heatmap_plotting(r, L_sub_vector,Q_total_subvol, Q_density_subvol, show_axes, output, filePath_st.main,ind_bulk,r_1,r_2,N);
	subvol_heatmap_plotting(r, L_sub_vector,Q_t_subvol, Q_density_subvol, show_axes, output, filePath_st.main,ind_bulk,r_1,r_2,N);
    
    elseif strcmp('user_defined',discretization_type)
        %subvol_heatmap_plotting_user_defined(r, L_sub_vector, Q_total_subvol, Q_density_subvol, show_axes, output, filePath_st.main, N, N1);
	subvol_heatmap_plotting_user_defined(r, L_sub_vector, Q_t_subvol, Q_density_subvol, show_axes, output, filePath_st.main, N, N1);
    
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate spectral and total conductance at temperature, T_conductance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If flag signals that conductance should be calculated between bulk
% objects.
if output.conductance

    for i = 1: length(T_conductance)
        
    %[ G_omega_bulk_12(:,i), G_bulk_12(:,i) ] = conductance_bulk( Trans_omega_12, T_conductance(i), omega );
    [ G_w_AB(:,i), G_t_AB(:,i) ] = conductance_bulk( Trans_omega_12, T_conductance(i), omega );
    
    end
    
    spectral_conductance_plot(omega, G_w_AB(:,:), output, filePath_st.main); % Updated on June 28, 2024
    %spectral_conductance_plot(omega, G_w_AB(:,3), output, filePath_st.main); % used in v1
    %spectral_conductance_plot(omega, G_omega_bulk_12(:,3), output, filePath_st.main); % used in v0
    
    %{
    [ G_omega_bulk_12, G_bulk_12 ] = conductance_bulk( Trans_omega_12, T_conductance, omega );

	spectral_conductance_plot(omega, G_omega_bulk_12, output, filePath_st.main);
    %}

end % End conductance calculations


%%%%%%%%%%%%%%%%
% Save results %
%%%%%%%%%%%%%%%%

clear G_sys_2D;

% Save all workspace variables
if output.save_workspace
    save([filePath_st.main, '/', file_name_saved])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print status to Command Window %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_memory_consumption(struct2cell(whos));

end % End function
