function [] = DSGF_main(description, discretization, L_char, origin, material, T, T_conductance, epsilon_ref, omega, observation_point, output)

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
show_axes = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set calculation options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculation approach
% 'direct' = direct matrix inversion using 'mldivide' operator
% 'iterative' =  iterative solver from Martin et al.
calc_approach = CalculationOption.direct;

% Conduct convergence analysis?
convergence_analysis = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************DISCRETIZATION*******************************%

% Number of bulk objects
N_bulk = length(discretization);

% Determine file structure of bulk object discretization and extract
% discretization information.
N_each_object = zeros(N_bulk,1);       % Preallocate
volume = zeros(N_bulk,1);              % Preallocate
r_each_object = cell(N_bulk,1);        % Preallocate
ind_bulk = 2.^[0:N_bulk];            % Preallocate
delta_V_each_object = cell(N_bulk,1);  % Preallocate
L_sub_each_object = cell(N_bulk,1);    % Preallocate
for ii = 1:N_bulk % Loop through all bulk objects
    if isenum(discretization{ii}) % Sample discretization is specified

	[r_each_object{ii}, N_each_object(ii), delta_V_each_object{ii}, L_sub_each_object{ii}, volume(ii)] = read_sample_discretization(discretization{ii}, L_char(ii));

    else % User-defined discretization is specified
        discFile = discretization{ii}; % File name of discretization
        discDir = "Input_parameters/Discretizations/User_defined"; % Directory where discretization is stored

        % Import unscaled discretization of each object
        r_each_object{ii} = readmatrix(append(append(append(discDir, '/'), discFile), '.txt'));

        % Scale discretization
        r_each_object{ii} = L_char(ii).*r_each_object{ii};

        % Number of subvolumes in each bulk object
        [N_each_object(ii),~] = size(r_each_object{ii});

        % Subvolume size for each object (uniform discretization)
        [L_sub, delta_V] = calculate_Lsub_uniform(r_each_object{ii});
        delta_V_each_object{ii} = ones(N_each_object(ii), 1).*delta_V; % Volume of subvolumes in each object (uniform discretization)
        L_sub_each_object{ii} = ones(N_each_object(ii), 1).*L_sub;     % Length of side of a cubic subvolume in each object (uniform discretization)

    end % End if sample or user-defined discretization

    % Move the center-of-mass of each object to the origin [0,0,0]
    r_each_object{ii} = center_of_mass(r_each_object{ii});

    % Move each discretization to its user-specified origin
    r_each_object{ii} = r_each_object{ii} + repmat(origin(ii,:), N_each_object(ii), 1);

end % End loop through bulk objects

% Discretized lattice including subvolumes of all objects in one matrix (N x 3 matrix)
r = cell2mat(r_each_object);

% Total number of subvolumes
[N,~] = size(r);

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



%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%

discretization_plotting(r, L_sub_vector, N, show_axes, output, filePath_st.main);


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
	case Material.SiO_2

	    epsilon = SiO2_dielectric_function(omega, constants); % (N x 1) vector of all dielectric functions for every frequency

	case Material.SiC

	    epsilon = SiC_dielectric_function(omega, constants); % (N x 1) vector of all dielectric functions for every frequency


	case Material.SiN

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
G_omega = zeros(N_omega,1);
G_omega_eV = zeros(N_omega,1);
Trans_omega_12 = zeros(N_omega,1);
Trans_omega_13 = zeros(N_omega,1);
Trans_lambda_12 = zeros(N_omega,1);
Trans_lambda_13 = zeros(N_omega,1);
Q_omega_bulk = zeros(N_omega, length(ind_bulk));
Q_omega_subvol = zeros(N_omega, N);

for omega_loop = 1:N_omega % Loop through all frequencies
    
    % Add space in Command Window
    disp(' ')
    
    % Record time at start of frequency loop
    t1 = toc;
    
    if calc_approach == CalculationOption.direct
        
        [ G_sys_2D, Trans, Q_omega_bulk(omega_loop, :), Q_omega_subvol(omega_loop, :) ] = direct_function(omega(omega_loop), r, epsilon(omega_loop)*ones(N,1), epsilon_ref, delta_V_vector, T_vector, ind_bulk);
        
    elseif calc_approach == CalculationOption.iterative
        
        % THIS IS STILL UNDER CONSTUCTION
        
    end % End direct vs. iterative approach
    

    % Calculate spectral transmission coefficient between object 1 and object 2 [dimensionless]
    Trans_omega_12(omega_loop) = trans_coeff_function_bulk( Trans, ind_bulk, 1, 2 );
    
    % Convert transmission coefficient for bulk objects into units of per
    % wavelength (rather than frequency) [1/(s*m)]
    %Trans_lambda_12(omega_loop) = Trans_omega_12(omega_loop)*(omega(omega_loop)^2)/(((2*pi)^2)*c_0);
    %Trans_lambda_13(omega_loop) = Trans_omega_13(omega_loop)*(omega(omega_loop)^2)/(((2*pi)^2)*c_0);
    
        
    %%%%%%%%%%%%%%%%
    % Save results %
    %%%%%%%%%%%%%%%%
    
    % Save all workspace variables
    if output.save_workspace
        save([filePath_st.main, '/', file_name_saved])
    end

    % Export DSGF matrix for this frequency loop
    if output.DSGF_matrix
        % File name where results will be saved (based on what frequency band is chosen)
        fileName_DSGF = ['DSGFmatrix_omega' num2str(omega_loop)];
        t2 = toc;
        % Export DSGF matrix
        writematrix(G_sys_2D, [filePath_st.DSGF, '/', fileName_DSGF, '.csv']);
        t3 = toc;
        disp(['Time to save DSGF matrix as a .csv file  = ' num2str(t3-t2) ' s = ' num2str((t2-t1)/60) ' minutes'])
    end % End save DSGF matrix for this frequency loop


    % Export transmission coefficient matrix for this frequency loop
    if output.transmission_coefficient_matrix
        % File name where results will be saved (based on what frequency band is chosen)
        fileName_Trans = ['TransMatrix_omega' num2str(omega_loop)];

        t4 = toc;
        % Export spectral transmission coefficient matrix
        writematrix(Trans, [filePath_st.Trans, '/', fileName_Trans, '.csv']);
        %save(filePath_st.Trans, 'omega', 'Trans_12_omega_DSGF')
        t5 = toc;
        disp(['Time to save transmission coefficient matrix as a .csv file  = ' num2str(t5-t4) ' s = ' num2str((t4-t3)/60) ' minutes'])
    end



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


[Q_total_subvol, Q_density_subvol, Q_total_subvol_matrix] = total_heat_dissipation_in_subvol(N, omega, Q_omega_subvol, delta_V_vector, r);


%% 
%%%%%%%%%%%%%%%%%%
% Plot heat maps %
%%%%%%%%%%%%%%%%%%

subvol_heatmap_plotting(r, L_sub_vector, Q_total_subvol, show_axes, output, filePath_st.main);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate spectral and total conductance at temperature, T_conductance %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If flag signals that conductance should be calculated between bulk
% objects.
if output.conductance

    [ G_omega_bulk_12, G_bulk_12 ] = conductance_bulk( Trans_omega_12, T_conductance, omega );

	spectral_conductance_plot(omega, G_omega_bulk_12, output, filePath_st.main);

end % End conductance calculations


%%%%%%%%%%%%%%%%
% Save results %
%%%%%%%%%%%%%%%%

% Save all workspace variables
if output.save_workspace
    save([filePath_st.main, '/', file_name_saved])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print status to Command Window %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print out workspace memory required
vars = struct2cell(whos);
workspaceMem = sum(cell2mat(vars(3,:)))/1e9;
disp('----------------------------------------------------------------')
disp(['Workspace memory required is ', num2str(workspaceMem), ' GB'])
disp('')
disp('----------------------------------------------------------------')
%memory


end % End function
