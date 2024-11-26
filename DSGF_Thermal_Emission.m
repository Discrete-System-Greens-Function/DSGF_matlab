
function [] = DSGF_Thermal_Emission(description, discretization, material, T_equilibrium, epsilon_ref, omega, output, L_char, d);

%Thermal emission code using the discrete system Green's function method
%Created by Joseph C. McKay
%See the manuscript "Thermal emission from agglomerated polaritonic dielectric particles" for the derivation 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates:
% System Green's functions
% Spectral, hemispherical emissivity
% The total power emitted
%
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
constants.epsilon_0 = 8.85418782e-12; % Permittivity of free space [F/m]
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
    
    
    [N_each_object, volume, r_each_object, ind_bulk, delta_V_each_object, L_sub_each_object,origin,A_c] = read_discretization_emissivity(discretization, L_char,d);
    
    % Discretized lattice including subvolumes of all objects in one matrix (N x 3 matrix)
    r = cell2mat(r_each_object);
    
    % Total number of subvolumes
    [N,~] = size(r);
    
    r_1 = r(1:ind_bulk(2)-1,:);
    r_2 = r(ind_bulk(2):N,:);
    
    % Subvolume size for all N subvolumes (N x 1 vectors)
    delta_V_vector = cell2mat(delta_V_each_object); % Volume of subvolumes for all N subvolumes
    L_sub_vector = cell2mat(L_sub_each_object);     % Length of side of a cubic subvolume for all N subvolumes
    
    %These deltaV and A_c are based on uniformly discretized particles of
    %the same size and shape. If the user modifies the discretization or
    %shape, the user must update these two variables (delta_V_subvolume and A_c)
    delta_V_subvolume = L_sub_vector(2)^3;% volume of a single subvolume to calculate emissivity
    A_c = A_c(1);                         % cross-sectional area to calculate emissivity
    N_objects = length(L_char);           % number of objects in the system


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate DSGF and the spectral, hemispherical emissivity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k_0 = omega*sqrt(constants.epsilon_0*constants.mu_0);% Vacuum wave vector [1/m]

%Preallocate vectors
emissivity = zeros(N_omega,1);


for omega_loop = 1:N_omega % Loop through all frequencies
    
    % Add space in Command Window
    disp(' ')
    
    % Record time at start of frequency loop
    t1 = toc;
    
    switch calc_approach
	    
	    case CalculationOption.direct
       
            %Calculate the system Green's function
            [G_sys_2D] = G_sys_function(omega(omega_loop), r, epsilon(omega_loop)*ones(1,N), epsilon_ref, delta_V_vector');
            
            %Calculate the emissivity
            [emissivity(omega_loop)] = Q_abs_eff(G_sys_2D,k_0(omega_loop), epsilon(omega_loop),epsilon_ref, delta_V_subvolume,A_c);

	     case CalculationOption.iterative
        
			% not available for MATLAB version
        
    end % End direct vs. iterative approach
    

   
     
    %%%%%%%%%%%%%%%%
    % Save results %
    %%%%%%%%%%%%%%%%
%{   	
    % Save all workspace variables
    if output.save_workspace
        save([filePath_st.main, '/', file_name_saved])
    end
%}    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print status to Command window %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % End time for one frequency loop
    t6 = toc;

    % Print time for one frequency loop
    %disp(['t = ' num2str(t6-t1) ' seconds = ' num2str((t6-t1)/60) ' minutes for one frequency loop'])

    % Print frequencies remaining to Command Window
    %disp([num2str(length(omega) - omega_loop) ' frequencies remaining'])

end % End loop through all frequencies













%%%%%%%%%%%%%%%%%%%%%%%
% Plot discretization %
%%%%%%%%%%%%%%%%%%%%%%%


    discretization_plotting(r, L_sub_vector, N, show_axes, output, filePath_st.main,ind_bulk,r_1,r_2);



% % String name describing simulation geometry
% simulation_geometry = [description '_R' num2str((1e9)*radius_1) 'nm_R' num2str((1e9)*radius_2) 'nm_dc' num2str((1e9)*d_center) 'nm_N' num2str(N1 + N2)];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot dielectric function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dielectric_function_plotting(omega, epsilon, material, N_omega, output, filePath_st.main);
   

% % Save discretization figure file
% if save_fig == 1
%     fig_path_2 = [saveDir '/' file_name_saved '_discretization.fig'];
%     saveas(Fig_discretization, fig_path_2)
%     clear Fig_dielectric_function Fig_discretization % Remove previous plot handles
% end
emissivity = -emissivity;%sign convention in all derivations were - for energy leaving the objects   

  

%------Calculating Emissivity per particle--------------
emissivity_per_particle = emissivity./N_objects;
 

%-------Calculate the power emitted per particle--------
% Mean energy of an electromagnetic state for all subvolumes
theta = (constants.h_bar*omega)./(exp(constants.h_bar*omega./(constants.k_b*T_equilibrium)) - 1);
 

%Calculate black body intensity
Ib = ((k_0.^2).*theta./(4*pi^3));
 
 
%Total power emitted per particle [W]
P_em_spectral =  emissivity_per_particle.*Ib.*A_c;
P_em_total = (4*pi).*trapz(omega,P_em_spectral);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the spectral, hemispherical emissivity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------Display emissivity per particle--------------
figure;
plot(omega, emissivity_per_particle, 'k');hold on;
xlabel('Radial frequency, \omega [rad/s]');
ylabel('Spectral emissivity');
title(['Spectral Emissivity of ' string(material) ' per Object'], 'fontsize', 14);





%%%%%%%%%%%%%%%%
% Save results %
%%%%%%%%%%%%%%%%

clear G_sys_2D;%to save memory, remove the large system Green's functions


% Save all workspace variables
if output.save_workspace
    save([filePath_st.main, '/', file_name_saved])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print status to Command Window %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display_memory_consumption(struct2cell(whos));

end
