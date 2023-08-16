%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSGF: User Inputs
%
%
% DESCRIPTION: This script should be edited with all user inputs for DSGF 
%              simulations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace and close all figures
clear, clc, close all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************START OF USER INPUTS***************************%

%*******************************DESCRIPTION*******************************%

% Short description of system you are modeling (this will be used to name
% the saved files)

description = 'Sample_2_dipoles_r1_5nm_r2_5nm_d_10nm_N_2';        %  Sample_2_spheres_r1_50nm_r2_50nm_d_10nm_N_16  '2cubes_total' User-defined_2_cubes_L1_500nm_L2_500nm_d_500nm_N_72

%********************SELECTION OF TYPE OF SIMULATION *********************%

% Choose between sample (spheres) or user-defined (cubes)

discretization_type = 'sample'; 

if strcmp('sample',discretization_type) 

    %********************DISCRETIZATION OF EACH OBJECT********************%
    % 
    % Define the discretization for each bulk object. In the sample, each bulk object needs 
    % its own discretization. The discretizations can be taken from the
    % pre-made samples or defined by the user.  
    %
    % The number at the end of the chosen discretization represents the 
    % number of subvolumes in that discretization. 
    %
    % Pre-made sample discretization options:
    %     Discretization.sphere_*
    %     Discretization.cube_*
    %
    % Example with two sample discretizations chosen:
    %      discretization = {Discretization.sphere_8, Discretization.sphere_8};
    %

    discretization = {Discretization.cube_125, Discretization.cube_125};
    
    %**************************SCALE EACH OBJECT**************************%

    % Characteriztic length for scaling the discretized lattice of each bulk
    % object.
    %
    % If a pre-made sample is chosen, the characteristic length is:
    %     sphere: radius
    %     dipole: radius
    %     cube: side length
    %
    % If a user-defined input is chosen, the characteristic length is the
    % scaling factor of the user-input cubic lattice.
    %

    L_char = [500.e-9, 500.e-9]; % [m]
   
    %**********************DISTANCE BETWEEN OBJECTS***********************%
    
    % Distance between the objects
    d =500.e-9; %[m]

    %************************ORIGIN OF EACH OBJECT************************%

    % Matrix containing Cartesian coordinates of the origin of each object.
    
    %origin = [0,0,0; 110e-9, 0, 0]; % [m]
    %origin = [0,0,0; (L_char(1) + d + L_char(2)), 0, 0]; % [m]  for sphere
    origin = [0,0,0; (L_char(1)/2 + d + L_char(2)/2), 0, 0]; % [m] for cube     

    
elseif strcmp('user-defined',discretization_type)
    
    %********************DISCRETIZATION OF THE SYSTEM*********************%
    %
    % Define the discretization for the system (2 group of objects).
    % The user should modify the discretization and delta_V parameters
    % according to the name of the file with the desired user-defined discretization. 
    % These files are generated using matlab scripts.
        
    discretization = "2_thin_films_Lx500nm_Ly500nm_Lz500nm_d500nm_N72_discretization";
    delta_V = "2_thin_films_Lx500nm_Ly500nm_Lz500nm_d500nm_N72_delta_V_vector";

end

%********************************MATERIAL*********************************%
% Options:
%     'SiO_2'
%     'SiC'
%     'SiN'
%     'user_defined'

material = Material.SiO_2;


%********************DIELECTRIC FUNCTION OF BACKGROUND********************%

% The dielectric function of the background reference medium must be purely
% real-valued.

epsilon_ref = 1;


%************************FREQUENCY DISCRETIZATION*************************%

% Vector of angular frequencies at which simulations will be run.
% Vector is of dimension (N_omega x 1)
% Suggestions of wavelength range:
%           SiO2: uniform_omega(5e-6, 25e-6, 100);
%           SiC: uniform_omega(9.92e-6, 13.42e-6, 200);
%           SiN: uniform_omega(8e-6, 90e-6, 300);

% Wavelength [lambda] limits are provided
[omega] = uniform_omega(5e-6, 25e-6, 100); 

%Frequencies in [rad/s] limits are provided
%[omega] = non_uniform_omega(material);


%***********************TEMPERATURE OF EACH OBJECT************************%

T = [300, 400]; % [K]


%****************TEMPERATURE FOR CONDUCTANCE CALCULATIONS*****************%

% Temperature at which the spectral conductance will be calculated.

T_cond = [200,250,300,350,400]; % [K]


%************************** WAVE TYPE (OPTIONAL)**************************%

% Specify the contribution to be accounted: total, propagating only or
% evanescent only. Choose between "total" or "prop". 

wave_type = "total"; 


%*****************************DESIRED OUTPUTS*****************************%

% Output the power dissipated in every subvolume?
output.power_dissipated_subvol = true;

% Output the power dissipated in each bulk object?
output.power_dissipated_bulk = true;

% Output the heatmap into slices?
output.heatmap_sliced = false;

% Output the total and spectral conductance for each bulk object?
output.conductance = true;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = false;

% Output the DSGF matrices for every frequency?
output.DSGF_matrix = false;

% Output the heat transfer coefficient?
output.heat_transfer_coefficient = true;

% Save figures?
output.save_fig = true;

% figure format
output.figure_format = FigureFormat.fig;

% Save all Workspace variables in .mat file?
output.save_workspace = true;

                           
%***************************END OF USER INPUTS****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp('sample',discretization_type) 
    delta_V = '';
elseif strcmp('user-defined',discretization_type)
    L_char = '';
    origin = '';
end    


DSGF_main(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, discretization_type, L_char, origin, delta_V);

