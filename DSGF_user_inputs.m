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

description = '2_cubes_SiO2_Lchar_500nm_d_100nm_N_2000';% 

%********************SELECTION OF TYPE OF SIMULATION *********************%

% Choose between sample or user_defined 

spatial_discretization_type = 'sample'; 

if strcmp('sample',spatial_discretization_type) 
    
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

    discretization = {Discretization.cube_1000, Discretization.cube_1000};
    
    %**************************SCALE EACH OBJECT**************************%

    % Characteristic length for scaling the discretized lattice of each bulk
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

    L_char = [500.e-9, 500.e-9]; % [m] %[50.e-9, 50.e-9]
   
    %**********************DISTANCE BETWEEN OBJECTS***********************%
    
    % Distance between the objects
    d =100.e-9; %[m]


    
elseif strcmp('user_defined',spatial_discretization_type)
    
    %********************DISCRETIZATION OF THE SYSTEM*********************%
    %
    % Define the discretization for the system (2 group of objects).
    % The user should modify the discretization and delta_V parameters
    % according to the name of the file with the desired user-defined discretization. 
    % These files are generated using matlab scripts.
        
    discretization = "2_membranes_Lx1000nm_Ly1000nm_Lz100nm_d100nm_N1600_discretization"; %"2_cubes_Lx500nm_Ly500nm_Lz500nm_d500nm_N576_discretization"; %
    delta_V = "2_membranes_Lx1000nm_Ly1000nm_Lz100nm_d100nm_N1600_delta_V_vector"; %"2_cubes_Lx500nm_Ly500nm_Lz500nm_d500nm_N576_delta_V_vector"; %

end

%********************************MATERIAL*********************************%
% Options:
%     'SiO2'
%     'SiC'
%     'Si3N4'
%     'user_defined'

material = Material.SiO2;


%********************DIELECTRIC FUNCTION OF BACKGROUND********************%

% The dielectric function of the background reference medium must be purely
% real-valued.

epsilon_ref = 1;


%************************FREQUENCY DISCRETIZATION*************************%

% Vector of angular frequencies at which simulations will be run.
% Vector is of dimension (N_omega x 1)

% Insert the name of the file with the spectral discretization.
spectral_discretization = "SiO2_100_uniform_5_25_um.csv"; % "SiC_100_nonuniform_140_190_Trad_s.csv"

%***********************TEMPERATURE OF EACH OBJECT************************%

T = [300, 800]; % [K]


%****************TEMPERATURE FOR CONDUCTANCE CALCULATIONS*****************%

% Temperature at which the spectral conductance will be calculated.

T_cond = [200,250,300,350,400]; % [K]


%*****************************DESIRED OUTPUTS*****************************%

% Save all Workspace variables in .mat file?
output.save_workspace = true;

% Save figures?
output.save_fig = true;

% figure format
output.figure_format = FigureFormat.fig;

% Output the total and spectral conductance for each bulk object?
output.conductance = true;

% Output the power dissipated in every subvolume?
output.power_dissipated_subvol = true;

% Output the power dissipated in each bulk object?
output.power_dissipated_bulk = true;

% Output the heatmap into slices?
output.heatmap_sliced = false;

% Output the DSGF matrices for every frequency?
output.DSGF_matrix = false;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = false;

                           
%***************************END OF USER INPUTS****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[omega] = readmatrix(append("Library/spectral_discretizations/", spectral_discretization));

if strcmp('sample',spatial_discretization_type) 
    delta_V = '';
    wave_type = "total";
elseif strcmp('user_defined',spatial_discretization_type)
    L_char = '';
    d = '';
    wave_type = "total";
end    


%DSGF_main(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, discretization_type, L_char, origin, delta_V);
DSGF_main(description, discretization, material, T, T_cond, epsilon_ref, omega, wave_type, output, spatial_discretization_type, L_char, delta_V, d);

