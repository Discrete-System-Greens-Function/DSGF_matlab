%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSGF: User Inputs for emissivity 
%
% Created by: Joseph C. McKay
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

description = '2_cubes_SiO2_Lchar_500nm_d_100nm_N_2000';% Modify this name for your problem 

%********************SELECTION OF TYPE OF SIMULATION *********************%

% Here is a sample 

spatial_discretization_type = 'sample'; 


    
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
    %     Discretization.sphere_*# of subvolumes
    %     Discretization.cube_*# of subvolumes
    %
    % Example with two sample discretizations chosen:
    %      discretization = {Discretization.sphere_8, Discretization.sphere_8};
    %

    discretization = {Discretization.sphere_8, Discretization.sphere_8};
    
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

    L_char = [50.e-9, 50.e-9]; % [m] %[50.e-9, 50.e-9]
    N_objects = length(L_char);% number of objects in the system
    %**********************DISTANCE BETWEEN OBJECTS***********************%
    
    % Distance between the objects
    d =100.e-9; %[m] edge-to-edge gap distance 

   

   


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



%****************TEMPERATURE FOR EMISSIVITY CALCULATIONS*****************%

% Define the equilibrium temperature
T_equilibrium = 300; %[K]

%*****************************DESIRED OUTPUTS*****************************%
%The thermal conductance has different boundary conditions than the thermal emissivity calculations 

% Save all Workspace variables in .mat file?
output.save_workspace = true;

% Save figures?
output.save_fig = true;

% figure format
output.figure_format = FigureFormat.fig;


% Output the heatmap into slices?
output.heatmap_sliced = false;

% Output the DSGF matrices for every frequency?
output.DSGF_matrix = false;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = false;

                           
%***************************END OF USER INPUTS****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[omega] = readmatrix(append("Library/spectral_discretizations/", spectral_discretization));




%------For thermal emissivity and power emitted calculations only----------'
    %These deltaV and A_c are based on uniformly discretized particles of
    %the same size and shape. If the user modifies the discretization or
    %shape, the user must update these two variables (delta_V_subvolume and A_c)
    
    %Define the cross-sectional area of the objects to normalize the absorption cross-section
    %The objects are assumed to have the same geometrical
    %A_c = (L_char.^2);            %cube - geometrical cross-section of the bulk object
    %A_c = pi*(L_char.^2);         %dipole - geometrical cross-section of the bulk object
    A_c = pi*(L_char(1).^2);       %sphere - geometrical cross-section of the bulk object
    
   
    


DSGF_Thermal_Emission(description, discretization, material, T_equilibrium, epsilon_ref, omega, output, L_char, d,A_c,N_objects);


