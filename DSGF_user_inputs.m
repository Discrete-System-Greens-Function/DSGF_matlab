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

description = '2spheres';


%**********************DISCRETIZATION OF EACH OBJECT**********************%
% 
% Define the discretization for each bulk object. Each bulk object needs 
% its own discretization. The discretizations can be taken from the
% pre-made samples or defined by the user.  
%
% If a pre-made sample is chosen, the user can choose the discretization 
% from a list of pre-defined files. Files are included for spheres, cubes,
% and thin films with variable number of subvolumes and for diples defined
% by a single subvolume.  The number at the end of the chosen 
% discretization represents the number of subvolumes in that
% discretization. 
%
% Pre-made sample discretization options:
%     Discretization.sphere_*
%     Discretization.cube_*
%     Discretization.thin_film_*
%     Discretization.dipole
%
% Example with two sample discretizations chosen:
%      discretization = {Discretization.sphere_8, Discretization.sphere_8};
%
% If a user-defined input is chosen, the user-defined .txt file of the 
% discretization should be stored in the directory 
% Input_parameters/Discretizations/User_defined.  The inputs are then
% strings of the file name of the user-defined discretization.
%
% Example of two user-defined discretizations:
%      discretization = {"user_defined_discretization_1", "user_defined_discretization_2"}
%
% Sample discretizations and user-defined discretizations can be
% mixed-and-matched.
%
% Example of one sample discretization and one user-defined discretization
%      discretization = {Discretization.sphere_8 "user_defined_discretization_2"}
%

discretization = {Discretization.sphere_8, Discretization.sphere_8};


%**************************VOLUME OF EACH OBJECT**************************%

% Characteriztic length (e.g., radius, side length, film thickness)
L_char = [50e-9, 50e-9]; % [m]



%**************************ORIGIN OF EACH OBJECT**************************%

% Matrix containing Cartesian coordinates of the origin of each object.

origin = [0,0,0;
          110e-9, 0, 0]; % [m]


%********************************MATERIAL*********************************%
% Options:
%     'SiO_2'
%     'SiC'
%     'SiN'
%     'user_defined'

material = Material.SiO_2;


%***********************TEMPERATURE OF EACH OBJECT************************%

T = [300, 800]; % [K]


%****************TEMPERATURE FOR CONDUCTANCE CALCULATIONS*****************%

% Temperature at which the spectral conductance will be calculated.

T_cond = 300; % [K]


%********************DIELECTRIC FUNCTION OF BACKGROUND********************%

% The dielectric function of the background reference medium must be purely
% real-valued.

epsilon_ref = 1;


%************************FREQUENCY DISCRETIZATION*************************%

% Vector of angular frequencies at which simulations will be run.
% Vector is of dimension (N_omega x 1)

N_omega = 100;
lambda = linspace(5e-6, 25e-6, N_omega); % [m]
c_0 = 299792458;            % Speed of light in vacuum [m/s]
omega = (2*pi*c_0./lambda).'; % [rad/s]


%**********************OBSERVATION POINT (OPTIONAL)***********************%

% Cartesian coordinates of observation point at which the LDOS will be
% calculated

observation_point = [0,0,60e-9]; % [m]


%*****************************DESIRED OUTPUTS*****************************%

% Output the power dissipated in every subvolume?
output.power_dissipated_subvol = true;

% Output the power dissipated in each bulk object?
output.power_dissipated_bulk = true;

% Output the total and spectral conductance for each bulk object?
output.conductance = true;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = true;

% Output the DSGF matrices for every frequency?
output.DSGF_matrix = true;

% Output the local density of states (LDOS)?
% NOTE: This only works when an observation point is given as an input.
output.LDOS = true;

% Output the heat transfer coefficient?
output.heat_transfer_coefficient = true;

% Save figures?
output.save_fig = true;

% Save all Workspace variables in .mat file?
output.save_workspace = true;

                           
%***************************END OF USER INPUTS****************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DSGF_main(description, discretization, L_char, origin, material, T, T_cond, epsilon_ref, omega, observation_point, output);
