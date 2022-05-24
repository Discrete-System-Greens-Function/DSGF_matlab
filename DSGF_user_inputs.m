%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lindsay Walter, Livia Correa, Joseph McKay, Jan Cas
% Discrete System Green's Function user inputs
% Updated 5/23/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code implements the Discrete System Green's Function (DSGF) approach
% to model near-field radiative heat transfer between two objects.


% Clear Workspace and close all figures
clear, clc, close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%**************************START OF USER INPUTS**************************%

%******************DESCRIPTION******************%

% Short description of system you are modeling (this will be used to name
% the saved files)

description = '2spheres';


%******************DISCRETIZATION OF EACH OBJECT******************%
% Options:
%     'sphere_*'
%     'cube_*'
%     'thin_film_*'
%     'dipole'
%     'user_defined'

discretization = ["sphere_8", "sphere_8"];


%******************VOLUME OF EACH OBJECT******************%

% Characteriztic length (e.g., radius, side length, film thickness)
L_char = [50e-9, 50e-9]; % [m]

% Volume 
volume = (4/3)*pi*(L_char.^3);  % [m^3]



%******************ORIGIN OF EACH OBJECT******************%

% Matrix containing Cartesian coordinates of the origin of each object.

origin = [0,0,0;
          110e-9, 0, 0]; % [m]


%******************MATERIAL******************%
% Options:
%     'SiO_2'
%     'SiC'
%     'SiN'
%     'user_defined'

material = 'SiO_2';


%******************TEMPERATURE OF EACH OBJECT******************%

T = [300, 800]; % [K]


%******************TEMPERATURE FOR CONDUCTANCE CALCULATIONS******************%

% Temperature at which the spectral conductance will be calculated.

T_cond = 300; % [K]


%******************DIELECTRIC FUNCTION OF BACKGROUND******************%

% The dielectric function of the background reference medium must be purely
% real-valued.

epsilon_ref = 1;


%******************FREQUENCY DISCRETIZATION******************%

% Vector of angular frequencies at which simulations will be run.
% Vector is of dimension (N_omega x 1)

N_omega = 100;
lambda = linspace(5e-6, 25e-6, N_omega); % [m]
c_0 = 299792458;            % Speed of light in vacuum [m/s]
omega = (2*pi*c_0./lambda).'; % [rad/s]


%******************OBSERVATION POINT (OPTIONAL)******************%

% Cartesian coordinates of observation point at which the LDOS will be
% calculated

observation_point = [0,0,60e-9]; % [m]


%******************DESIRED OUTPUTS******************%

% 0: Do NOT output
% 1: Do output

% Output the power dissipated in every subvolume?
output.power_dissipated_subvol = 1;

% Output the power dissipated in each bulk object?
output.power_dissipated_bulk = 1;

% Output the total and spectral conductance for each bulk object?
output.conductance = 1;

% Output the transmission coefficient matrix?
output.transmission_coefficient_matrix = 1;

% Output the DSGF matrices for every frequency?
output.DSGF_matrix = 1;

% Output the local density of states (LDOS)?
% NOTE: This only works when an observation point is given as an input.
output.LDOS = 1;

% Output the heat transfer coefficient?
output.heat_transfer_coefficient = 1;

% Save figures?
output.save_fig = 1;

% Save all Workspace variables in .mat file?
output.save_workspace = 1;

                           
%**************************END OF USER INPUTS**************************%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DSGF_main(description, discretization, volume, origin, material, T, T_cond, epsilon_ref, omega, observation_point, output);











