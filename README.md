MATLAB DSGF code
Last Update: February 20, 2024

1. SYSTEM REQUIREMENT
A MATLAB license from MathWorks is required. The DSGF code has been tested on MATLAB version R2021a. 

2. INSTALLATION
No installation is required. Download the DSGF_matlab source code.

3. DEMO AND INSTRUCTIONS FOR USE

The only file to be modified is: DSGF_user_inputs.

Step 1: Write a description of your simulation. The description will be available in the results table.

Step 2: Select the discretization_type between a ‘sample’ or a ‘user-defined’ simulation. ‘Sample’ is used for simulations with spheres, dipoles, and cubes. ‘User-defined’ is used for simulations with membranes.

Step 3: Define discretization for your simulation. Different parameters need to be modified depending on your selection in Step 2.
        For ‘sample’, modify: discretization, L_char, and d. 
        For ‘user-defined’, modify: discretization and delta_V with the names of the discretization files. These files are located in Library/Discretizations/User_defined. 

Step 4: Select material. 

Step 5: Define dielectric function of the background reference medium.

Step 6: Define the frequency discretization. The option uniform_lambda is defined in terms of wavelength [m] while the options uniform_omega and non_uniform_omega are defined in angular frequency [rad/s]. 

Step 7: Define temperatures for each object. These temperatures are used for power dissipation calculations.

Step 8: Define temperatures for conductance calculations.

Step 9: Select outputs of the simulation.

Step 10: In the MATLAB editor, press the Run button.

The results are stored in a folder named with the time the simulation was launched.
