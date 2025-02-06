MATLAB DSGF code

Last Update: Feb 6, 2025

The DSGF (discrete system Green's function) is a numerical framework to simulate radiative heat transfer between micro/nanostructured objects. The DSGF method is based on fluctuational electrodynamics and is applicable to near-field and far-field problems. The DSGF method and solver were developed at the University of Utah. 


1. SYSTEM REQUIREMENT

   A MATLAB license from MathWorks is required. The DSGF code has been tested on MATLAB version R2021a. 

2. INSTALLATION

   No installation is required. Download the DSGF_matlab source code.

3. DEMO AND INSTRUCTIONS FOR USE
   
   The only file to be modified is: DSGF_user_inputs_emissivity.

   Step 1: Write a description of your simulation. The description will be available in the results table.

   Step 2: Define discretization for your simulation. Different parameters need to be modified depending on your selection in Step 2.
        For the ‘sample’ provided, modify: discretization, L_char, and d. 
    
   Step 3: Select material. 

   Step 4: Define dielectric function of the background reference medium.

   Step 5: Define the frequency discretization. The option uniform_lambda is defined in terms of wavelength [m] while the options uniform_omega and non_uniform_omega are defined in angular frequency [rad/s]. 

   Step 6: Define the only one temperature required for the power emitted calculation.

   Step 7: Select outputs of the simulation.

   Step 8: Select the cross-sectional area to normalize the absorption cross-section. 

   Step 9: In the MATLAB editor, press the Run button.
   
   Key Assumptions: This code solves for the spectral, hemispherical emissivity for objects uniformly discretized, having the same size, shape, equal edge-to-edge spacing, and normalized to the same cross-sectional area. 
		    If the user wants to run emissivity calculations on objects of different sizes, shapes, and cross-sectional areas, the user must modify the DSGF_Thermal_Emission code to account for the differences. 


   The results are stored in a folder named with the time the simulation was launched.
