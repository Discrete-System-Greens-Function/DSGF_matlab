function [emissivity] = Q_abs_eff(G_sys_2D,k_0, epsilon,epsilon_ref, delta_V,A_c)
% This function calculates spectral hemispherical emissivity for all thermal objects
% at a given frequency omega using the discrete system Green's function (DSGF) approach.

% INPUTS:  G_sys_2D      system Green's function matrix for all subvolume interactions (3N x 3N)
%          k_0           the vacuum wave vector for a single frequency
%          epsilon       dielectric function at a single frequency 
%          epsilon_ref   the background reference medium
%          delta_V       a single subvolume sizes [m^3]
%          A_c           the geometrical cross-section of each thermal object 
%
% OUTPUTS: emissivity    The spectral hemispherical emissivity for
%                        all thermal objects
 

%Create a diagonal matrix using the extracted diagonal elements
diagonalMatrix = diag(diag(G_sys_2D));

%spectral hemispherical emissivity  
emissivity = 2*pi*imag(sum(diag(imag(epsilon)*delta_V*((diagonalMatrix')+(epsilon-epsilon_ref)*(k_0^2)*delta_V*G_sys_2D*(G_sys_2D')))))/A_c;





