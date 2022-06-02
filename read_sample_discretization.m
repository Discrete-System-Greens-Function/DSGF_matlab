function [r_each_object, N_each_object, delta_V_each_object, L_sub_each_object, volume] = read_sample_discretization(discretization, L_char)
% reads the sample discretizations
%
%	Inputs:
% 		discretization - enum of the discretization
%		L_char - characteristic length of bulk object
%
%	Outputs:
%		r_each_object - radius of the bulk object
%		N_each_object - number of subvolumes in bulk object
%		delta_v_each_object - subvolume size in bulk object
%		L_sub_each_object - Length of side of a cubic subvolume in bulk object
%		volume - volume of bulk object


	discFile = string(discretization);       % File name of discretization
	geometry = extractBefore(discFile, '_');     % Geometry of bulk object
	discDir = append('Input_parameters/Discretizations/', geometry); % Directory where discretization is stored

	% Import unscaled discretization of each object
	r_each_object = xlsread(append(append(append(discDir, '/'), discFile), '.xlsx'));

	% Number of subvolumes in each bulk object
	[N_each_object,~] = size(r_each_object);

	% Scale sample discretizations based on the geometry
	switch (geometry)
	    case "sphere"
		volume = (4/3)*pi*(L_char.^3);  % Volume of sphere [m^3]
	    case "cube"
		volume = (L_char.^3);           % Volume of cube [m^3]
	    case "thin_film"

	    case "dipole"
		volume = (4/3)*pi*(L_char.^3);  % Volume of spherical dipole [m^3]
	end % End switch-case through geometries

	% Subvolume size for each object (uniform discretization)
	delta_V_each_object = ones(N_each_object, 1).*(volume/N_each_object); % Volume of subvolumes in each object (uniform discretization)
	L_sub_each_object = delta_V_each_object.^(1/3);  % Length of side of a cubic subvolume in each object (uniform discretization)

	% Scale discretization
	r_each_object = L_sub_each_object(1).*r_each_object;

end
