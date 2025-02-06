function [r_each_object, N_each_object, delta_V_each_object, L_sub_each_object] = read_user_discretization(discretization, L_char)
% reads the user discretization file
% this file must be a txt file
%
%	Inputs:
%		discretization - filename of the custom discretization
%		L_char - characteristic length
%	
%	Outputs:
%		r_each_object - radius of bulk object
%		N_each_object - number of subvolumes in object
%		delta_V_each_object - volume of subvolumes in object
%		L_sub_each_object - length of side of a cubic subvolume in object


        discFile = discretization; % File name of discretization
        discDir = "Library/Discretizations/User_defined"; % Directory where discretization is stored

        % Import unscaled discretization of each object
        r_each_object = readmatrix(append(append(append(discDir, '/'), discFile), '.txt'));

        % Scale discretization
        r_each_object = L_char.*r_each_object;

        % Number of subvolumes in each bulk object
        [N_each_object,~] = size(r_each_object);

        % Subvolume size for each object (uniform discretization)
        [L_sub, delta_V] = calculate_Lsub_uniform(r_each_object);
        delta_V_each_object = ones(N_each_object, 1).*delta_V; % Volume of subvolumes in each object (uniform discretization)
        L_sub_each_object = ones(N_each_object, 1).*L_sub;     % Length of side of a cubic subvolume in each object (uniform discretization)
end
