function [N_each_object, volume, r_each_object, ind_bulk, delta_V_each_object, L_sub_each_object,origin,A_c] = read_discretization_emissivity(discretization, L_char, d)
% reads all the discretizations for all the bulk objects
%
%	Inputs:
%		discretization - vector of the all the bulk objects and their discretizations
%		L_char - vector of the characteristic lengths of the objects
%		origin - vector of the origins of all the bulk objects
%
%	Outputs:
%		N_each_object
%		volume
%		r_each_object
%		ind_bulk
%		delta_V_each_object
%		L_sub_each_object


	% Number of bulk objects
	N_bulk = length(discretization);
    
    
	discFile = string(discretization);       % File name of discretization
	geometry = extractBefore(discFile, '_');     % Geometry of bulk object
    
    % Matrix containing Cartesian coordinates of the origin of each object.
    % Origin sample discretizations based on the geometry
    switch (geometry)
        case "sphere"
            origin = [0,0,0; (L_char(1) + d + L_char(2)), 0, 0]; % [m]  for sphere
            A_c = pi*(L_char.^2);            %geometrical cross-section of the bulk object
            
        case "cube"
            origin = [0,0,0; (L_char(1)/2 + d + L_char(2)/2), 0, 0]; % [m] for cube
            A_c = (L_char.^2);            %geometrical cross-section of the bulk object
        case "thin_film"
            
        case "dipole"
            origin = [0,0,0; (L_char(1) + d + L_char(2)), 0, 0]; % [m]  for sphere
            A_c = pi*(L_char.^2);            %geometrical cross-section of the bulk object
    end % End switch-case through geometries

	% Determine file structure of bulk object discretization and extract
	% discretization information.
	N_each_object = zeros(N_bulk,1);       % Preallocate
	volume = zeros(N_bulk,1);              % Preallocate
	r_each_object = cell(N_bulk,1);        % Preallocate
	ind_bulk = ones(N_bulk,1);            % Preallocate
	delta_V_each_object = cell(N_bulk,1);  % Preallocate
	L_sub_each_object = cell(N_bulk,1);    % Preallocate
	for ii = 1:N_bulk % Loop through all bulk objects
	    if isenum(discretization{ii}) % Sample discretization is specified

		%[r_each_object{ii}, N_each_object(ii), delta_V_each_object{ii}, L_sub_each_object{ii}, volume(ii)] = read_sample_discretization(discretization{ii}, L_char(ii));
        [r_each_object{ii}, N_each_object(ii), delta_V_each_object{ii}, L_sub_each_object{ii}, volume(ii)] = read_sample_discretization(geometry{ii}, L_char(ii),discFile(ii));

	    else % User-defined discretization is specified

		[r_each_object{ii}, N_each_object(ii), delta_V_each_object{ii}, L_sub_each_object{ii}] = read_user_discretization(discretization{ii}, L_char(ii));

	    end % End if sample or user-defined discretization

	    if ii ~= 1
	    	ind_bulk(ii) = sum(N_each_object(1:ii-1))+1;
	    end

	    % Move the center-of-mass of each object to the origin [0,0,0]
	    r_each_object{ii} = center_of_mass(r_each_object{ii});

	    % Move each discretization to its user-specified origin
	    r_each_object{ii} = r_each_object{ii} + repmat(origin(ii,:), N_each_object(ii), 1);

	end % End loop through bulk objects

end
