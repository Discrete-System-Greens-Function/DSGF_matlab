function [N_each_object, volume, r_each_object, ind_bulk, delta_V_each_object, L_sub_each_object] = read_discretization(discretization, L_char, origin)
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

		[r_each_object{ii}, N_each_object(ii), delta_V_each_object{ii}, L_sub_each_object{ii}, volume(ii)] = read_sample_discretization(discretization{ii}, L_char(ii));

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
