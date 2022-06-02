function [filePath_st, file_name_saved] = result_setup(output, description)
% function that sets up the results folder
%
% 	Inputs:
% 		output - struct containing what things will be saved
% 		description - description to include in the filename for the .mat variables
%	
%	Outputs:
% 		filePath_st - struct with 3 fields (main, DSGF, Trans)
%		file_name_saved - filename for the .mat variables that get saved

	% Set time stamp
	time_stamp = datestr(now,'yyyy-mm-dd_HH-MM-SS');

	% Directory where results files will be saved (string)
	saveDir = sprintf('Results_%s', time_stamp);  % Folder name.
	if not(isfolder(saveDir)) % If the results directory doesn't exist, create it.
	    mkdir(saveDir)
	end

	% Make new directory where DSGF matrix files will be saved
	if output.DSGF_matrix
	    mkdir(saveDir, 'DSGF_matrices')
	end

	% Make new directory where transmission coefficient matrix files will be saved
	if output.transmission_coefficient_matrix
	    mkdir(saveDir, 'Trans_matrices')
	end

	% File paths where different results will be saved
	filePath_st.main = saveDir;
	filePath_st.DSGF = [saveDir '/' 'DSGF_matrices'];
	filePath_st.Trans = [saveDir '/' 'Trans_matrices'];

	% File name of .mat saved variables
	file_name_saved = ['results_' description '_' time_stamp]; % File name where results will be saved

end
