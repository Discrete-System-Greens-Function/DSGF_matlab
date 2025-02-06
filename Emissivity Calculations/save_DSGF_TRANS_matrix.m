function save_DSGF_TRANS_matrix(output, filePath_st, index_omega, DSGF_matrix, Trans_coeff_matrix)
% saves the DGSF and the transmission coeff matrix
%
%	Inputs:
%		output - struct containing which of the 2 matrices should be saved
%		filePath_st - struct containing the file path for the matrices
%		index_omega - index denoting at which angular frequency the matrices are calculated
%		DSGF_matrix - the DSGF_matrix
%		Trans_coeff_matrix - the matrix for the transmission coefficient
%


    if output.DSGF_matrix
        % File name where results will be saved (based on what frequency band is chosen)
        fileName_DSGF = [filePath_st.DSGF '/DSGF_matrix_omega_' num2str(index_omega) '.csv'];
	    %save_matrix(fileName_DSGF, DSGF_matrix, [filePath_st.DSGF, '/',fileName_DSGF, '.csv'], 'DSGF') % no label on rows and columns
        save_matrix(fileName_DSGF, DSGF_matrix, fileName_DSGF, 'DSGF') % no label on rows and columns
    end

    if output.transmission_coefficient_matrix
        % File name where results will be saved (based on what frequency band is chosen)
        fileName_Trans = [filePath_st.Trans '/TransMatrix_omega_' num2str(index_omega) '.csv'];
        %save_matrix(fileName_Trans, Trans_coeff_matrix, [filePath_st.Trans, '/', fileName_Trans, '.csv'], 'transmission coefficient')
        columnNames = arrayfun(@(x) sprintf(' j = %d', x), 1:length(Trans_coeff_matrix), 'UniformOutput', false); % Define column names (matching the number of columns)
        rowNames = arrayfun(@(x) sprintf('i = %d', x), 1:length(Trans_coeff_matrix), 'UniformOutput', false);% Define row names (matching the number of rows)
        Transtable = array2table(Trans_coeff_matrix,'VariableNames', columnNames);% Create a table with column names
        Transtable.Properties.RowNames = rowNames;
        writetable(Transtable,fileName_Trans, 'WriteRowNames', true);
    end

end
