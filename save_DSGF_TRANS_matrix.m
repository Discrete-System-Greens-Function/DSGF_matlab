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
        fileName_DSGF = ['DSGFmatrix_omega' num2str(index_omega)];
	save_matrix(fileName_DSGF, DSGF_matrix, [filePath_st.DSGF, '/', fileName_DSGF, '.csv'], 'DSGF')
    end

    if output.transmission_coefficient_matrix
        % File name where results will be saved (based on what frequency band is chosen)
        fileName_Trans = ['TransMatrix_omega' num2str(index_omega)];

	save_matrix(fileName_Trans, Trans_coeff_matrix, [filePath_st.Trans, '/', fileName_Trans, '.csv'], 'transmission coefficient')
    end

end
