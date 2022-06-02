function save_matrix(filename, matrix, filepath, name)
% function to save a matrix and display how long it took
%
%	Inputs:
%		filename - name of the file to save the matrix to
%		matrix - matrix to save
%		filepath - filepath to save it to
%		name - name for display purposes


	% File name where results will be saved (based on what frequency band is chosen)
	t1 = toc;
	% Export DSGF matrix
	writematrix(matrix, filepath);
	t2 = toc;
	disp(['Time to save ' name ' matrix as a .csv file  = ' num2str(t2-t1) ' s = ' num2str((t2-t1)/60) ' minutes'])

end
