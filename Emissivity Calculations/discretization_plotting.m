function discretization_plotting(r, L_sub, N, show_axes, output, saveDir,ind_bulk,r_1,r_2)
% Plot the voxel discretization and the normal discretization lattice
%	
%	Inputs:
%		r - NEED HELP
%		L_sub - NEED HELP
%		N - NEED HELP
%		show_axes - show the axes on the graphs
%		output - struct which contains wether the figures should be saved
%		saveDir - where to store the figures
%
	
	% Plot voxel image of discretization
	FIG_voxel = figure(1);
	%[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'on', [], [] );
 	[vert, fac] = voxel_image( r_1, L_sub(ind_bulk(1)), [], [], [], [], 'on', [], [] );
    	[vert, fac] = voxel_image( r_2, L_sub(ind_bulk(2)), [], [], [], [], 'on', [], [] );
	xlabel('x-axis (m)');
	ylabel('y-axis (m)');
	zlabel('z-axis (m)');
	if ~show_axes
	    grid off
	    axis off
	    set(gca, 'fontsize', 30)
	end
	view(5,20)

	% Visualize discretized lattice
	FIG_discretization = figure(2);
	plot3(r(:,1), r(:,2), r(:,3), 'x')
	title(['Location of each subvolume for N = ', num2str(N), ' total subvolumes'], 'fontsize', 14)
	xlabel('x-axis [m]')
	ylabel('y-axis [m]')
	zlabel('z-axis [m]')
	grid on
	if ~show_axes
	    grid off
	    axis off
	    set(gca, 'fontsize', 30)
	end
	set(gca,'DataAspectRatio',[1 1 1])
	view(3)


	% Save discretization figure file
	if output.save_fig
	    fig_path_dielectricFunction = [saveDir '/voxel'];
	    saveas(FIG_voxel, fig_path_dielectricFunction, string(output.figure_format))
	    clear FIG_voxel % Remove previous plot handles

	    fig_path_discretization = [saveDir '/discretization'];
	    saveas(FIG_discretization, fig_path_discretization, string(output.figure_format))
	    clear FIG_discretization % Remove previous plot handles
	end

end
