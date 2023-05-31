function discretization_plotting_user_defined(r, L_sub, N, show_axes, output, saveDir, N1)
% Plot the voxel discretization and the normal discretization lattice
%	
%	Inputs:
%		r - Discretized lattice including subvolumes of all objects in one matrix (N x 3 matrix)
%		L_sub - Length of side of a cubic subvolume
%		N - total number of subvolumes 
%		show_axes - show the axes on the graphs
%		output - struct which contains wether the figures should be saved
%		saveDir - where to store the figures
%
	N_ref = 1;
    for i = 1:N1 
       if L_sub(i+1) ~= L_sub(i)
           N_ref = i
       end    
    end 
    
	% Plot voxel image of discretization
	FIG_voxel = figure(1);
	%[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'on', [], [] );
    
	[vert, fac] = voxel_image( r(1:N_ref,:), L_sub(1), [], [], [], [], 'on', [], [] );
    [vert, fac] = voxel_image( r(N_ref+1:N1,:), L_sub(N1), [], [], [], [], 'on', [], [] );
    [vert, fac] = voxel_image( r(N1+1:N-N_ref,:), L_sub(N1+1), [], [], [], [], 'on', [], [] );
    [vert, fac] = voxel_image( r(N-N_ref+1:N,:), L_sub(N), [], [], [], [], 'on', [], [] );
    
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
