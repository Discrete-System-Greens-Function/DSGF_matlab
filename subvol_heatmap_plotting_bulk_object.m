function subvol_heatmap(r, L_sub, Q_total_subvol, c_limits, show_axes, output, saveDir)
% plots the subvolume heatmap in 2 different views
%	
%	Inputs:
%		r - 
%		L_sub - 
%		Q_total_subvol -
%		c_limits -
%		show_axes - show the axes in the graphs
%		output - save the figure or not

	% Subvolume heat map for full particles (VIEW 1)
	heatmap_view_1 = figure(4);
	%subplot(1,2,1)
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
	xlabel('x-axis (m)');
	ylabel('y-axis (m)');
	zlabel('z-axis (m)');
	if ~show_axes
	    grid off
	    axis off
	    colorbar off
	end
	%view(2)
	view(-30,35)

	% Subvolume heat map for full particles (VIEW 2)
	heatmap_view_2 = figure(5);
	%subplot(1,2,2)
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
	xlabel('x-axis (m)');
	ylabel('y-axis (m)');
	zlabel('z-axis (m)');
	if ~show_axes
	    grid off
	    axis off
	    %cb = colorbar;
	    %colorbar('east')
	    %set(cb,'position',[0.2 0.2 .05 .5]) % [xposition yposition width height].
	    set(gca, 'fontsize', 30)
	end
	view(35,20)
	
	if output.save_fig
		
	    heatmap_view_1_path = [saveDir '/heatmap_full_view1'];
	    heatmap_view_2_path = [saveDir '/heatmap_full_view2'];

	    saveas(heatmap_view_1, heatmap_view_1_path, string(output.figure_format))
	    saveas(heatmap_view_2, heatmap_view_2_path, string(output.figure_format))

	end

end
