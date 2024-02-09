function subvol_heatmap(r, L_sub, Q_total_subvol, Q_density_subvol, c_limits, c_limits_density, abs_limit_density, show_axes, output, saveDir,ind_bulk,r_1,r_2,N)
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
    axis_y = 'Heat dissipated in a subvolume, Q [W]';
	%subplot(1,2,1)
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	%[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
    
    [vert, fac] = voxel_image( r_1, L_sub(ind_bulk(1)), [], [], [], [], 'heatmap', Q_total_subvol(ind_bulk(1):ind_bulk(2)-1).', c_limits, axis_y);
    [vert, fac] = voxel_image( r_2, L_sub(ind_bulk(2)), [], [], [], [], 'heatmap', Q_total_subvol(ind_bulk(2):N).', c_limits, axis_y); 
	xlabel('x-axis (m)');
	ylabel('y-axis (m)');
	zlabel('z-axis (m)');
	if ~show_axes
	    grid off
	    axis off
	    colorbar off
    end

    norm_Q_density_subvol = Q_density_subvol./abs_limit_density;
    
	% Subvolume heat map for full particles (VIEW 2)
	heatmap_view_2 = figure(5);
    axis_y = 'Normalized heat density [W/m^3] in a subvolume';
	%subplot(1,2,2)
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	%[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
    
    [vert, fac] = voxel_image( r_1, L_sub(ind_bulk(1)), [], [], [], [], 'heatmap', norm_Q_density_subvol(ind_bulk(1):ind_bulk(2)-1), c_limits_density, axis_y );
    [vert, fac] = voxel_image( r_2, L_sub(ind_bulk(2)), [], [], [], [], 'heatmap', norm_Q_density_subvol(ind_bulk(2):N), c_limits_density, axis_y );
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
    %view(2)
    %dim = [.6 .5 .3 .3];
    %str = ['Q_{Max} = '  num2str(abs_limit_density) 'W/m^3'];
    %annotation('textbox',dim,'String',str,'FitBoxToText','on');
	view(35,20)
	
	if output.save_fig
		
	    heatmap_view_1_path = [saveDir '/heatmap_full_power'];
	    heatmap_view_2_path = [saveDir '/heatmap_full_density'];

	    saveas(heatmap_view_1, heatmap_view_1_path, string(output.figure_format))
	    saveas(heatmap_view_2, heatmap_view_2_path, string(output.figure_format))

	end

end
