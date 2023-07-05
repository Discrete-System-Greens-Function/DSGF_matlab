function subvol_heatmap_plotting_bulk_object_user_defined(r, L_sub, Q_total_subvol,  Q_density_subvol, c_limits, c_limits_density, abs_limit_density, show_axes, output, saveDir, N, N1)
%subvol_heatmap_user_defined

% plots the subvolume heatmap in 2 different views
%	
%	Inputs:
%		r - 
%		L_sub - 
%		Q_total_subvol -
%		c_limits -
%		show_axes - show the axes in the graphs
%		output - save the figure or not

    N_ref = 1;
    for i = 1:N1 
       if L_sub(i+1) ~= L_sub(i)
           N_ref = i;
       end    
    end 

	% Subvolume heat map for full particles (VIEW 1)
	heatmap_view_1 = figure(4);
    axis_y = 'Heat dissipated in a subvolume, Q [W]';
	%subplot(1,2,1)
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	
    %[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
        
	[vert, fac] = voxel_image( r(1:N_ref,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(:,1:N_ref), c_limits, axis_y );
    [vert, fac] = voxel_image( r(N_ref+1:N1,:), L_sub(N1), [], [], [], [], 'heatmap', Q_total_subvol(:,N_ref+1:N1), c_limits, axis_y );
    [vert, fac] = voxel_image( r(N1+1:N-N_ref,:), L_sub(N1+1), [], [], [], [], 'heatmap', Q_total_subvol(:,N1+1:N-N_ref), c_limits, axis_y );
    [vert, fac] = voxel_image( r(N-N_ref+1:N,:), L_sub(N), [], [], [], [], 'heatmap', Q_total_subvol(:,N-N_ref+1:N), c_limits, axis_y );
    
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
    
    
    norm_Q_density_subvol = Q_density_subvol'./abs_limit_density;

	% Subvolume heat map for full particles (VIEW 2)
	heatmap_view_2 = figure(5);
    axis_y = 'Normalized heat density in a subvolume';
	%subplot(1,2,2)
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
    
    %[vert, fac] = voxel_image( r, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol.', c_limits );
    
    [vert, fac] = voxel_image( r(1:N_ref,:), L_sub(1), [], [], [], [], 'heatmap', norm_Q_density_subvol(:,1:N_ref), c_limits_density, axis_y );
    [vert, fac] = voxel_image( r(N_ref+1:N1,:), L_sub(N1), [], [], [], [], 'heatmap', norm_Q_density_subvol(:,N_ref+1:N1), c_limits_density, axis_y );
    [vert, fac] = voxel_image( r(N1+1:N-N_ref,:), L_sub(N1+1), [], [], [], [], 'heatmap', norm_Q_density_subvol(:,N1+1:N-N_ref), c_limits_density, axis_y );
    [vert, fac] = voxel_image( r(N-N_ref+1:N,:), L_sub(N), [], [], [], [], 'heatmap', norm_Q_density_subvol(:,N-N_ref+1:N), c_limits_density, axis_y );
        
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
		
	    heatmap_view_1_path = [saveDir '/heatmap_full_power'];
	    heatmap_view_2_path = [saveDir '/heatmap_full_density'];

	    saveas(heatmap_view_1, heatmap_view_1_path, string(output.figure_format))
	    saveas(heatmap_view_2, heatmap_view_2_path, string(output.figure_format))

	end

end
