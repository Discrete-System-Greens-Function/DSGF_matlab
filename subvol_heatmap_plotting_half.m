function subvol_heatmap_plotting_half(r_st, L_sub, Q_total_subvol_st, c_limits, saveDir, output)
% plots the slice of the particles on the XY & XZ plane
% 
%	Inputs:
%		r_st - struct containing the r's for the slices
%		L_sub -
%		Q_total_subvol_st - struct containing the total heat dissipation in the slices
%		c_limits - 
%		saveDir - where to save the figures
%		output - struct containing wheter the figure should be saved
%

	
	% XY-PLANE CUT: Subvolume heat map for half particles
	heatmap_xy_half_fig = figure(6);
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	[vert, fac] = voxel_image( r_st.half.xy, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol_st.half.xy.', c_limits );
	%view(2)

	% XZ-PLANE CUT: Subvolume heat map for half particles
	heatmap_xz_half_fig = figure(8);
	%[vert, fac] = voxel_image( r(1:N1,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(1:N1).' ); % Absorber (T = 0 K)
	%[vert, fac] = voxel_image( r(N1+1:end,:), L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol(N1+1:end).' ); % Emitter (T = 300 K)
	[vert, fac] = voxel_image( r_st.half.xz, L_sub(1), [], [], [], [], 'heatmap', Q_total_subvol_st.half.xz.', c_limits );
	%view(2)

	% Save figure files
	if output.save_fig
	    heatmap_xy_half_path = [saveDir '/heatmap_xy_half.fig'];
	    heatmap_xz_half_path = [saveDir '/heatmap_xz_half.fig'];

	    saveas(heatmap_xy_half_fig, heatmap_xy_half_path)
	    saveas(heatmap_xz_half_fig, heatmap_xz_half_path)

	end

end
