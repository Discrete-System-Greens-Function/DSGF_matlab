function subvol_heatmap_plotting(r, L_sub, Q_total_subvol, show_axes, output, saveDir)
% plots the heatmaps for the total heat dissipation per subvolume
% plots generated:
% 	heatmap for bulk object (2 different views)
% 	heatmap for slices of particles (XY & XZ)
% 	heatmap for half particles (XY & XZ) %% this might need correction %%


	[r_st, ind_st, Q_total_subvol_st] = planar_cut(r, L_sub, Q_total_subvol);


	% Set heatmap color axis limits
	abs_limit = max(abs(Q_total_subvol));
	c_limits = [-abs_limit, abs_limit];
	%c_limits = [min(Q_total_subvol), max(Q_total_subvol)];


	% plotting the subvolume heatmap in 2 seperate views
	subvol_heatmap_plotting_bulk_object(r, L_sub, Q_total_subvol, c_limits, show_axes, output, saveDir);


	subvol_heatmap_plotting_slice(r_st, L_sub, Q_total_subvol_st, c_limits, saveDir, output);


	subvol_heatmap_plotting_half(r_st, L_sub, Q_total_subvol_st, c_limits, saveDir, output);

end
