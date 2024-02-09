function subvol_heatmap_plotting_user_defined(r, L_sub, Q_total_subvol, Q_density_subvol, show_axes, output, saveDir,  N, N1)
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
    
    
    abs_limit_density = max(abs(Q_density_subvol));
    c_limits_density = [-abs_limit_density/abs_limit_density, abs_limit_density/abs_limit_density];


	% plotting the subvolume heatmap in 2 seperate views
	subvol_heatmap_plotting_bulk_object_user_defined(r, L_sub, Q_total_subvol, Q_density_subvol, c_limits, c_limits_density, abs_limit_density, show_axes, output, saveDir, N, N1);
    %subvol_heatmap_plotting_bulk_object(r, L_sub, Q_total_subvol, Q_density_subvol, c_limits, c_limits_density, abs_limit_density, show_axes, output, saveDir,ind_bulk,r_1,r_2,N);


    
    if output.heatmap_sliced
        subvol_heatmap_plotting_slice(r_st, L_sub, Q_total_subvol_st, c_limits, saveDir, output);


        subvol_heatmap_plotting_half(r_st, L_sub, Q_total_subvol_st, c_limits, saveDir, output);
    end 
    
end
