function [r_st, ind_st, Q_total_subvol_st] = planar_cut(r, L_sub, Q_total_subvol)
% calculates the indices and total heat dissipation for half & slice on XY-PLANE
% calculates the indices and total heat dissipation for half & slice on XZ-PLANE
%
%	Inputs:
%		r - 
%		L_sub - 
%		Q_total_subvol -
%
%	Outputs:
%	
%	all the outputs have 4 fields:
%	half.xy | half.xz | slice.xy | slice.xz
%
%		r_st - 
%		ind_st - contains all the indices for the slices 
%		Q_total_subvol_st - Q_total_subvol indexed using the respective ind_st

	r_st = struct();
	ind_st = struct();
	q_total_subvol_st = struct();

	% XY-PLANE CUT: Find locations and indices for particle halves cut down the middle (cut along the xy-plane)
	r_z_avg = (max(r(:,3)) + min(r(:,3)))/2;
	ind_st.half.xy = find(r(:,3) >= r_z_avg);
	r_st.half.xy = r(ind_st.half.xy, :);
	Q_total_subvol_st.half.xy = Q_total_subvol(ind_st.half.xy);

	% XZ-PLANE CUT: Find locations and indices for particle halves cut down the middle (cut along the xz-plane)
	r_y_avg = (max(r(:,2)) + min(r(:,2)))/2;
	ind_st.half.xz = find(r(:,2) >= r_y_avg);
	r_st.half.xz = r(ind_st.half.xz, :);
	Q_total_subvol_st.half.xz = Q_total_subvol(ind_st.half.xz);

	% XY-PLANE CUT: Find locations and indices of a single slice (one subvolume thick)
	ind_st.slice.xy = find((r(:,3) >= r_z_avg) & (r(:,3) <= r_z_avg + L_sub(1)));
	r_st.slice.xy = r(ind_st.slice.xy, :);
	Q_total_subvol_st.slice.xy = Q_total_subvol(ind_st.slice.xy);

	% XZ-PLANE CUT: Find locations and indices of a single slice (one subvolume thick)
	ind_st.slice.xz = find((r(:,2) >= r_y_avg) & (r(:,2) <= r_y_avg + L_sub(1)));
	r_st.slice.xz = r(ind_st.slice.xz, :);
	Q_total_subvol_st.slice.xz = Q_total_subvol(ind_st.slice.xz);

end
