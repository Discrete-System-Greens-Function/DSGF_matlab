function dielectric_function_plotting(omega, epsilon, material, N_omega, output, saveDir)
% Plots the results from the dielectric function
% 
%	Inputs:
%		omega - NEED HELP 
%		epsilon - NEED HELP
%		material - which material is being plotted
%		N_omega - NEED HELP
%		output - struct containing whether the figure should be saved
%		saveDir - where to store the figure


	FIG_dielectric_function = figure(3);
	plot(omega, real(epsilon), 'x-', omega, imag(epsilon), 'o--', 'linewidth', 2)
	xlabel('Frequency [rad/s]', 'fontsize', 12)
	ylabel('Dielectric function, \epsilon', 'fontsize', 12)
	title(['Dielectric function of ' string(material) ', N_o_m_e_g_a = ' num2str(N_omega)], 'fontsize', 16)
	legend('Real part', 'Imaginary part', 'location', 'best')
	set(gca, 'fontsize', 16)
	axis tight
	grid on

	% Save dielectric function figure file
	if output.save_fig
	    fig_path_dielectricFunction = [saveDir '/dielectric_function.fig'];
	    saveas(FIG_dielectric_function, fig_path_dielectricFunction)
	    clear FIG_dielectric_function % Remove previous plot handles
	end
