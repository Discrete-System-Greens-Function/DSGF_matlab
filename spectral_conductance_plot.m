function spectral_conductance_plot(omega, G_omega_bulk_12, output, saveDir)
% plots the spectral conductance against the frequency
%
% 	Inputs:
%		omega - ?
%		G_omega_bulk_12 - ?
%		output - struct containing whether the figure should be saved
%		saveDir - location for where the figure should be saved

    % Plot spectral conductance vs. frequency
    spectral_conductance_fig = figure(10);
    semilogy(omega, G_omega_bulk_12, '-x', 'linewidth', 2)
    %plot(E_eV, G_omega_eV.', '-x', 'linewidth', 2)
    xlabel('Frequency [rad/s]')
    ylabel('Spectral conductance, G_\omega  [WK^-^1(rad/s)^-^1]', 'fontsize', 12)
%     title([material ', ' geometry ', R_1 = ' num2str(radius_1*(10^9)) 'nm, R_2 = ' num2str(radius_2*(10^9)) 'nm, d_c = '...
%         num2str(d_center*(10^9)) 'nm, ' num2str(N_omega) ' frequencies, N = ' num2str(N) ' total subvolumes'], 'fontsize', 16)
    axis tight
    set(gca, 'fontsize', 22)
    %legend([num2str(N(1)/2) ' subvolumes'], 'location', 'best')
    grid on

    % Save figure files
    if output.save_fig
        spectral_conductance_fig_path = [saveDir '/spectralConductance.fig'];
        saveas(spectral_conductance_fig, spectral_conductance_fig_path)
    end

end
