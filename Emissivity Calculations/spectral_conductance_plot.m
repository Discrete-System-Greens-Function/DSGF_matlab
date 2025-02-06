function spectral_conductance_plot(omega, G_w_AB, T_cond, output, saveDir)
% plots the spectral conductance against the frequency
%
% 	Inputs:
%		omega - Radial frequency
%		G_w_AB - Spectral conductance
%		output - struct containing whether the figure should be saved
%		saveDir - location for where the figure should be saved

%{
    spectral_conductance_fig = figure(10);
    semilogy(omega, G_w_AB(:,3), '-x', 'linewidth', 2)
    xlabel('Frequency [rad/s]')
    ylabel('Spectral conductance, G_\omega  [WK^-^1(rad/s)^-^1]', 'fontsize', 12)
%     title([material ', ' geometry ', R_1 = ' num2str(radius_1*(10^9)) 'nm, R_2 = ' num2str(radius_2*(10^9)) 'nm, d_c = '...
%         num2str(d_center*(10^9)) 'nm, ' num2str(N_omega) ' frequencies, N = ' num2str(N) ' total subvolumes'], 'fontsize', 16)
    axis tight
    set(gca, 'fontsize', 22)
    tempstr = sprintf('T = %d K',T_cond(3));
    legend(tempstr , 'location', 'best')
    grid on
%}
    
    % Plot multiple spectral conductance vs. frequency
    multiple_spectral_conductance_fig = figure(11);
    tempstrmult = cell(length(T_cond),1);
    hold on
    for i = 1: length(T_cond)
        semilogy(omega, G_w_AB(:,i), '-x', 'linewidth', 2)
        tempstrmult{i} = sprintf('T = %d K',T_cond(i)); %strcat
       
    end
    %plot(E_eV, G_omega_eV.', '-x', 'linewidth', 2)
    xlabel('Frequency [rad/s]')
    ylabel('Spectral conductance, G_\omega  [WK^-^1(rad/s)^-^1]', 'fontsize', 12)
%     title([material ', ' geometry ', R_1 = ' num2str(radius_1*(10^9)) 'nm, R_2 = ' num2str(radius_2*(10^9)) 'nm, d_c = '...
%         num2str(d_center*(10^9)) 'nm, ' num2str(N_omega) ' frequencies, N = ' num2str(N) ' total subvolumes'], 'fontsize', 16)
    axis tight
    set(gca, 'fontsize', 22)
    %legend(string(T_cond(1)),string(T_cond(2)),string(T_cond(3)),string(T_cond(4)),string(T_cond(5)),'location', 'best')
    %tempstrmult = sprintf('T = %d K',T_cond(:));
    %legend(string(T_cond(:)))
    legend(tempstrmult , 'location', 'best')
    %legend(string(tempstrmult(:)))
    grid on
    hold off
    %}
    
    % Save figure files
    if output.save_fig
       % spectral_conductance_fig_path = [saveDir '/spectralConductance'];
       % saveas(spectral_conductance_fig, spectral_conductance_fig_path, string(output.figure_format))
        
        multiple_spectral_conductance_fig_path = [saveDir '/multipleSpectralConductance'];
        saveas(multiple_spectral_conductance_fig, multiple_spectral_conductance_fig_path, string(output.figure_format))
    end
    
end
