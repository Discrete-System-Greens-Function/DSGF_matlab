function [omega] = non_uniform_omega(material)
switch (material)
	case Material.SiO_2

	    omega_i = 7.5e13;
        N_ref_1 = 20;
        omega_ref_1 = 1.e14;
        N_ref_2 = 20; 
        omega_ref_2 = 1.5e14 ;
        N_ref_3 = 20; 
        omega_ref_3 = 2.e14 ;
        N_ref_4 = 30; 
        omega_ref_4 = 2.5e14 ;
        N_ref_5 = 10; 
        omega_ref_5 = 3e14;
        N_ref_6 = 5;
        omega_f = 3.5e14 ;

        N_ref_total = N_ref_1+N_ref_2+N_ref_3+N_ref_4+N_ref_5+N_ref_6;
        omega_range_ref_1 = linspace(omega_i, omega_ref_1, N_ref_1);
        omega_range_ref_2 = linspace(omega_ref_1, omega_ref_2, N_ref_2);
        omega_range_ref_3 = linspace(omega_ref_2, omega_ref_3, N_ref_3);
        omega_range_ref_4 = linspace(omega_ref_3, omega_ref_4, N_ref_4);
        omega_range_ref_5 = linspace(omega_ref_4, omega_ref_5, N_ref_5);
        omega_range_ref_6 = linspace(omega_ref_5, omega_f, N_ref_6);


        omega = [omega_range_ref_1,omega_range_ref_2,omega_range_ref_3,omega_range_ref_4,omega_range_ref_5,omega_range_ref_6];
        omega = unique (omega');
        N_omega = length(omega);

	case Material.SiC

	    omega_i = 1.4e14;
        N_ref_1 = 5;
        omega_ref_1 = 1.43e14;
        N_ref_2 = 60; 
        omega_ref_2 = 1.55e14 ;
        N_ref_3 = 25; %15, 10
        omega_ref_3 = 1.76e14 ;
        N_ref_4 = 5; %10
        omega_ref_4 = 1.8e14 ;
        N_ref_5 = 5; %10
        omega_ref_5 = 1.84e14;
        N_ref_6 = 5; %10
        omega_f = 1.9e14 ;

        N_ref_total = N_ref_1+N_ref_2+N_ref_3+N_ref_4+N_ref_5+N_ref_6;
        omega_range_ref_1 = linspace(omega_i, omega_ref_1, N_ref_1);
        omega_range_ref_2 = linspace(omega_ref_1, omega_ref_2, N_ref_2);
        omega_range_ref_3 = linspace(omega_ref_2, omega_ref_3, N_ref_3);
        omega_range_ref_4 = linspace(omega_ref_3, omega_ref_4, N_ref_4);
        omega_range_ref_5 = linspace(omega_ref_4, omega_ref_5, N_ref_5);
        omega_range_ref_6 = linspace(omega_ref_5, omega_f, N_ref_6);


        omega = [omega_range_ref_1,omega_range_ref_2,omega_range_ref_3,omega_range_ref_4,omega_range_ref_5,omega_range_ref_6];
        omega = unique (omega');
        N_omega = length(omega);


	case Material.SiN

	    omega_i = 0.2e14;%7.5e13;
        N_ref_1 = 5;
        omega_ref_1 = 0.5e14;
        N_ref_2 = 20; 
        omega_ref_2 = 1.e14 ;
        N_ref_3 = 20; 
        omega_ref_3 = 1.5e14 ;
        N_ref_4 = 20; 
        omega_ref_4 = 2.e14 ;
        N_ref_5 = 20; 
        omega_ref_5 = 2.5e14;
        N_ref_6 = 5;
        omega_f = 3.e14 ;

        N_ref_total = N_ref_1+N_ref_2+N_ref_3+N_ref_4+N_ref_5+N_ref_6;
        omega_range_ref_1 = linspace(omega_i, omega_ref_1, N_ref_1);
        omega_range_ref_2 = linspace(omega_ref_1, omega_ref_2, N_ref_2);
        omega_range_ref_3 = linspace(omega_ref_2, omega_ref_3, N_ref_3);
        omega_range_ref_4 = linspace(omega_ref_3, omega_ref_4, N_ref_4);
        omega_range_ref_5 = linspace(omega_ref_4, omega_ref_5, N_ref_5);
        omega_range_ref_6 = linspace(omega_ref_5, omega_f, N_ref_6);


        omega = [omega_range_ref_1,omega_range_ref_2,omega_range_ref_3,omega_range_ref_4,omega_range_ref_5,omega_range_ref_6];
        omega = unique (omega');
        N_omega = length(omega);

end

