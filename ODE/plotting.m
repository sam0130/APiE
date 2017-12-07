clc; close all; clear variables;
load('plotting.mat')
%% Euler %%%
figure()
hold on
plot(t/T_period,pot_euler);
plot(t/T_period,kin_euler);
plot(t/T_period,tot_euler);
plot(t/T_period,tot_exact);
legend('PE_{num}','KE_{num}','Tot_{num}','Tot_{exact}')
xlabel('No. of cycles')
ylabel('Energy (J)')
title('Euler vs Exact')
%% Verlet %%%
figure()
hold on
plot(t/T_period,pot_verlet);
plot(t/T_period,kin_verlet);
plot(t/T_period,tot_verlet);
plot(t/T_period,tot_exact);
legend('PE_{num}','KE_{num}','Tot_{num}','Tot_{exact}')
xlabel('No. of cycles')
ylabel('Energy (J)')
title('Verlet Leapfrog (corrected) vs Exact')

% energy using corrected velocities----------------------------------------
figure()
hold on
plot(t/T_period,tot_verlet)
plot(t/T_period,tot_verlet_w,'k--')
plot(t/T_period, tot_exact)
legend('Verlet Corrected','Verlet LeapFrogging','Exact')
xlabel('No. of cycles')
ylabel('Energy (J)')
title('Comparison of corrected and normal Verlet leapfrog')
%% ode45 %%Exact
figure()
xlabel('Time (s)')
ylabel('Energy (J)')
hold on
plot(t_ode/T_period,pot_ode45);
plot(t_ode/T_period,kin_ode45);
plot(t_ode/T_period,tot_ode45);
plot(t/T_period,tot_exact);
legend('PE_{num}','KE_{num}','Tot_{num}','Tot_{exact}')
xlabel('No. of cycles')
ylabel('Energy (J)')
title('ode45 Vs Exact')
%% comparing displacements
figure()
hold on
plot(t_ode/T_period, xy_ode45(:,1));
plot(t/T_period,x_euler)
plot(t/T_period,x_leap)
plot(t/T_period, x_anal);
xlabel('No. of cycles')
ylabel('Displacement (m)')
title('Displacement Comparison')
legend('ode45', 'Euler', 'Verlet', 'Exact')
%% comparing energies
figure()
hold on
plot(t/T_period,tot_euler)
plot(t_ode/T_period,tot_ode45);
plot(t/T_period,tot_verlet)
plot(t/T_period, tot_exact)
xlabel('No. of cycles')
ylabel('Energy (J)')
legend('Euler','ode45','Verlet Leapfrog','Exact')
title('Comparison of total energies by different methods')
%% amplitude vs omega
figure() 
loglog(omg_f,A_steady/f)
hold on
loglog(omg_f,A_Euler/f)
loglog(omg_f, A_Verlet/f)
loglog(omg_f, A_ode45/f)
xlabel('\omega_f')
ylabel('A/f')
legend('Exact','Euler','Verlet','ode45')
title('Amplitude variation with driving frequency')
