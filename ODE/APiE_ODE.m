clc; clear variables; close all;
%%%%%%%%%%%%%% Constants  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m k gamma f omega
m = 2;                              % mass;   
k = 5;                              % spring   
gamma = 6;                          % friction parameter
f = 0;                              % force amplitude
omega = 0;                          % force frequency

omg = (k/m)^0.5;                    % frequency
T_period = 2*pi/omg;
omg_1 = sqrt( (omg^2) - (( (gamma/m)^2)/4) );   % modified frequency for friction

n = 10;                             % No. of cycles
T = n*T_period;                     % Total time
N = 10000;                          % No. of time steps    
DeltaT = T/N;
t = 0:DeltaT:T;
%% Exact Solution
% initial conditions ==>> 
%x(0) = 1; v(0) = 1;
%A = sqrt(1 + 1/omg^2);                              % Amplitude
%phi = atan(-1/omg);
A = sqrt(1 + ((1 + gamma/2/m)^2)/(omg_1^2));        % modified amplitude with friction
phi = atan( - (1 + gamma/2/m)/omg_1 );              % modified phase factor with friction
%x_anal = @(t) A*cos(omg*t + phi);   
x_anal = @(t) A*exp((-gamma/2/m).*t).*cos(omg_1*t + phi);  % modified solution with friction
%v_anal = @(t) -A*omg*sin(omg*t + phi);
v_anal = @(t) A*( (-gamma/2/m) * exp((-gamma/2/m).*t).*cos(omg_1*t + phi) ...
                - omg_1 * sin(omg_1*t + phi) .* exp((-gamma/2/m).*t));        % modified solution with friction  
% syms x(t1) 
% Du = diff(x);
% cond1 = x(0) == 1;
% cond2 = Du(0) == 1;
% conds = [cond1 cond2];
% x_anal = dsolve(diff(x,t1,2) + gamma*(diff(x,t1,1)) + omg^2 * (x) == 0, conds);
% v_anal = diff(x_anal);
% x_anal = eval(subs(x_anal,t));
% v_anal = eval(subs(v_anal,t));
%% Euler Algorithm
x_euler = zeros(length(t),1);         
v = zeros(length(t),1);
x_euler(1) = 1;                       
v(1) = 1;              

for i = 1:(length(t)-1)
    x_euler(i+1) = x_euler(i) + DeltaT*v(i);
    v(i+1) = v(i) + (-k*x_euler(i))/m * DeltaT;
end

[ kin_exact, pot_exact, tot_exact] = calcEnergy(x_anal(t),v_anal(t));
[ kin_euler, pot_euler, tot_euler ] = calcEnergy(x_euler,v);

figure()
hold on
plot(t,x_euler)
plot(t, x_anal(t))
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Euler vs Exact')
legend('Numerical', 'Exact')

figure()
hold on
plot(t,pot_euler);
plot(t,kin_euler);
plot(t,tot_euler);
plot(t,tot_exact);
legend('PE_{num}','KE_{num}','Tot_{num}','Tot_{exact}')
xlabel('Time (s)')
ylabel('Energy (J)')
title('Euler vs Exact')
%% Verlet Leapfrog Algorithm
% w leap-frogs
% v_leap corrects the w by taking average of -1/2 and +1/2 velocity

x_leap = zeros(length(t),1); 
w = zeros(length(t),1);                      % 1/2 velocities
v_leap = zeros(length(t),1);                 % Corrected velocity 

x_leap(1) = 1;       
v_leap(1) = 1;
w(1) = (v_leap(1) - 0.5*DeltaT*(-k*x_leap(1))/m )/(1 - ((gamma/m)*DeltaT/2));     % Euler step for -1/2 velocity        
for i = 1:(length(t)-1)
    w(i+1) = w(i) + DeltaT*(-k*x_leap(i))/m - gamma*w(i)*DeltaT/m ;               % the velocity uses the normal points' displacements for the force term, while staggered points are being used for the friction term
    x_leap(i+1) = x_leap(i) + DeltaT*w(i+1);
    v_leap(i) = 0.5 * (w(i+1) + w(i));
end
v_leap(end) = w(end) + 0.5*DeltaT*(-k*(x_leap(end)))/m - gamma*w(i)*0.5*DeltaT/m;

[ kin_verlet_w, pot_verlet_w, tot_verlet_w] = calcEnergy(x_leap, w);        % Energies using the leap frogging velocities
[ kin_verlet, pot_verlet, tot_verlet] = calcEnergy(x_leap, v_leap);         % Energies using corrected velocities

figure()
hold on
plot(t,x_leap)
plot(t, x_anal(t))
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Verlet Leapfrog vs Exact')
legend('Numerical', 'Exact')

figure()
hold on
plot(t,pot_verlet);
plot(t,kin_verlet);
plot(t,tot_verlet);
plot(t,tot_exact);
legend('PE_{num}','KE_{num}','Tot_{num}','Tot_{exact}')
xlabel('Time (s)')
ylabel('Energy (J)')
title('Verlet Leapfrog (corrected) vs Exact')

figure()
hold on
plot(t,tot_verlet)
plot(t,tot_verlet_w,'k--')
plot(t, tot_exact)
legend('Verlet Corrected','Verlet LeapFrogging','Exact')
xlabel('Time (s)')
ylabel('Energy (J)')
title('Comparison of corrected and normal Verlet leapfrog')
xlabel('Time (s)')
ylabel('Energy (J)')
%% ode45
x0 = [1 1];                                            % initial conditions
tspan = [0; T];
[t_ode, xy_ode45] = ode45(@harmOscill, 0:0.05:tspan(2), x0);
figure()
hold on
plot(t_ode, xy_ode45(:,1));
plot(t,x_leap)
plot(t, x_anal(t));
title('ode45 vs Exact')
legend('ode45', 'Verlet', 'Exact')

[ kin_ode45, pot_ode45, tot_ode45] = calcEnergy(xy_ode45(:,1), xy_ode45(:,2)); 

figure()
xlabel('Time (s)')
ylabel('Energy (J)')
hold on
plot(t_ode,pot_ode45);
plot(t_ode,kin_ode45);
plot(t_ode,tot_ode45);
plot(t,tot_exact);
legend('PE_{num}','KE_{num}','Tot_{num}','Tot_{exact}')
xlabel('Time (s)')
ylabel('Energy (J)')
title('ode45 Vs Exact')
%% comparing energies
figure()
hold on
plot(t,tot_euler)
plot(t_ode,tot_ode45);
plot(t,tot_verlet)
plot(t, tot_exact)
xlabel('Time (s)')
ylabel('Energy (J')
legend('Euler','ode45','Verlet Leapfrog','Exact')
title('Comparison of total energies by different methods')
%% Force using ODE
x0 = [1 1];                             % initial conditions
tspan = [0; T];
[t_ode, xy_ode45] = ode45(@harmOscill_fric_force, tspan, x0);
figure()
plot(t_ode, xy_ode45(:,1));