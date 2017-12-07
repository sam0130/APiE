clc; clear variables; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 2;                                                                     % mass   
k = 5;                                                                     % spring   
gamma = 0.3;                                                               % friction parameter
f = 3*m;                                                                   % force amplitude
omg = (k/m)^0.5;                                                           % frequency
%omg_f = 0.9*omg;                                                           % force frequency
omg_f = 0.1*omg:0.2:10*omg;                                                % force freq vector for amplitude analysis

T_period = 2*pi/omg;
n = 30;                                                                    % No. of cycles
T = n*T_period;                                                            % Total time
N = 10000;                                                                 % No. of time steps    
DeltaT = T/N;                                                               
t = 0:DeltaT:T;
cycles = ceil(t/T_period);                                      
cycles_steady = find(cycles>15);
%% Exact Solution
x_anal0 = 1;
v_anal0 = 1;
A_steady = zeros(length(omg_f),1);
for count = 1:length(omg_f)
[x_trans,x_steady,x_anal,v_anal,A_steady(count), A_total]=...
          exactSolutionFn(x_anal0,v_anal0,t, m, k, gamma, f, omg_f(count));
end                                 
%% Euler Algorithm
x_euler = zeros(length(t),1);         
v = zeros(length(t),1);
x_euler(1) = 1;                       
v(1) = 1;              
A_Euler = zeros(length(omg_f),1);
for count = 1:length(omg_f)
    for i = 1:(length(t)-1)
        x_euler(i+1) = x_euler(i) + DeltaT*v(i);
        v(i+1) = v(i) + (-k*x_euler(i))/m * DeltaT + (-gamma*v(i)*DeltaT/m) ...
            + f*cos(omg_f(count)*i*DeltaT)*DeltaT/m;
    end
    A_Euler(count) = max(x_euler(cycles_steady));
end
[ kin_exact, pot_exact, tot_exact] = calcEnergy(x_anal,v_anal,m,k);
[ kin_euler, pot_euler, tot_euler ] = calcEnergy(x_euler,v,m,k);
%% Verlet Leapfrog Algorithm
% w leap-frogs
% v_leap corrects the w by taking average of -1/2 and +1/2 velocity

x_leap = zeros(length(t),1); 
w = zeros(length(t),1);                                                    % 1/2 velocities
v_leap = zeros(length(t),1);                                               % Corrected velocity 

x_leap(1) = 1;       
v_leap(1) = 1;
w(1) = (v_leap(1) - 0.5*DeltaT*(-k*x_leap(1))/m - ...
     DeltaT*0.5*f/m*cos(omg_f(count)*0.5*DeltaT))/(1-((gamma/m)*DeltaT/2));% Euler step for -1/2 velocity        
    
%%% the velocity uses the normal points' displacements for the force term
%%% while staggered points are being used for the friction term
A_Verlet = zeros(length(omg_f),1);
for count = 1:length(omg_f)
    for i = 1:(length(t)-1)
        w(i+1) = w(i) + DeltaT*(-k*x_leap(i))/m - gamma*w(i)*DeltaT/m + ...
            DeltaT*(f/m)*cos(((i+1)*DeltaT*omg_f(count)));
        x_leap(i+1) = x_leap(i) + DeltaT*w(i+1);
        v_leap(i) = 0.5 * (w(i+1) + w(i));
    end
    A_Verlet(count) = max(x_leap(cycles_steady));
end
v_leap(end) = w(end) + 0.5*DeltaT*(-k*x_leap(end))/m - gamma*w(end)*0.5*...
          DeltaT/m + 0.5*DeltaT*(f/m)*cos(((i+1)*0.5*DeltaT*omg_f(count))); 

[ kin_verlet_w, pot_verlet_w, tot_verlet_w] = calcEnergy(x_leap, w, m,k);       % Energies using the leap frogging velocities
[ kin_verlet, pot_verlet, tot_verlet] = calcEnergy(x_leap, v_leap,m,k);        % Energies using corrected velocities
%% ode45
x0=[1 1];                                                                  % initial conditions
tspan=[0; T];
A_ode45 = zeros(length(omg_f),1);
for count = 1:length(omg_f)
    odefun = @(t,x) [x(2);
        -k/m * x(1) - gamma/m * x(2) + f/m*cos(omg_f(count)*t)];
    [t_ode,xy_ode45]=ode45(odefun, 0:DeltaT:tspan(2), x0);
    A_ode45(count) = max(xy_ode45(cycles_steady,2));    
end
[kin_ode45,pot_ode45,tot_ode45]=calcEnergy(xy_ode45(:,1),xy_ode45(:,2),m,k); 
%% saving data
filename = 'plotting.mat';
save(filename)