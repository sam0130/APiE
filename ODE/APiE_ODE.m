clc; clear; close all;

m = 2; k = 5; gamma = 0.3;      % mass; spring constant; friction parameter
omg = (k/m)^0.5;                                                % frequency
T_period = 2*pi/omg;

n = 10;                                                      % No. of cycles
T = n*T_period;                                                % Total time
N = 10000;
DeltaT = T/N;
t = 0:DeltaT:T;

%---------------Initialise-------------------------------------------------
x = zeros(length(t),1);         v = zeros(length(t),1);
x(1) = 1;                       v(1) = 1;              % Initial Conditions
%--------------------------------------------------------------------------

%---------------Forward Euler----------------------------------------------
for i = 1:(length(t)-1)
    x(i+1) = x(i) + DeltaT*v(i);
    v(i+1) = v(i) + (-k*x(i))/m * DeltaT;
end
%--------------------------------------------------------------------------

x_anal = @(t) cos(omg*t) + (1/omg)*sin(omg*t);             % Exact Solution
v_anal = @(t) -omg*sin(omg*t) + cos(omg*t);
A = sqrt(1+1/omg^2);                                            % Amplitude

figure()
hold on
plot(t,x,'ko')
plot(t, x_anal(t))

pot = 0.5*k*x.^2;
kin = 0.5*m*v.^2;
tot_numerical = pot + kin;
total_euler = tot_numerical;
tot_exact = 0.5*k*A^2*ones(length(t),1)';

figure()
hold on
plot(t,pot);
plot(t,kin);
plot(t,tot_numerical);
plot(t,tot_exact);
legend('potential_{numerical}','kinetic_{numerical}','total_{numerical}','total_{exact}')

%%
%-----------Leapfrog Verlet------------------------------------------------
x_leap = zeros(length(t),1);         w = zeros(length(t),1);

x_leap(1) = 1;       w(1) = 1 - 0.5*DeltaT*(-k*x_leap(1))/m;
for i = 1:(length(t)-1)
    w(i+1) = w(i) + DeltaT*(-k*x_leap(i))/m;
    x_leap(i+1) = x_leap(i) + DeltaT*w(i+1);
end
%--------------------------------------------------------------------------

figure()
hold on
plot(t,x_leap,'ko')
plot(t, x_anal(t))

pot = 0.5*k*x_leap.^2;
kin = 0.5*m*w.^2;
tot_numerical = pot + kin;
total_leap_normal = tot_numerical;

figure()
hold on
plot(t,pot);
plot(t,kin);
plot(t,tot_numerical);
plot(t,tot_exact);
legend('potential_{numerical}','kinetic_{numerical}','total_{numerical}','total_{exact}')

v_leap = zeros(length(t),1);
for i = 1:length(t)-1
    v_leap(i) = 0.5 * (w(i+1) + w(i));
end
v_leap(end) = w(end) + 0.5*DeltaT*(-k* (x_leap(end)))/m;

pot = 0.5*k*x_leap.^2;
kin = 0.5*m*v_leap.^2;
tot_numerical = pot + kin;
total_leap_avg = tot_numerical;

figure()
hold on
plot(t,pot);
plot(t,kin);
plot(t,tot_numerical);
plot(t,tot_exact);
legend('potential_{numerical}','kinetic_{numerical}','total_{numerical}','total_{exact}')

figure()
hold on
plot(t,total_leap_normal)
plot(t,total_leap_avg)
plot(t, tot_exact)
%% ODE MATLAB
x0 = [1 1];                                            % initial conditions
tspan = [0; T];
[t_ode, xy] = ode45(@harmOscill, tspan, x0);
figure()
plot(t_ode, xy(:,1));

pot = 0.5*k*xy(:,1).^2;
kin = 0.5*m*xy(:,2).^2;
tot_numerical = pot + kin;
total_ode45 = tot_numerical;

figure()
hold on
plot(t_ode,pot);
plot(t_ode,kin);
plot(t_ode,total_ode45);
%plot(t,tot_exact);
legend('potential_{numerical}','kinetic_{numerical}','total_{numerical}')
%% comparing energies
figure()

hold on
plot(t,total_euler)
plot(t_ode,total_ode45,'k*');
plot(t,total_leap_avg)
plot(t, tot_exact)
legend('Euler','ODE45','Leapfrog','Exact')

%% Friction using ODE
x0 = [1 1];                             % initial conditions
tspan = [0; T];
[t_ode, xy] = ode45(@harmOscill_fric, tspan, x0);
figure()
plot(t_ode, xy(:,1));

pot = 0.5*k*xy(:,1).^2;
kin = 0.5*m*xy(:,2).^2;
tot_numerical = pot + kin;
total_fric_ode45 = tot_numerical;

figure()
hold on
plot(t_ode,pot);
plot(t_ode,kin);
plot(t_ode,total_fric_ode45);
legend('potential_{numerical}','kinetic_{numerical}','total_{numerical}')

%% Friction using Leapfrog
x_leap = zeros(length(t),1);         w = zeros(length(t),1);

x_leap(1) = 1;       
w(1) = (1 - 0.5*DeltaT*(-k*x_leap(1))/m )/(1 - (gamma/m*DeltaT/2));        % leapfrog 1/2 velocity         
for i = 1:(length(t)-1)
    w(i+1) = w(i) + DeltaT*(-k*x_leap(i))/m - gamma*w(i)*DeltaT/m ;
    x_leap(i+1) = x_leap(i) + DeltaT*w(i+1);
end

figure()
hold on
plot(t,x_leap)

v_leap = zeros(length(t),1);
for i = 1:length(t)-1
    v_leap(i) = 0.5 * (w(i+1) + w(i));
end
v_leap(end) = w(end) + 0.5*DeltaT*(-k*(x_leap(end)))/m - gamma*w(i)*0.5*DeltaT/m;

pot = 0.5*k*x_leap.^2;
kin = 0.5*m*v_leap.^2;
tot_numerical = pot + kin;
total_leap_avg = tot_numerical;

figure()
hold on
plot(t,pot);
plot(t,kin);
plot(t,tot_numerical);
plot(t,tot_exact);
legend('potential_{numerical}','kinetic_{numerical}','total_{numerical}','total_{exact}')
%% Force using ODE
x0 = [1 1];                             % initial conditions
tspan = [0; T];
[t_ode, xy] = ode45(@harmOscill_fric_force, tspan, x0);
figure()
plot(t_ode, xy(:,1));
