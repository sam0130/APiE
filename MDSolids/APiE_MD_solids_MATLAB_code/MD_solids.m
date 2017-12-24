%% MD solids

clc; clear; close all;

global k x_e

m = 2;
k = 5;
gamma = 0;
x_e = 0.5;
omg0 = sqrt(k/m);
omg = sqrt(omg0^2 - gamma^2);

n = 10;          % no. of contact cycles
tc = pi/omg;    % contact duration
Nt = 50;        % no of steps within contact
deltaT = tc/Nt;
T = n*tc;

steps = Nt*n+1;
x = zeros(steps,2);
v = zeros(steps,2);
fc = zeros(steps,2);

x(1,:) = [0 x_e];
v(1,:) = [0 0.1];

n = (x(1,1) - x(1,2))/abs(x(1,1) - x(1,2));
delta = (abs(x(1,1) - x(1,2)) - x_e)*n;
fc(1,1) = -k*delta;
fc(1,2) = -1*fc(1,1);

% first step
x_prelim = x(1,:) - v(1,:)*deltaT;
x(2,:) = 2*x(1,:) - x_prelim + (deltaT^2)*(fc(1,:)/m);

% dynamics
for i=2:steps-1
    
    % force
    fc(i,:) = force(x(i,:));
    
    % verlet
    x(i+1,:) = 2*x(i,:) - x(i-1,:) + (deltaT^2)*(fc(i,:)/m);
    
    % central difference
    v(i,:) = (x(i+1,:) - x(i-1,:))/(2*deltaT);
end

time = 0:deltaT:T;
% plot 
plot(x(:,1), time, '.')
hold on;
plot(x(:,2), time, '.')
legend('particle 1', 'particle 2')