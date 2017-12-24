%% MD solids

clc; clear; close all;


m = 2;                                                                     % Equal Mass of particles
k = 5;                                                                     % Spring Constant
gamma = 0;                                                                 % Friction parameter                      
x_e = 0.5;                                                                 % Equillibrium length of spring
omg_0 = (k/m)^0.5;                                                         % Elastic frequency 
eta = gamma/2/m;                                                           % Viscous dissipation                                                
omg = (omg_0^2 - eta^2)^0.5;                                               % Eigen frequency 

N = 10;                                                                    % No. of contact cycles
tc = pi/omg;                                                               % Contact duration
Nt = 50;                                                                   % No of steps within contact
deltaT = tc/Nt;                                                            % Time step, DT = Tc/Nt 
t = 0:deltaT:N*tc;                                                         % Time vector

x = zeros(size(t,2),2);
v = zeros(size(t,2),2);
fc = zeros(size(t,2),2);


x(1,:) = [0 x_e];                                                          % Initial positions 
v(1,:) = [0 0.1];                                                          % Initial velocities
fc(1,:) = force_2particles_1D(x(1,:),k,x_e);


% first step
x_prelim = x(1,:) - v(1,:)*deltaT;
x(2,:) = 2*x(1,:) - x_prelim + (deltaT^2)*(fc(1,:)/m);

% dynamics
for i=2:size(t,2)-1
    
    % force
    fc(i,:) = force_2particles_1D(x(i,:),k, x_e);
    
    % verlet
    x(i+1,:) = 2*x(i,:) - x(i-1,:) + (deltaT^2)*(fc(i,:)/m);
    
    % central difference
    v(i,:) = (x(i+1,:) - x(i-1,:))/(2*deltaT);
end
v(size(t,2), :) = (x(size(t,2),:) - x(size(t,2)-1,:))/(deltaT);    


% plot 
plot(x(:,1), t, '.')
hold on;
plot(x(:,2), t, '.')
legend('particle 1', 'particle 2')