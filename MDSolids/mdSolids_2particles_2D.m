%% MD solids

clc; clear; close all;

m = 2;                                                                     % Equal Mass of particles
k = 5;                                                                     % Spring Constant
gamma = 0;                                                                 % Friction parameter                      
r_e.x = 0.5;                                                                 % Equillibrium length of spring
r_e.y = 0.5;
re = (r_e.x^2 + r_e.y^2)^0.5;
omg_0 = (k/m)^0.5;                                                         % Elastic frequency 
eta = gamma/2/m;                                                           % Viscous dissipation                                                
omg = (omg_0^2 - eta^2)^0.5;                                               % Eigen frequency 

N = 10;                                                                    % No. of contact cycles
tc = pi/omg;                                                               % Contact duration
Nt = 50;                                                                   % No of steps within contact
deltaT = tc/Nt;                                                            % Time step, DT = Tc/Nt 
t = 0:deltaT:N*tc;                                                         % Time vector
Np = 2;                                                                    % No. of particles

% Allocating positions and velocities
position.x = zeros(size(t,2),Np);                                          
position.y = zeros(size(t,2),Np);
velocity.x = zeros(size(t,2),Np);
velocity.y = zeros(size(t,2),Np);
f.x = zeros(size(t,2),Np);
f.y = zeros(size(t,2),Np);


C = [0 1; 1 0];                                                            % Connectivity matrix 

% Initial Conditions
position.x(1,:) = [0 r_e.x];                        
position.y(1,:) = [0 r_e.y];                        
velocity.x(1,:) = [0 0.1];                        
velocity.y(1,:) = [0 0]; 

% forces at the first step
[f.x(1,:), f.y(1,:)] = forces2D(position.x(1,:),position.y(1,:), k,re, C);


% first step
x_prelim = position.x(1,:) - velocity.x(1,:)*deltaT;
y_prelim = position.y(1,:) - velocity.y(1,:)*deltaT;
position.x(2,:) = 2*position.x(1,:) - x_prelim + (deltaT^2)*(f.x(1,:)/m);
position.y(2,:) = 2*position.y(1,:) - y_prelim + (deltaT^2)*(f.y(1,:)/m);


% dynamics
for i=2:size(t,2)-1
    
    % forces
    [f.x(i,:), f.y(i,:)] = forces2D(position.x(i,:),position.y(i,:),k, re, C );
       
    % verlet
    position.x(i+1,:) = 2*position.x(i,:) - position.x(i-1,:) + (deltaT^2)*(f.x(i,:)/m);
    position.y(i+1,:) = 2*position.y(i,:) - position.y(i-1,:) + (deltaT^2)*(f.y(i,:)/m);
    
    % central difference
    velocity.x(i,:) = (position.x(i+1,:) - position.x(i-1,:))/(2*deltaT);
    velocity.y(i,:) = (position.y(i+1,:) - position.y(i-1,:))/(2*deltaT);
end
velocity.x(size(t,2), :) = (position.x(size(t,2),:) - position.x(size(t,2)-1,:))/(deltaT);    
velocity.y(size(t,2), :) = (position.y(size(t,2),:) - position.y(size(t,2)-1,:))/(deltaT);    



% plot 
plot(position.x(:,1), t, '.')
hold on;
plot(position.x(:,2), t, '.')
legend('particle 1', 'particle 2')

% save data for animation
filename =  'plotting.mat';
save(filename)
