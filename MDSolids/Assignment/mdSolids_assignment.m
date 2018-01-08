%% MD solids

clc; clear; close all;

m = 2;                                                                     % Equal Mass of particles
k = 5;                                                                    % Spring Constant
omg_0 = (k/m)^0.5;                                                         % Elastic frequency 
eta = ((log(2)*omg_0/pi)^2)/( 1 + (log(2)/pi)^2);                          % coeffecient of rest = 0.5
eta = sqrt(eta);
re = 0.5;                                                                  % Equillibrium length of spring 
gamma = 0;  %2*m*eta;
omg = (omg_0^2 - eta^2)^0.5;                                               % Eigen frequency 

N = 50;                                                                    % No. of contact cycles
tc = pi/omg;                                                               % Contact duration
Nt = 50;                                                                  % No of steps within contact
deltaT = tc/Nt;                                                            % Time step, DT = Tc/Nt 
t = 0:deltaT:N*tc;                                                         % Time vector
Np = 3;                                                                   % No. of particles in the x and y directions                                               
Vo = 0.1;                                                                  % Supplied velocity 

% Allocating positions and velocities
position.x = zeros(size(t,2),Np^2);                                          
position.y = zeros(size(t,2),Np^2);
velocity.x = zeros(size(t,2),Np^2);
velocity.y = zeros(size(t,2),Np^2);
f.x = zeros(size(t,2),Np^2);
f.y = zeros(size(t,2),Np^2);

% Initialize pot and kin energies
U = zeros(size(t,2),1);
E_k = zeros(size(t,2),1);                                                           

% Initial Conditions 
xp = 0:re:(Np-1)*re;
yp = 0:re:(Np-1)*re;

[Xp, Yp] = meshgrid(xp,yp);

% Initial Positions
for i = 1:Np
    for j = 1:Np
        k_part = (i-1)*Np + j;                   % particle index (for Np = 4, k = 16....) 
        position.x(1,k_part) = Xp(i,j);
        position.y(1,k_part) = Yp(i,j);
    end
end

% Initial Velocities
velocity.y(1,Np^2) = Vo;                         % velocity for the top-right corner particle 

% Connectivity matrix for 2D square lattice
[ C ] = cMatrix_2D_distMethod( Np, position ); 

% Initial Energies
E_k(1) = sum( 0.5*m*(velocity.x(1,:).^2 + velocity.y(1,:).^2));
[U(1),~,~] = forces2D(position.x(1,:),position.y(1,:),k, re, C);

% constraint and no constraint indices
no_constraint = 1:Np^2;                         % vector of particle nos [1 2 3 ....]
bottom_constraint = 1:Np;
left_constraint = 1:Np:(Np*(Np-1))+1;
net_constraint = union(bottom_constraint, left_constraint);
no_constraint = setdiff(no_constraint, net_constraint);

% forces at the first step
[~, f.x(1,:), f.y(1,:)] = forces2D(position.x(1,:),position.y(1,:), k,re, C);

% first step
x_prelim = position.x(1,:) - velocity.x(1,:)*deltaT;
y_prelim = position.y(1,:) - velocity.y(1,:)*deltaT;

% x and y position update for all non-constrained particles
%position.x(2,no_constraint) = 2*position.x(1,no_constraint) - x_prelim(no_constraint) + (deltaT^2)*(f.x(1,no_constraint)/m);
%position.y(2,no_constraint) = 2*position.y(1,no_constraint) - y_prelim(no_constraint) + (deltaT^2)*(f.y(1,no_constraint)/m);

denom = (1 + gamma*deltaT/2/m);

position.x(2,no_constraint) = (2*position.x(1,no_constraint) + (-1 + (gamma*deltaT/2/m))*x_prelim(no_constraint) + (deltaT^2)*(f.x(1,no_constraint)/m))/denom;
position.y(2,no_constraint) = (2*position.y(1,no_constraint) + (-1 + (gamma*deltaT/2/m))*y_prelim(no_constraint) + (deltaT^2)*(f.y(1,no_constraint)/m))/denom;


% y update for left constraint
%position.y(2,left_constraint) = 2*position.y(1,left_constraint) - y_prelim(left_constraint) + (deltaT^2)*(f.y(1,left_constraint)/m);
position.y(2,left_constraint) = (2*position.y(1,left_constraint) + (-1 + (gamma*deltaT/2/m))*y_prelim(left_constraint) + (deltaT^2)*(f.y(1,left_constraint)/m))/denom;


% x update for bottom constraint
%position.x(2,bottom_constraint) = 2*position.x(1,bottom_constraint) - x_prelim(bottom_constraint) + (deltaT^2)*(f.x(1,bottom_constraint)/m);
position.x(2,bottom_constraint) = (2*position.x(1,bottom_constraint) + (-1 + (gamma*deltaT/2/m))*x_prelim(bottom_constraint) + (deltaT^2)*(f.x(1,bottom_constraint)/m))/denom;


% constraints
position.y(2,bottom_constraint) = 0;
position.x(2,left_constraint) = 0;

% dynamics
for i=2:size(t,2)-1
    % forces
    [~,f.x(i,:), f.y(i,:)] = forces2D(position.x(i,:),position.y(i,:),k, re, C );
  
    % -----verlet ------%
    % x and y position update for all non-constrained particles
    %position.x(i+1,no_constraint) = 2*position.x(i,no_constraint) - position.x(i-1,no_constraint) + (deltaT^2)*(f.x(i,no_constraint)/m);
    position.x(i+1,no_constraint) = (2*position.x(i,no_constraint) + (-1 + (gamma*deltaT/2/m))*position.x(i-1,no_constraint) + (deltaT^2)*(f.x(i,no_constraint)/m))/denom;
    %position.y(i+1,no_constraint) = 2*position.y(i,no_constraint) - position.y(i-1,no_constraint) + (deltaT^2)*(f.y(i,no_constraint)/m);
    position.y(i+1,no_constraint) = (2*position.y(i,no_constraint) + (-1 + (gamma*deltaT/2/m))*position.y(i-1,no_constraint) + (deltaT^2)*(f.y(i,no_constraint)/m))/denom;

    
    % y update for left constraint
    %position.y(i+1,left_constraint) = 2*position.y(i,left_constraint) - position.y(i-1,left_constraint) + (deltaT^2)*(f.y(i,left_constraint)/m);
    position.y(i+1,left_constraint) = (2*position.y(i,left_constraint) + (-1 + (gamma*deltaT/2/m))*position.y(i-1,left_constraint) + (deltaT^2)*(f.y(i,left_constraint)/m))/denom;

    
    % x update for bottom constraint
    %position.x(i+1,bottom_constraint) = 2*position.x(i,bottom_constraint) - position.x(i-1,bottom_constraint) + (deltaT^2)*(f.x(i,bottom_constraint)/m);
    position.x(i+1,bottom_constraint) = (2*position.x(i,bottom_constraint) + (-1 + (gamma*deltaT/2/m))*position.x(i-1,bottom_constraint) + (deltaT^2)*(f.x(i,bottom_constraint)/m))/denom;

    
    % constraints
    position.y(i+1,bottom_constraint) = 0;
    position.x(i+1,left_constraint) = 0;
    
    % central difference
    velocity.x(i,:) = (position.x(i+1,:) - position.x(i-1,:))/(2*deltaT);
    velocity.y(i,:) = (position.y(i+1,:) - position.y(i-1,:))/(2*deltaT);
    
    E_k(i) = sum( 0.5*m*(velocity.x(i,:).^2 + velocity.y(i,:).^2));
    [U(i),~,~] = forces2D(position.x(i,:),position.y(i,:),k, re, C);
end
velocity.x(size(t,2), :) = (position.x(size(t,2),:) - position.x(size(t,2)-1,:))/(deltaT);    
velocity.y(size(t,2), :) = (position.y(size(t,2),:) - position.y(size(t,2)-1,:))/(deltaT);    

E_k(size(t,2)) = sum( 0.5*m*(velocity.x(size(t,2),:).^2 + velocity.y(size(t,2),:).^2));
[U(size(t,2)),~,~] = forces2D(position.x(size(t,2),:),position.y(size(t,2),:),k, re, C);

% plot 
figure()
for i = 1:Np^2
    hold on;
    plot(position.x(:,i), t, '.')
end

figure()
hold on
plot(t, U)
plot(t, E_k)
plot(t, U+E_k)
legend('Potential', 'Kinetic', 'Total')

% save data for animation
filename =  'plotting.mat';
save(filename)
