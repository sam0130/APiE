%% MD solids

clc; clear; close all;

m = 2;                                                                     % Equal Mass of particles
k = 5;                                                                     % Spring Constant
gamma = 0;                                                                 % Friction parameter                      
re = 0.5;                                                                  % Equillibrium length of spring 
omg_0 = (k/m)^0.5;                                                         % Elastic frequency 
eta = gamma/2/m;                                                           % Viscous dissipation                                                
omg = (omg_0^2 - eta^2)^0.5;                                               % Eigen frequency 

N = 100;                                                                    % No. of contact cycles
tc = pi/omg;                                                               % Contact duration
Nt = 100;                                                                   % No of steps within contact
deltaT = tc/Nt;                                                            % Time step, DT = Tc/Nt 
t = 0:deltaT:N*tc;                                                         % Time vector
Np = 5;                                                                    % No. of particles

% Allocating positions and velocities
position.x = zeros(size(t,2),Np);                                          
position.y = zeros(size(t,2),Np);
velocity.x = zeros(size(t,2),Np);
velocity.y = zeros(size(t,2),Np);
f.x = zeros(size(t,2),Np);
f.y = zeros(size(t,2),Np);

% Initialize pot and kin energies
U = zeros(size(t,2),1);
E_k = zeros(size(t,2),1);

% Connectivity matrix for 1D linear system
C = zeros(Np,Np);
for i = 1:Np-1
    for j = i+1 : Np
        if(abs(i-j)==1) 
        C(i,j) = 1;
        C(j,i) = 1;
        end
    end
end                                                            

% Initial Conditions
position.x(1,:) = 0:re:(Np-1)*re;                        
position.y(1,:) = [0 re re 0.5*re 0];                       
velocity.x(1,(Np+1)/2) = 0.03;                        
velocity.y(1,(Np+1)/2) = 0.03; 

E_k(1) = sum( 0.5*m*(velocity.x(1,:).^2 + velocity.y(1,:).^2));
[U(1),~,~] = forces2D(position.x(1,:),position.y(1,:),k, re, C);

% forces at the first step
[~, f.x(1,:), f.y(1,:)] = forces2D(position.x(1,:),position.y(1,:), k,re, C);

% first step
x_prelim = position.x(1,:) - velocity.x(1,:)*deltaT;
y_prelim = position.y(1,:) - velocity.y(1,:)*deltaT;

% Left wall                         
position.x(2,1) = 0;           
position.y(2,1) = 0;            

% Right wall
position.x(2,end) = (Np-1)*re;
position.y(2,end) = 0; 

position.x(2,2:end-1) = 2*position.x(1,2:end-1) - x_prelim(2:end-1) + (deltaT^2)*(f.x(1,2:end-1)/m);
position.y(2,2:end-1) = 2*position.y(1,2:end-1) - y_prelim(2:end-1) + (deltaT^2)*(f.y(1,2:end-1)/m);

% dynamics
for i=2:size(t,2)-1
    % forces
    [~,f.x(i,:), f.y(i,:)] = forces2D(position.x(i,:),position.y(i,:),k, re, C );
  
    % verlet
    
    position.x(i+1,1) = 0;
    position.y(i+1,1) = 0; 
    
    position.x(i+1,2:end-1) = 2*position.x(i,2:end-1) - position.x(i-1,2:end-1) + (deltaT^2)*(f.x(i,2:end-1)/m);
    position.y(i+1,2:end-1) = 2*position.y(i,2:end-1) - position.y(i-1,2:end-1) + (deltaT^2)*(f.y(i,2:end-1)/m);
    
    position.x(i+1,end) = (Np-1)*re;
    position.y(i+1,end) = 0; 
    
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
for i = 1:Np
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
