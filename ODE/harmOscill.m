function xv = harmOscill( t,x)
% Define the ODEs to solve through ode45
global k m gamma f omega

xv = [x(2);
      -k/m * x(1) - gamma/m * x(2) + f/m*cos(omega*t)];


end

