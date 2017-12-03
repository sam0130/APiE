function [ x_trans, x_steady, x_anal, v_anal, A_total ] = exactSolutionFn( x_anal0, v_anal0,t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global m k gamma f omg_f

omg = (k/m)^0.5;                    % frequency

omg_prime = sqrt( omg^2 - (gamma/2/m)^2 );                                 % damping reduced frequency
A = (f/m)/(sqrt( (omg^2 - omg_f^2)^2  + 4*(gamma/2/m)^2*omg_f^2 ));        % steady state amplitude

phi_d = 0;                                                                 % force phase
phi = atan( (gamma*omg_f)/(k - m*omg_f^2)) - phi_d;                        % steady solution phase
phi_h = atan((omg_prime* (x_anal0 - A*cos(phi)))/(v_anal0 ...              % transient solution phase
        + (gamma/2/m) * (x_anal0 - A*cos(phi)) - A*omg_f*sin(phi)) );     

A_h = (x_anal0 - A*cos(phi))/sin(phi_h);                                   % transient amplitude  

x_trans = A_h * exp(-gamma/2/m*t) .* sin(omg_prime*t + phi_h);             % transient solution
x_steady = A*cos(omg_f*t - phi);                                           % steady state solution

x_anal = x_trans + x_steady;                                               % net solution 
v_anal = A_h*((-gamma/2/m)*exp(-gamma*t/2/m).*sin(omg_prime*t + phi_h)...
              + omg_prime*cos(omg_prime*t + phi_h).*exp(-gamma*t/2/m)) ...
              - A*omg_f*sin(omg_f*t - phi);

A_total = sqrt((A_h * exp(-gamma/2/m*t)).^2 + A^2);


%A_steady = (f/m)/( (omg^2 - omg_f^2));

end

