function [ x_trans, x_steady, x_anal, v_anal, A_steady, A_total ] = ...
               exactSolutionFn( x_anal0, v_anal0,t, m, k, gamma, f, omg_f)
%exactSolutionFn Exact solution of the harmonic oscillator
%   x_trans, v_trans: homogeneous solution
%   x_steady, v_steady: particular solution
%   x_anal, v_anal: total solution
%   A_total:  amplitude

omg = (k/m)^0.5;                    % frequency

omg_prime = sqrt( omg^2 - (gamma/2/m)^2 );                                 % damping reduced frequency
A_steady = (f/m)/(sqrt( (omg^2 - omg_f^2)^2  + 4*(gamma/2/m)^2*omg_f^2 ));        % steady state amplitude

phi_d = 0;                                                                 % force phase
phi = atan( (gamma*omg_f)/(k - m*omg_f^2)) - phi_d;                        % steady solution phase
phi_h = atan((omg_prime* (x_anal0 - A_steady*cos(phi)))/(v_anal0 ...              % transient solution phase
        + (gamma/2/m) * (x_anal0 - A_steady*cos(phi)) - A_steady*omg_f*sin(phi)) );     

A_h = (x_anal0 - A_steady*cos(phi))/sin(phi_h);                                   % transient amplitude  

x_trans = A_h * exp(-gamma/2/m*t) .* sin(omg_prime*t + phi_h);             % transient solution
x_steady = A_steady*cos(omg_f*t - phi);                                           % steady state solution

x_anal = x_trans + x_steady;                                               % net solution 
v_anal = A_h*((-gamma/2/m)*exp(-gamma*t/2/m).*sin(omg_prime*t + phi_h)...
              + omg_prime*cos(omg_prime*t + phi_h).*exp(-gamma*t/2/m)) ...
              - A_steady*omg_f*sin(omg_f*t - phi);

A_total = sqrt((A_h * exp(-gamma/2/m*t)).^2 + A_steady^2);
end

