function [ U] = pot_energy( x, x_e, k )
%POT_ENERGY Calculates the potential energy of the system per time step
%   Supply the position vector at the current time step
%   U = 0.5*k*((x(i+1) - x(i) - x_e).^2

Np = size(x,2);
U = 0;
for i = 1:Np-1
    U = U + 0.5*k*((x(i+1) - x(i)) - x_e)^2;
end

