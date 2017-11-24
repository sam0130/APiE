function xv = harmOscill_fric( t,x)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

k = 5; m = 2; gamma = 0.3;

xv = [x(2);
      -k/m * x(1) - gamma/m * x(2)];


end

