function xv = harmOscill( t,x)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

k = 5; m = 2;

xv = [x(2);
      -k/m * x(1)];


end

