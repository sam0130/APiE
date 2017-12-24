function [fc] = force(x)
%FORCE Summary of this function goes here

global k x_e

n = (x(1) - x(2))/abs(x(1) - x(2));
delta = (abs(x(1) - x(2)) - x_e)*n;

fc(1,1) = -k*delta;
fc(1,2) = -1*fc(1,1);

end

