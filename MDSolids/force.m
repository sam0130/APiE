function [fc] = force(x, k, x_e)
%FORCE Summary of this function goes here

n = [1 -1];
fc = k*((x(2) - x(1))-x_e)*n;

end
