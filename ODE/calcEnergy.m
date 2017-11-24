function [ kin, pot, tot ] = calcEnergy(x,v)
%CALCENERGY calculates the kinetic, potential and total energies

global m k 

kin = 0.5*k*x.^2;
pot = 0.5*m*v.^2;
tot = pot + kin;

end

