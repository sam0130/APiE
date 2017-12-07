function [ kin, pot, tot ] = calcEnergy(x,v,m,k)
%CALCENERGY calculates the kinetic, potential and total energies
 

kin = 0.5*k*x.^2;
pot = 0.5*m*v.^2;
tot = pot + kin;

end

