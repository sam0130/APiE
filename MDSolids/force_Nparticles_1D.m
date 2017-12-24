function [fc] = force_Nparticles_1D(x, k, x_e)
%FORCE Summary of this function goes here
fc = zeros(size(x,2),1);
for i = 1:size(x,2)
    if(i==1)
        fc(i) = k*((x(i+1) - x(i))-x_e);
    elseif(i==size(x,2))
        fc(i) = -k*((x(i) - x(i-1))-x_e);
    else
        fc(i) = k*(x(i+1) - 2*x(i) + x(i-1));
    end
end
