clear variables; close all; clc;

Np = 4;
%% Square lattice
C_sq = zeros(Np,Np);
for i = 1:Np-1
    for j = i+1 : Np
        C_sq(i,j) = 1;
        C_sq(j,i) = 1;
    end
end

%% Linear Lattice
C_lin = zeros(Np,Np);
for i = 1:Np-1
    for j = i+1 : Np
        if(abs(i-j)==1) 
        C_lin(i,j) = 1;
        C_lin(j,i) = 1;
        end
    end
end
