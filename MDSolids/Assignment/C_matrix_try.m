clear variables; close all; clc;

Np = 9;
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
%%
N = 4;  
r = N; c = N;                        %# Get the matrix size
diagVec1 = repmat([ones(c-1,1); 0],r,1);  %# Make the first diagonal vector
                                          %#   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);             %# Remove the last value
diagVec2 = [0; diagVec1(1:(c*(r-1)))];    %# Make the second diagonal vector
                                          %#   (for anti-diagonal connections)
diagVec3 = ones(c*(r-1),1);               %# Make the third diagonal vector
                                          %#   (for vertical connections)
diagVec4 = diagVec2(2:end-1);             %# Make the fourth diagonal vector
                                          %#   (for diagonal connections)
adj = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
      diag(diagVec2,c-1)+...
      diag(diagVec3,c)+...
      diag(diagVec4,c+1);
adj = adj+adj.'; 

spy(adj)
