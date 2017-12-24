function [ f_x, f_y ] = forces2D( x,y,k,re, C )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Np = size(x,2);

f_x = zeros(Np,1);
f_y = zeros(Np,1);
for i = 1:Np-1
    for j = i+1 : Np
        if(C(i,j)==1)
            d = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);    % distance b/w particles
            del = -(d - re);                                % spring compression            
            f_ij = k*del;                                   % force magnitude
            
            n_x = (x(i) - x(j))/((x(i) - x(j))^2);          % x unit vector
            n_y = (y(i) - y(j))/((y(i) - y(j))^2);          % y unit vector
            
            f_ij_x = f_ij*n_x;
            f_ij_y = f_ij*n_y;
            
            % force on particle i
            f_x(i) = f_x(i) + f_ij_x;                                       
            f_y(i) = f_y(i) + f_ij_y;
            
            % force on particle j
            f_x(j) = f_x(j) - f_ij_x;                                       
            f_y(j) = f_y(j) - f_ij_y;
        end
    end
end


end

