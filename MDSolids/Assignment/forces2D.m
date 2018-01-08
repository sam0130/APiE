function [U, f_x, f_y ] = forces2D( x,y,k,re, C )
%FORCES2D Calculates forces and potential energy
%   Each particle pair considered only once
%   Requires the C connectivity matric
%   Uses Newton's third law to calculate equal opposite force on the jth
%   paricle
Np = size(x,2);

f_x = zeros(Np,1);
f_y = zeros(Np,1);
U = 0;
for i = 1:Np-1
    for j = i+1 : Np
        if(C(i,j)~=0)
            d = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);           % distance b/w particles
            del = -(d - re*C(i,j));                                % spring compression            
            f_ij = C(i,j)*k*del;                                   % force magnitude
            n_x = (x(i) - x(j))/d;          % x unit vector
            n_y = (y(i) - y(j))/d;          % y unit vector
            
            f_ij_x = f_ij*n_x;
            f_ij_y = f_ij*n_y;
            
            % force on particle i
            f_x(i) = f_x(i) + f_ij_x;                                       
            f_y(i) = f_y(i) + f_ij_y;
            
            % force on particle j
            f_x(j) = f_x(j) - f_ij_x;                                       
            f_y(j) = f_y(j) - f_ij_y;
            
            % potential energy
            U = U + 0.5*C(i,j)*k*del^2;
            
        end
    end
end


end

