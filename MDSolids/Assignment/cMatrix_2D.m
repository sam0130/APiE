function [ C ] = cMatrix_2D( Np )
%CMATRIX_2D connectivity matrix for 2D square lattice
%   Considers diagonal and vertical/horizontal springs
% Diagonal

% r = Np; c = Np;                        %# Get the matrix size
% diagVec1 = repmat([ones(c-1,1); 0],r,1);  %# Make the first diagonal vector
%                                           %#   (for horizontal connections)
% diagVec1 = diagVec1(1:end-1);             %# Remove the last value
% diagVec2 = [0; diagVec1(1:(c*(r-1)))];    %# Make the second diagonal vector
%                                           %#   (for anti-diagonal connections)
% diagVec3 = ones(c*(r-1),1);               %# Make the third diagonal vector
%                                           %#   (for vertical connections)
% diagVec4 = diagVec2(2:end-1);             %# Make the fourth diagonal vector
%                                           %#   (for diagonal connections)
% adj = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
%       diag(diagVec2,c-1)+...
%       diag(diagVec3,c)+...
%       diag(diagVec4,c+1);
% C = adj+adj.'; 

% Non-diagonal
r = Np; c = Np;                        %# Get the matrix size
diagVec1 = repmat([ones(c-1,1); 0],r,1);  %# Make the first diagonal vector
                                          %#   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);             %# Remove the last value
diagVec2 = ones(c*(r-1),1);               %# Make the second diagonal vector
                                          %#   (for vertical connections)
adj = diag(diagVec1,1)+...                %# Add the diagonals to a zero matrix
      diag(diagVec2,c);
C = adj+adj.';  


end

