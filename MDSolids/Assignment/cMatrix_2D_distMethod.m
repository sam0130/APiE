function [ C ] = cMatrix_2D_distMethod( Np, position )

Px = position.x(1,:);
Py = position.y(1,:);
C = zeros(Np^2);
for i = 1:((Np^2)-1)
    for j = i+1:Np^2
        
        d = sqrt((Px(i) - Px(j))^2 + (Py(i) - Py(j))^2);    % distance b/w particles
        if(d ==0.5  || d == sqrt(2)*0.5)
             C(i,j) = 1;
         end
    end
end
C = C + C';
end

