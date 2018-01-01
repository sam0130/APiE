% Initial Conditions 
clear variables;
Np = 3;
re = 0.5;
xp = 0:re:(Np-1)*re;
yp = 0:re:(Np-1)*re;
D = zeros(Np^2);
[Xp, Yp] = meshgrid(xp,yp);

for i = 1:Np
    for j = 1:Np
        k_part = (i-1)*Np + j;                   % particle index (for Np = 4, k = 16....) 
        Px(k_part) = Xp(i,j);
        Py(k_part) = Yp(i,j);
    end
end

for i = 1:((Np^2)-1)
    for j = i+1:Np^2
        
        d = sqrt((Px(i) - Px(j))^2 + (Py(i) - Py(j))^2);    % distance b/w particles
        if(d ==0.5 )
             D(i,j) = 1;
         end
    end
end


C = D+D';
