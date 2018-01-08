



% p = [1 0 -omg_0^2 0 (log(2)/pi)^2];
% r = roots(p)

p = [pi^2+(log(2))^2    0   -(((log(2))^2)*omg_0^2)];
eta_roots = roots(p)


eta = ((log(2)*omg_0/pi)^2)/( 1 + (log(2)/pi)^2);
eta = sqrt(eta)