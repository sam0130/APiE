function xyz = lorenz(t,x);
%LORENZ: Defines the Lorenz equations.
sig=10;
beta=8/3;
rho=28;
xyz=[-sig*x(1) + sig*x(2); rho*x(1) - x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];