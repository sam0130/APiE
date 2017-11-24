clc; clear; close all;

% dydx = @(x,y) cos(x);
% x = [0 10];
% y0 = 0;
% [T,Yt] = ode45(dydx, x, y0);
% figure()
% plot(T,Yt);
% 
% %%
% ic = 0;
% f = @(t,y) t*y^2 + y;
% [t,y] = ode45(f, [0 0.5], ic);
% figure()
% plot(t,y);
% %%
% f = @(x,y) x*y^2 + y;
% [x,y] = ode23(f,[0 1], 1);
% figure()
% plot(x,y)
% %%
% function main
% [T, Y] = ode15s(@f, [0 3000], [2,0]);
% plot(T,Y(:,1),'-o');
% 
% end
%% lorenz
x0 = [-8 8 27];
tspan = [0,20];
[t, xyz] = ode45(@lorenz, tspan, x0);
plot3(xyz(:,1),xyz(:,2),xyz(:,3))