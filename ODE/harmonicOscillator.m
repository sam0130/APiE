% Harmonic Oscillator
clc; close all; clear;
T = 4; %
n = 100;
A = 1;      % Amplitude
k = 5;      % Spring constant
m = 2;      % mass
omg = sqrt(k/m);
DeltaT = T/n;
t = linspace(0,T,n);
v(1) = A*omg;
x(1) = 0;
x_anal = A*sin(omg*t);
for i = 1:(length(t) - 1)
    x(i+1) = x(i) + v(i)*DeltaT;
    v(i+1) = v(i) + (-k*x(i))*DeltaT/m;
end

%% 
% PostProcessing
figure()
hold on
plot(t(1:20:end),x(1:20:end),'ko');
plot(t,x_anal);
%err(:,1) = zeros(length(t),1);
%err(:,2) = x_anal - x;
%errorbar(t(1:20:end),x(1:20:end),err(1:20:end));

pot  = 0.5 * k * x.^2;
kin  = 0.5 * m * v.^2;
figure()
hold on
plot(t,pot);
plot(t,kin);
plot(t,kin+pot);
plot(t,0.5*m*v(1)^2*ones(length(t),1)')
legend('Potential','Kinetic','Total','Initial')





