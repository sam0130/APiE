clc; clear variables; close all;
%%%%%%%%%%%%%% Constants  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global m k gamma f omg_f
m = 2;                              % mass;   
k = 5;                              % spring   
gamma = 0.3;                        % friction parameter
f = 2;                              % force amplitude

omg = (k/m)^0.5;                    % body frequency

%omg_f = linspace(0,2*omg,100);                            % force frequency

omg_f = omg+1;

T_period = 2*pi/omg;
omg_1 = sqrt( (omg^2) - (( (gamma/m)^2)/4) );   % modified frequency for friction

n = 15;                             % No. of cycles
T = n*T_period;                     % Total time
N = 10000;                          % No. of time steps    
DeltaT = T/N;
t = 0:DeltaT:T;
%%
x_anal0 = 1;
v_anal0 = 1;



[ x_trans, x_steady, x_anal, v_anal, A_total ] = exactSolutionFn( x_anal0, v_anal0,t );

figure()
hold on
plot(t/T_period,x_trans)
plot(t/T_period,x_steady)
plot(t/T_period, x_anal)
legend('transient','steady state','total')

%disp(gamma/2/m/omg);

% for i = 1:length(omg_f)
% [ ~, ~, ~, A_steady(i) ] = exactSolutionFn( x_anal0, v_anal0,t, omg_f(i) );
% end
% 
% figure()
% hold on
% plot(omg_f,A_steady)

