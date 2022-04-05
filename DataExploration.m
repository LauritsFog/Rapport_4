addpath(genpath('Data'));

% 1: time sampling points (minutes). 2: Tracer in arterial blood (kBq / ml). 3..7: Tracer in 5
% different ROI (kBq / ml). 

data1 = table2array(readtable("patient1.csv"));

%%

figure
for i = 1:6
    subplot(3,2,i)
    plot(data1(:,1),data1(:,i+1))
end

%%

% Trying to solve with some random constant values and intitial condition.

C0 = ones(2,1)*4;
K = ones(4,1);
tspan = [0 60];
CA = 1;

[tout,Cout] = ode45(@(t,C) FDGModeldF(t,C,K,CA),tspan,C0);

figure
plot(tout,Cout(:,1),tout,Cout(:,2))

%%


