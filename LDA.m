addpath(genpath('Data'));

% 1: time sampling points (minutes). 2: Tracer in arterial blood (kBq / ml). 3..7: Tracer in 5
% different ROI (kBq / ml). 

data = cell(10,1);

for i = 1:10
    data{i} = table2array(readtable("patient"+i+".csv"));
end

%%

P = 6;

K = ComputeRateConstants(data(1:P));

%%

pHealthy = 0.5;
pSick = 0.5;

[Sf_healthy,Sf_sick] = computeSFunctions(K(:,1:3)',K(:,4:6)',pHealthy,pSick);

%%


