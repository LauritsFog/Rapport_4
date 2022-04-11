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

% Used for x-axis values. 

xstr = cell(20,1);
iter = 1;
for r = 1:5
    for k = 1:4
        xstr{iter} = string("r" + r + ", k" + k);
        iter = iter + 1;
    end
end
xaxis = categorical(cellstr(xstr));

%%

c = cool(P);
str = cell(P,1);

figure
for p = 1:P
    plot(xaxis,K(:,p),'Color',c(p,:))
    hold on
    str{p} = string("Patient " + p);
end
legend(str)
ylabel('Rate constants values')
xlabel('Region #, rate constant #')

%%

% Plotting without patient numbers. 

figure
b = bar(xaxis,K);
for i = 1:3
    b(i).FaceColor = 'b';
    b(i+3).FaceColor = 'r';
end
legend([b(1),b(4)],{'Healthy patients','Sick patients'})
ylabel('Rate constants values')
xlabel('Region #, rate constant #')

%%

% Plotting with patient numbers. 

c = cool(8);
str = cell(P,1);

figure
b = bar(xaxis,K);
for i = 1:3
    ii = i+3;
    b(i).FaceColor = c(i,:);
    b(ii).FaceColor = c(i+5,:);
    str{i} = string("Patient " + i);
    str{ii} = string("Patient " + ii);
end
legend(str)
ylabel('Rate constants values')
xlabel('Region #, rate constant #')

