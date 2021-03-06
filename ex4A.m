addpath(genpath('Data'));

% 1: time sampling points (minutes). 2: Tracer in arterial blood (kBq / ml). 3..7: Tracer in 5
% different ROI (kBq / ml). 

data = cell(10,1);

for i = 1:10
    data{i} = table2array(readtable("patient"+i+".csv"));
end
%%
% Using composite trapezoid rule to integrate data for well patients. 
int_well = nan(158,5,3);

for p = 1:3
for i = 3:7
int_well(:,i-2,p) = cumtrapz(data{p}(:,1)',data{p}(:,i));
end
end

%%
% Using composite trapezoid rule to integrate data for ill patients.
int_ill = nan(158,5,3);

for p = 4:6
for i = 3:7
int_ill(:,i-2,p-3) = cumtrapz(data{p}(:,1)',data{p}(:,i));
end
end

%%
% takes mean of each region for well and ill patients an puts into vector
well_mean = [mean(int_well(end,1,1:3)),mean(int_well(end,2,1:3)),mean(int_well(end,3,1:3)),mean(int_well(end,4,1:3)),mean(int_well(end,5,1:3))];
ill_mean = [mean(int_ill(end,1,1:3)),mean(int_ill(end,2,1:3)),mean(int_ill(end,3,1:3)),mean(int_ill(end,4,1:3)),mean(int_ill(end,5,1:3))];

% Making some kind of baseline metod where making a threshold with the mean
% of all the regions
(mean(well_mean)+mean(ill_mean))/2
% a person that is ill has in general lower integrals than a person that is
% well
%so everything under 1.1776e+03 is ill and everything over this is well

%%
% trying on one of the patients we do not know 
for i = 3:7
int_p1(:,i-2) = cumtrapz(data{7}(:,1)',data{7}(:,i));
end
mean_p1 = mean(int_p1(end,1:5))
%This patient is ill (close to well)

% trying on one of the patients we do not know 
for i = 3:7
int_p2(:,i-2) = cumtrapz(data{8}(:,1)',data{8}(:,i));
end
mean_p2 = mean(int_p2(end,1:5))
%This patient is ill

% trying on one of the patients we do not know 
for i = 3:7
int_p3(:,i-2) = cumtrapz(data{9}(:,1)',data{9}(:,i));
end
mean_p3 = mean(int_p3(end,1:5))
%This patient is ill

% trying on one of the patients we do not know 
for i = 3:7
int_p4(:,i-2) = cumtrapz(data{10}(:,1)',data{10}(:,i));
end
mean_p4 = mean(int_p4(end,1:5))
%This patient is well