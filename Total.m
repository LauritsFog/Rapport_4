clear all , clc, close all

%% Inizializing the data

addpath(genpath('data'));

data = cell(10,1);

for i = 1:10
    data{i} = table2array(readtable("patient"+i+".csv"));
end

%% Data for every patient

% Healthy patients 
patient1    = data{1};
patient2    = data{2};
patient3    = data{3};

% Ill patients
patient4    = data{4};
patient5    = data{5};
patient6    = data{6};

% Unknown patients
patient7    = data{7};
patient8    = data{8};
patient9    = data{9};
patient10   = data{10};

%% Plotting activity curves for all patients

% Number of patients
numpatients = 10;
Regions = 5;

% Patients 1 for 5 regions
figure(1) 

for i=3:7
    subplot(1,5,i-2)
    plot(patient2(:,1),patient2(:,i))
end

figure(2) 

for i=3:7
    subplot(1,5,i-2)
    plot(patient4(:,1),patient4(:,i))
end

%% Integral curves for patient 1 

% The integral curves are found given equations 4 and solved using cumtrapz

P=patient1;
Integrals = 4;
Observations = length(P);

A=zeros(5,4); % 5 regions, 4 integrals
B=zeros(Observations,Integrals,Regions);

for r=3:7
            % The integrals involving C_A. 

            int1 = cumtrapz(P(:,1)',P(:,2));
            int2 = cumtrapz(P(:,1)',int1);

            % The integral involving C_T.

            int3 = cumtrapz(P(:,1)',P(:,r));
            int4 = cumtrapz(P(:,1)',int3);
            
            % Saving all the cummulative integrals for each region
            B(:,:,r-2) = [int1,int2,int3,int4];
            
            % Saving the integral values for each region (extra)
            A(r-2,:)=[int1(end),int2(end),int3(end),int4(end)];
end

% This is extra for understanding
% Defining the integral values for all the regions on patient 1
integrals_reg1=A(1,:);
integrals_reg2=A(2,:);
integrals_reg3=A(3,:);
integrals_reg4=A(4,:);
integrals_reg5=A(5,:);

%% Plotting some of the integral curves

% Plot of patient P with integral curve 3 on all regions
% The plots with integral curve 1 and 2 will be the same for all regions
% since it does not depend on CT. Therefore plotting curve 3 to see the
% difference

figure(3) 

for i=1:5
    
    subplot(1,5,i)
    plot(P(:,1),B(:,3,i))
    
end

%% Solving the ODE with least squares
% This is done for one patient on one region

Region=1;

SystemMatrix = B(:,:,Region);
y = patient1(:,Region+2); % since column 3 is region 1

c = (SystemMatrix'*SystemMatrix) \ (SystemMatrix'*y);

% Solving the k-values for patient 1 region 1
k1 = c(1);
k2 = -c(3)-c(2)/c(1);
k4 = c(4)/k2;
k3 = c(2)/k1-k4;

%% Solving rate constants for each patient on each region

% Patient_matrix=[patient1,patient2,patient3,patient4,patient5,patient6,patient7,patient8,patient9,patient10];

numpatients = 10;
numintegral = 4;
numregions = 5;

K=zeros(numregions,numintegral,numpatients);

for p=1:numpatients

    for r=3:7
    
            % The integrals involving C_A. 

            int1 = cumtrapz(data{p}(:,1)',data{p}(:,2));
            int2 = cumtrapz(data{p}(:,1)',int1);

            % The integral involving C_T.

            int3 = cumtrapz(data{p}(:,1)',data{p}(:,r));
            int4 = cumtrapz(data{p}(:,1)',int3);
            
            % Defining the system matrix and rhs.

            A = [int1, int2, int3, int4];

            y = data{p}(:,r);

            % Solving for c and isolating rate constants. 

            b2 = (A'*A)\(A'*y);

            k1_all = b2(1);
            k2_all = -b2(3)-b2(2)/b2(1);
            k4_all = b2(4)/k2_all;
            k3_all = b2(2)/k1_all-k4_all;
            
            K(r-2,:,p)=[k1_all,k2_all,k3_all,k4_all];
            
    end

end

%% Constants for each patient

K_patient1=K(:,:,1);
K_patient2=K(:,:,2);
K_patient3=K(:,:,3);
K_patient4=K(:,:,4);
K_patient5=K(:,:,5);
K_patient6=K(:,:,6);
K_patient7=K(:,:,7);
K_patient8=K(:,:,8);
K_patient9=K(:,:,9);
K_patient10=K(:,:,10);

%% 

figure (4)

col=['g','y','b'];
line=[0.2,0.5,0.7];
plotStyle = {'g','y','b',0.5,0.2,0.8};
alphas = [0.2 0.2 0.2 0.2;
          0.5 0.5 0.5 0.5;
          0.8 0.8 0.8 0.8];

for i=1:5
    
    subplot(2,5,i)
        
    for p=1:3
        con = {'K1'; 'K2'; 'K3'; 'K4'};
        Axis = K(i,:,p)';
        bar(Axis,0.5,col(p));
        grid on;
        xticklabels(con);
        ylabel('The constants value')
        xlabel('The rate constants')
        hFig.WindowState = 'maximized'; 
        str=sprintf('Healthy patient, Region %d',i);
        title(str)
        hold on 
    end
   
    legend('patient1','patient2','patient3') 
end

for i=1:5
    
    subplot(2,5,i+5)
    
    for p=4:6
        con = {'K1'; 'K2'; 'K3'; 'K4'};
        Axis = K(i,:,p)';
        bar(Axis,0.5,col(p-3)); 
        grid on;
        xticklabels(con);
        ylabel('The constants value')
        xlabel('The rate constants')
        hFig.WindowState = 'maximized';
        str=sprintf('Sick patient, Region %d',i);
        title(str)
        hold on 
    end
    
    legend('patient4','patient5','patient6')
end








