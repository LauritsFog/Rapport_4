addpath(genpath('Data'));

% 1: time sampling points (minutes). 2: Tracer in arterial blood (kBq / ml). 3..7: Tracer in 5
% different ROI (kBq / ml). 

data = cell(10,1);

for i = 1:10
    data{i} = table2array(readtable("patient"+i+".csv"));
end

%%

P = 6;

K = computeRateConstants(data(1:P));

Ks = (K - mean(K,2))./std(K,[],2);

Kchoice = Ks;

true = [1,1,1,0,0,0];

%%

% Testing LDA. 

pHealthy = 0.5;
pSick = 0.5;

[Sf_healthy,Sf_sick] = computeLDAFunctions(Kchoice(:,1:3)',Kchoice(:,4:6)',pHealthy,pSick);

%%

% Testing SVM.

SVM = fitcsvm(Kchoice',true);

%%

% Testing Baseline.

threshold = (mean(mean(Kchoice(:,1:3)))+mean(mean(Kchoice(:,4:6))))/2;
% [threshold, idx] = computeBaselineThreshold(Healthy_train,Sick_train);

%%

pHealthy = 0.5;
pSick = 0.5;

test_true = zeros(1,6);

CV = cvpartition(6,'Leaveout');

LDAClass = zeros(1,6);
SVMClass = zeros(1,6);
BaselineClass = zeros(1,6);

for i = 1:CV.NumTestSets
    
    CV_train = CV.training(i);
    
    ExtractHealthy = logical([CV_train(1:3);0;0;0]);
    ExtractSick = logical([0;0;0;CV_train(4:6)]);
    
    Healthy_train = Kchoice(:,ExtractHealthy); % Extracting healthy training data. 
    Sick_train = Kchoice(:,ExtractSick); % Extracting sick training data. 
    
    K_test = Kchoice(:,CV.test(i)); % Extracting testing data. 
    
    test_true(:,i) = true(CV.test(i)); % Extracting true classifications of test data. 
    
    % Training LDA. 
    [Sf_healthy,Sf_sick] = computeLDAFunctions(Healthy_train',Sick_train',pHealthy,pSick);
    
    % Training SVM. 
    SVM = fitcsvm([Healthy_train, Sick_train]',true(CV_train));
    
    % Training Baseline. 
    threshold = (mean(mean(Healthy_train))+mean(mean(Sick_train)))/2;
    
    % Classifying with LDA. 
    if Sf_healthy(K_test)/Sf_sick(K_test) > 1
        LDAClass(i) = 1;
    end
    % Classifying with SVM. 
    SVMClass(i) = predict(SVM,K_test');
    % Classifyiong with Baseline. 
    if mean(K_test) > threshold
        BaselineClass(i) = 1;
    end
end
