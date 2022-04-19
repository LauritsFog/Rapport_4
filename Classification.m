addpath(genpath('Data'));

% 1: time sampling points (minutes). 2: Tracer in arterial blood (kBq / ml). 3..7: Tracer in 5
% different ROI (kBq / ml). 

data = cell(10,1);

for i = 1:10
    data{i} = table2array(readtable("patient"+i+".csv"));
end

%%

N = 6;

K = computeRateConstants(data(1:N));

Ks = (K - mean(K,2))./std(K,[],2);

Kchoice = Ks;

true = [1,1,1,0,0,0];

%%

% Testing LDA. 

pHealthy = 0.5;
pSick = 0.5;
alpha = 0.5;

[Sf_healthy,Sf_sick] = computeLDAFunctions(Kchoice(:,1:3)',Kchoice(:,4:6)',pHealthy,pSick,alpha);

%%

% Testing SVM.

SVM = fitcsvm(Kchoice',true);

CVSVMModel = crossval(SVM);

classLoss = kfoldLoss(CVSVMModel);

%%

% Testing Baseline.

threshold = (mean(mean(Kchoice(:,1:3)))+mean(mean(Kchoice(:,4:6))))/2;
% [threshold, idx] = computeBaselineThreshold(Healthy_train,Sick_train);

%%

alphas = [0,10.^(-3:0)];
alpha = 0;
outlier = linspace(0,0.5,5);

pHealthy = 0.5;
pSick = 0.5;

test_true = zeros(1,6);

CV = cvpartition(6,'Leaveout');

LDAClass = zeros(1,6);
SVMClass = zeros(1,6);
BaselineClass = zeros(1,6);

test_true2 = zeros(length(alphas),5);

for i = 1:CV.NumTestSets
    [Healthy_train, Sick_train] = extractClassData(Kchoice(:,CV.training(i)),true(CV.training(i))); % Extracting training data. 
    
    true_train = true(CV.training(i)); % Extracting true classifications of test observations. 
    
    CV2 = cvpartition(5,'Leaveout');
    
    % Classifications (used for computing Egen).
    
    LDAClass2 = zeros(length(alphas),5);
    SVMClass2 = zeros(length(alphas),5);
    
    for j = 1:CV2.NumTestSets
        InnerData = [Healthy_train, Sick_train]; % Combining training data for inner loop. 
        
        % Extracting healthy data and sick data based on inner training
        % partioning. 
        [Healthy_train2, Sick_train2] = extractClassData(InnerData(:,CV2.training(j)),true_train(CV2.training(j)));
        
        % Extracting test data. 
        K_test2 = InnerData(:,CV2.test(j));
                
        % Getting true class of testing data. 
        test_true2(:,j) = true_train(CV2.test(j));
        
        for s = 1:length(alphas) % Testing every model on testing data j. 
            [Sf_healthy2,Sf_sick2] = computeLDAFunctions(Healthy_train2',Sick_train2',pHealthy,pSick,alphas(s));

            SVM = fitcsvm([Healthy_train2, Sick_train2]',true_train(CV2.training(j)),'OutlierFraction',outlier(s));
        
            % Classifying with LDA. 
            if Sf_healthy2(K_test2) > Sf_sick2(K_test2)
                LDAClass2(s,j) = 1;
            end

            % Classifying with SVM. 
            SVMClass2(s,j) = predict(SVM,K_test2');
        end
    end
    
    % Computing generalization errors. 
    SVMEgen = sum(1/5*abs(SVMClass2-test_true2),2);
    LDAEgen = sum(1/5*abs(LDAClass2-test_true2),2);
    
    % Finding index of parameter with lowest Egen. 
    [minval,SVMidx] = min(SVMEgen);
    [minval,LDAidx] = min(LDAEgen);
    
    K_test = Kchoice(:,CV.test(i)); % Extracting testing data. 
    
    test_true(i) = true(CV.test(i)); % Extracting true classifications of test data. 
    
    % Training LDA. 
    [Sf_healthy,Sf_sick] = computeLDAFunctions(Healthy_train',Sick_train',pHealthy,pSick,alphas(LDAidx));
    
    % Training SVM. 
    SVM = fitcsvm([Healthy_train, Sick_train]',true(CV.training(i)),'OutlierFraction',outlier(SVMidx));
    
    % Training Baseline. 
    threshold = (mean(mean(Healthy_train))+mean(mean(Sick_train)))/2;
    
    % Classifying with LDA. 
    if Sf_healthy(K_test) > Sf_sick(K_test)
        LDAClass(i) = 1;
    end
    
    % Classifying with SVM. 
    SVMClass(i) = predict(SVM,K_test');
    
    % Classifying with Baseline. 
    if mean(K_test) > threshold
        BaselineClass(i) = 1;
    end
end

% Computing generalization errors. 

genErrLDA = nnz(LDAClass-test_true)/N;
genErrSVM = nnz(SVMClass-test_true)/N;
genErrBaseline = nnz(BaselineClass-test_true)/N;

%%

% Plotting confusion matrices. 

figure
subplot(3,1,1)
confusionchart(LDAClass,test_true);
title('LDA')
subplot(3,1,2)
confusionchart(SVMClass,test_true);
title('SVM')
subplot(3,1,3)
confusionchart(BaselineClass,test_true);
title('Baseline')

%%

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

errorindex = 1;

Kplot = [Kchoice(:,CV.training(errorindex)), Kchoice(:,CV.test(errorindex))];

trueplot = true(CV.training(errorindex));

figure
b = bar(xaxis,Kplot);
for i = 1:5
    if true(i) == 1
        b(i).FaceColor = 'b';
    else
        b(i).FaceColor = 'r';
    end
end
legend([b(1),b(5),b(6)],{'Healthy training data','Sick training data','Testing data'})
ylabel('Rate constants values')
xlabel('Region #, rate constant #')
