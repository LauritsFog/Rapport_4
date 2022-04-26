% Udregner S-funktionerne ligesom i Salami-opgaven. Den eneste forskel er
% ved sigmaerne. Ved alpha = 1 bruges pooledSigma ikke, ved alpha = 0 er
% det kun pooledSigma der bruges. Så alpha = [0,1]. 

function[Sf_C1, Sf_C2] = computeLDAFunctions(dataC1,dataC2,pC1,pC2,alpha)

    C1Sigma = cov(dataC1);
    C2Sigma = cov(dataC2);

    meanC1 = mean(dataC1,1);
    meanC2 = mean(dataC2,1);

    lenC1 = length(dataC1);
    lenC2 = length(dataC2);

    pooledSigma = (1/((lenC1-1)+(lenC2-1))).*((lenC1-1).*C1Sigma+(lenC2-1).*C2Sigma);

    alphaC1Sigma = alpha*C1Sigma+(1-alpha)*pooledSigma;
    alphaC2Sigma = alpha*C2Sigma+(1-alpha)*pooledSigma;
    
    alphaC1Sigmainv = pinv(alphaC1Sigma);
    alphaC2Sigmainv = pinv(alphaC2Sigma);

    Sf_C1 = @(x) x'*alphaC1Sigmainv*meanC1'-(1/2)*meanC1*alphaC1Sigmainv*meanC1'+log(pC1);

    Sf_C2 = @(x) x'*alphaC2Sigmainv*meanC2'-(1/2)*meanC2*alphaC2Sigmainv*meanC2'+log(pC2);
end