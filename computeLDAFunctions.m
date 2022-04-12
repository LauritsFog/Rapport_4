function[Sf_C1, Sf_C2] = computeLDAFunctions(dataC1,dataC2,pC1,pC2)

    C1Sigma = cov(dataC1);
    C2Sigma = cov(dataC2);

    meanC1 = mean(dataC1,1);
    meanC2 = mean(dataC2,1);

    lenC1 = length(dataC1);
    lenC2 = length(dataC2);

    pooledSigma = (1/((lenC1-1)+(lenC2-1))).*((lenC1-1).*C1Sigma+(lenC2-1).*C2Sigma);

    pooledSigmaInv = pinv(pooledSigma);

    Sf_C1 = @(x) x'*pooledSigmaInv*meanC1'-(1/2)*meanC1*pooledSigmaInv*meanC1'+log(pC1);

    Sf_C2 = @(x) x'*pooledSigmaInv*meanC2'-(1/2)*meanC2*pooledSigmaInv*meanC2'+log(pC2);
end