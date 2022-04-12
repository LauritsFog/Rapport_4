function [threshold, idx] = computeBaselineThreshold(dataC1, dataC2)
    meanDif = abs(mean(dataC1,2)-mean(dataC2,2));

    [maxVal, idx] = max(meanDif);

    meanThresholds = (mean(dataC1,2)+mean(dataC2,2))/2;

    threshold = meanThresholds(idx);
end