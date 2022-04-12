function [dataC1, dataC2] = extractClassData(data,classification)
    dataC2_temp = cell(1,1);
    dataC1_temp = cell(1,1);
    
    for i = 1:length(classification)
        if classification(i)
            dataC1_temp{end+1} = data(:,i);
        else
            dataC2_temp{end+1} = data(:,i);
        end
    end
    
    dataC1 = cell2mat(dataC1_temp);
    dataC2 = cell2mat(dataC2_temp);
end