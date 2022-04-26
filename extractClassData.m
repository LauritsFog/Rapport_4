% Tager data med tilhørende klassifikation som input. Giver dataen delt i
% to som output, hvor den ene del tilhører klasse 1, og den anden klasse 2.
% Classification inputtet består af 1'er og 0'er. 

function [dataC1, dataC2] = extractClassData(data,classification)
    dataC2_temp = cell(1,1);
    dataC1_temp = cell(1,1);
    
    for i = 1:length(classification)
        if classification(i) % If class = 1 the column is put into dataC1. 
            dataC1_temp{end+1} = data(:,i);
        else % If class = 0 the column is put into dataC2. 
            dataC2_temp{end+1} = data(:,i);
        end
    end
    
    % Cells are converted into matrices. 
    
    dataC1 = cell2mat(dataC1_temp);
    dataC2 = cell2mat(dataC2_temp);
end