function [K] = computeRateConstants(data)
    P = length(data); % From patient 1 to P.

    K = cell(5,P);

    for p = 1:P % For each patient
        for r = 3:7 % For each region

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

            c = (A'*A)\(A'*y);

            k1 = c(1);
            k2 = -c(3)-c(2)/c(1);
            k4 = c(4)/k2;
            k3 = c(2)/k1-k4;

            K_temp = [k1;k2;k3;k4];

            K{r-2,p} = K_temp;
        end
    end

    % Combining into one matrix (rows = regions, columns = patients).

    K = cell2mat(K);
end