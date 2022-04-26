function int = computeintegral(p)

int_well = nan(158,5,3);
int_ill = nan(158,5,3);

for p = 1:3
for i = 3:7
int_well(:,i-2,p) = cumtrapz(data{p}(:,1)',data{p}(:,i));
end
end

for p = 4:6
for i = 3:7
int_ill(:,i-2,p-3) = cumtrapz(data{p}(:,1)',data{p}(:,i));
end
end

int =[ int_well(end,1:5,1)', int_well(end,1:5,2)', int_well(end,1:5,3)', int_ill(end,1:5,1)',int_ill(end,1:5,2)',int_ill(end,1:5,3)'];
