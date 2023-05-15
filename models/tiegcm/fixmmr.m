% Correct total mixing ratio after assimilation
% see https://github.com/NCAR/DART/issues/469
%   input/output : tiegcm primary restart file
%         
function fixmmr(file)

O1 = ncread(file, 'O1');
O2 = ncread(file, 'O2');
HE = ncread(file, 'HE');


total = O1 + O2 + HE;
index = (total>1);

O1(index) = O1(index)./total(index);
O2(index) = O2(index)./total(index);
HE(index) = HE(index)./total(index);

fixed = sum(index(:));
fprintf('Number of total mmr fixed: %d\n', fixed);

if fixed > 0
    ncwrite(file, 'O1', O1);
    ncwrite(file, 'O2', O2);
    ncwrite(file, 'HE', HE);
end

end