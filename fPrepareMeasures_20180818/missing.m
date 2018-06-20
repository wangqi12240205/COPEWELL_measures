function missing
T = readtable('data_2014.csv');
Iskip = [23,25,26,28,32:45,49:57,67:69,72:81,95,97,118,120,121,125:128,141,143,152,153];
T_Iskip = T(:,Iskip);
M = (table2array(T_Iskip)==0);


zero = sum(M)';
II = zero>=628&zero<1037;
T_less = T_Iskip(:,II);
T_less.Properties.VariableNames
por = zero(II)./3142;


keyboard
end