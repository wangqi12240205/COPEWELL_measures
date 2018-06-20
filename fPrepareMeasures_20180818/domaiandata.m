function domaiandata
Name = {'small_2010.csv',...
        'small_2010domain_maxcurrent.csv',...
        'small_2010domain_nocomp.csv',...
        'small_2010domain_composite.csv',...
        'big_2010.csv',...
        'big_2010domain_maxcurrent.csv',...
        'big_2010domain_nocomp.csv',...
        'big_2010domain_composite.csv',};
    for i = 1:8
    filename = Name{i};
T = readtable(filename);
   
fips = T(:,1);
for i = 1:9
    temp = T(:,[1 i+4]);
    temp.Properties.VariableNames{2} = 'value';
    name = T.Properties.VariableNames {i+4};
    writetable(temp,['mapsdata/' name filename]);
end
    end
keyboard

end