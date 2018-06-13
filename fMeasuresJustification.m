function fMeasuresJustification
% 2018 June 13
% Qi Wang


%% similar as fPrepareMeasures.mat)
rawdata = readtable('measures2015.csv', 'ReadVariableNames',true);
Dat0 = rawdata{:,2:end};
Dat = readtable('measures_current2015.csv', 'ReadVariableNames',true);
Dat = Dat{:,2:end};
measureList = readtable('list_2016_12_13.csv');
 