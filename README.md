# COPEWELL_measures
Measures justification

original csv files:
list_2016_12_13_deleteduplicate.csv   : list of current measures without duplication in new measures
list_2018_06_05.csv.                  : list of new measures
measures_current2015.csv              : dataset of current measures (scaled)
measures2015.csv   				      : dataset of new measures (scaled)

step1 : run fMeasuresjustification.m to optimize signs of new measures and combined with current measures to maximize cronbach's alpha.
input files : 4 original csv files
output files: 3 mat files:
dataplot_original.mat               : data for plot correlation plots on original signs 
dataplot_currentmaxv2.mat           : data for plot correlation plots on optimized signs 
fPrepareMeasuresQiResults_180620.mat: data for combination plots 

step2 : run ffinalchoicedata.m 
input: fPrepareMeasuresQiResults_180620.mat
output: generatet the mat file:
dataplot_finalchoice.mat            : data for plot correlation plots of final choiced combination

step3: run fPlotResultsQi.m to generate the combination plots (grid plots)
input: fPrepareMeasuresQiResults_180620.mat 
output: grided combination plots, saved in combination folder

step4: run fplot.m 
input: dataplot_original.mat               
       dataplot_currentmaxv2.mat
       dataplot_finalchoice.mat

output: three correlation plots on 
1. original signs, save in folder figure_original
2. optimazed signs, save in folder figure_maxcurrent
1. combination with current measures, save in folder figure_finalchoice



