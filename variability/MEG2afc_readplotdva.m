%% load data, 

[megdatraw] = MEG2afc_mergedva( 'raw' ); % 'raw' or BLC

%% do stats 
[megdatraw] = MEG2afc_mergedva_stats(megdatraw);

%% and plot
MEG2afc_mergedva_plotstats(megdatraw)