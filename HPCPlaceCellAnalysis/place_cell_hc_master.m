hpc = HPC_PlaceCellPreprocessor(data,floating,'ForceTimeLock',1);

% this is just to match Will's data, and speed up processing time...
hpc.setForceTimeLock(true);

hpc.processData;

% Important notes about the HPC_PlaceCellAnalyzer:

% The following class is dependent on another class I created that I call a
% dataObject. The purpose of the dataObject class is to make it easy to store
% and retrieve large numbers of variables between functions, while also
% providing a nice way of organizing the data as you do analyses, to prevent 
% your workspace getting filled up with tons of variables. In this code
% I initialize 3 instances of the dataObject class: analysisData, plottingData,
% and workingData. As all three objects are instances of the dataObject class,
% they ALL can use the methods of dataObject. this includes general stuff
% like being able to add and remove properties from the dataObject. Or being
% able to export specific properties from the dataObject to the workspace.
% 
% You can think of this like a structure. I created a function (the dataObject
% class) which stores whatever variables you want into the structure (the 
% dataObject instance). Then, you can use other functions (dataObject methods)
% to add or retrieve data from the structure. The main benefit is better
% handling of these data, and more flexibility in the way you can work with it.
% I'll have some examples below.
% 
% Importantly, I'm trying to get some kind of standardization for my own code going
% on, so all dataObject classes within the analysis have a hierarchy. 
% 
% The top of the hierarchy is analysisData. Most of the output of your analyses 
% will go into analysisData. This includes stuff like boolean identifiers
% (isPlaceCell), calculated metrics (spatial_information), and others. Note that 
% I included some things like the spatial_information variable, even though we 
% will use it for plotting. I figured that the actual values of the spatial information
% would be useful for beyond plotting, so I organized it into analysisData.
% 
% Next is plottingData. These data consist of things that are created for the SOLE
% purpose of visualization. For example, the heatmaps (I) aren't used in much further
% analysis, but are used to visualize place cells. They don't hold much analysis
% value, but are necessary for visualization, so they go into plottingData.
% 
% Lastly, is workingData. These are variables that may take  along time to create
% so you only want to run that analysis once... or they're useful variables
% that need to be passed around often, but will never be seen or used outside
% of creating variables for analysis or plotting.
% 
% The hierarchy works in this way, if any data has any chance of being
% part of analysisData, it goes in there. Then, plottingData takes
% precedence over workingData, which is mainly just to temporarily store junk
    

hpa = HPC_PlaceCellAnalyzer(hpc.getImagingData,hpc.getNeurotarData);

% Now that all the preprocessing is finished, you can do lots of stuff...

% Find place cells and display plot:

hpa.findPlaceCells;

% LOok at a specific cell's heat map:
cellID = 20;
hpa.makeHeatMap(cellID);

% Let's say that we want to get the isPlaceCell boolean vector from
% our data so that we can maybe store it or use it elsewhere. We can
% call the exportData method of our dataObject to spit out specific variables
% into our workspace

hpa.analysisData.exportData('isPlaceCell')
