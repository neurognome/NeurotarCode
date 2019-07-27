hpc = HPC_PlaceCellPreprocessor(data,floating);

% this is just to match Will's data, and speed up processing time...
hpc.setForceTimeLock(true);
hpc.processData;

hpa = HPC_PlaceCellAnalyzer(hpc.getImagingData,hpc.getNeurotarData);

hpa.findPlaceCells