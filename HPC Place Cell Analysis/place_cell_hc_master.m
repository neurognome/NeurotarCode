hpc = HPC_PlaceCellPreprocessor(data,floating);
hpc.setForceTimeLock(true);
hpc.processData;

hpa = HPC_PlaceCellAnalyzer(hpc.getImagingData,hpc.getNeurotarData);

hpa.findPlaceCells