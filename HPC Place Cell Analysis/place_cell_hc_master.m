hpc = HPC_PlaceCellPreprocessor(data,floating)

data = hpc.getImagingData;
floating = hpc.getNeurotarData;

HPC_PlaceCellAnalyzer(data,floating)