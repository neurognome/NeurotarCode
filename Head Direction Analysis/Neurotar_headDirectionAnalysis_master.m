%Neurotar_headDirectionAnalysis

preprocessed_data = HeadDirectionPreprocessor(data,floating);

processedDFF = preprocessed_data.getProcessedDFF;

hda = HeadDirectionVisualize(processedDFF);



hda.smoothResponse(3)

processedDFF2 = hda.getProcessedDFF;