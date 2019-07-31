%Neurotar_headDirectionAnalysis

hdp = HeadDirectionPreprocessor(data,floating);


hda = HeadDirectionAnalysis(hdp.workingData.data,hdp.workingData.floating);


hda.smoothResponse(3)

processedDFF2 = hda.getProcessedDFF;