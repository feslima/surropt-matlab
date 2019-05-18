function [xSampled,fObjSampled,gSampled] = sampleModel(x, sample_function,dataForSurrogate)

sampledData = sample_function(x);
dataForSurrogate.doeBuild = sampledData;
[xSampled,fObjSampled,gSampled] = distributeDataForSurrogate(dataForSurrogate);