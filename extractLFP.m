function [lfpByChannel, allPowerEst, allPowerEstByBand,F, allPowerVar] = extractLFP(option,lfpFilename, lfpFs, nChansInFile, freqBand)

% Move mandatory info from options
BinWidth = options.BinWidth;
AnalysisTimeWindow = options.AnalysisTimeWindow;
importMode = options.importMode;
EPHYS_DATAPATH = options.EPHYS_DATAPATH;
PERIPHERALS_DATAPATH = options.PERIPHERALS_DATAPATH;
EYEDATA_DATAPATH = options.EYEDATA_DATAPATH;
TRIALDATA_DATAPATH = options.TRIALDATA_DATAPATH;


end