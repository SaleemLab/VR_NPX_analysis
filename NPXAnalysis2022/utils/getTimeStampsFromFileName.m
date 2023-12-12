function timestamp = getTimeStampsFromFileName(filelist)
% Datestamp in filename is in format '...*2022-02-22T13_30_30*...'
% This gets date
dateNames = regexp(filelist,'(?<year>\d+)-(?<month>\d+)-(?<day>\d+)','match');
% This gets time
timeNames = regexp(filelist,'(?<hour>\d+)_(?<min>\d+)_(?<sec>\d+)','match');
% Concatenate into one time
timestamp = datenum([dateNames{1}, '-', timeNames{1}],'yyyy-mm-dd-HH_MM_SS');
end
