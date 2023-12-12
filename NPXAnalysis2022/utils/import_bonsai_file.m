% Function to import arbitrary data from bonsai log files
% SGS 03/04/2022
function [data] = import_bonsai_file(filename,delimiter, startRow, endRow)

% Initialize variables
if nargin < 2
    delimiter = ',';
end
if nargin < 3
    startRow = 2;
    endRow = inf;
end

% Open the photodiode text file.
fileID = fopen(filename,'r');
% Check to see how many columns
tl = fgetl(fileID);
frewind(fileID);
Cols = textscan(tl,'%s','delimiter',','); 
Cols = Cols{1};
nCols = length(Cols);
% Format for each line of text:
formatSpec = [repmat('%f',1,nCols),'%[^\n\r]'];

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Bonsai can return '.' in the column header - replace with underscore if
% present
for thisCol = 1:length(Cols)
    Cols{thisCol} = strrep(Cols{thisCol},'.','_');
end

% Allocate imported array to column variable names
for thisCol = 1:length(Cols)
    eval(sprintf('data.%s = dataArray{:,%01d};',Cols{thisCol},thisCol))
end




