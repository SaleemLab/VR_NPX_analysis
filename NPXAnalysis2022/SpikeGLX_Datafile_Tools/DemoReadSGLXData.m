% =============================================================
% Simple helper functions and MATLAB structures demonstrating
% how to read and manipulate SpikeGLX meta and binary files.
%
% The most important part of the demo is ReadMeta().
% Please read the comments for that function. Use of
% the 'meta' structure will make your data handling
% much easier!
%
% function DemoReadSGLXData()

% Ask user for binary file
[binName,path] = uigetfile('*.bin', 'Select Binary File');

% Parse the corresponding metafile
meta = ReadMeta(fullfile(path,binName));
nChan = str2double(meta.nSavedChans);
nSamp = floor(10.0 * SampRate(meta));
dataType = 'D';         %set to 'A' for analog, 'D' for digital data

% For an analog channel: gain correct saved channel ch (1-based for MATLAB).
ch = 1;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;

% Read these lines in dw (0-based).
% For 3B2 imec data: the sync pulse is stored in line 6.
% May be 1 or more line indices.
dLineList = [0,1,6];

if dataType == 'A'
    if strcmp(meta.typeThis, 'imec')
        dataArray = GainCorrectIM(dataArray, [ch], meta);
    else
        dataArray = GainCorrectNI(dataArray, [ch], meta);
    end
    plot(dataArray(ch,:));
else
    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    nIter = ceil(nFileSamp/nSamp);
    tic
    % Get first one second of data
    tsync = zeros(1,nFileSamp);
    for thisFilePart = 0:nIter-1
        thisSampleStart = thisFilePart*nSamp;
        dataArray = ReadBin(thisSampleStart, nSamp, meta, binName, path);
        digArray = ExtractDigital(dataArray, meta, dw, dLineList);
        tsync(thisSampleStart+1:thisSampleStart+size(digArray,2)) = digArray(3,:);
    end
    toc
%     for i = 1:numel(dLineList)
%         plot(digArray(i,:));
%         hold on
%     end
%     hold off
end
%end % DemoReadSGLXData




