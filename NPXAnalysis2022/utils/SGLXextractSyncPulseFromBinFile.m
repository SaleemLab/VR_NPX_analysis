% Program to try to extract the sync pulse from NP1 or NP2 data sets
%   This function supercedes extractSyncPulseFromAPFile, which used the
%   Neuropixels memmap style, and instead uses Jennifer Collonel's toolbox
% Input - path to binary file
% Output - times of onset and offset of the sync pulse in a structure
%
% History: SGS 14th April 2022 Wrote it
%
function syncTimes_ephys = SGLXextractSyncPulseFromBinFile(binpath)

% Extract path and filename
[fpath,binName,extension] = fileparts(binpath);
binName = [binName,extension];
% Parse the corresponding metafile
meta = ReadMeta(binpath);
nChan = str2double(meta.nSavedChans);
nSamp = floor(10.0 * SampRate(meta)); % Read in 10s interval chunks - could improve this
nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;
% Read these lines in dw (0-based).
% For 3B2 imec data: the sync pulse is stored in line 6.
% May be 1 or more line indices.
dLineList = [0,1,6];

nIter = ceil(nFileSamp/nSamp);
% Get first one second of data
tsync = zeros(1,nFileSamp);
for thisFilePart = 0:nIter-1
    thisSampleStart = thisFilePart*nSamp;
    dataArray = ReadBin(thisSampleStart, nSamp, meta, binName, fpath);
    digArray = ExtractDigital(dataArray, meta, dw, dLineList);
    tsync(thisSampleStart+1:thisSampleStart+size(digArray,2)) = digArray(3,:);
end


%%%
% Store the pulse changes as time stamps
% Find upswings in sync pulse
vec_idx = find(diff(tsync)>=0.5);
% add 1 to vec_idx to compensate for diff and convert from samples to s 
syncTimes_ephys.on = (vec_idx+1)./SampRate(meta); 
% Find downswings in sync pulse
vec_idx = find(diff(tsync)<=-0.5); 
% add 1 to vec_idx to compensate for diff and convert from samples to s 
syncTimes_ephys.off = (vec_idx+1)./SampRate(meta); 
syncTimes_ephys.Sync = tsync;
