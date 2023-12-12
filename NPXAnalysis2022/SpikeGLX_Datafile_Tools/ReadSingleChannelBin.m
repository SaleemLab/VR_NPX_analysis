
% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0.
%
% In this case, read only for one selected channel - input channel (base 1
% ie channel number does not start at 0 but at 1)
%
% The returned array has dimensions [1,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.
%
% History
%       SGS 14/4/2022 Adapted from ReadBin
%
function dataArray = ReadSingleChannelBin(ch, meta, binName, path)
% Get total number of channels
nChan = str2double(meta.nSavedChans);
% Number of samples is size of file / 2*nChan
nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
nFileSamp = 100*60; % 100s * 60 Hz
% Establish size of return array (desired channels / nSamp)
sizeA = [1, nFileSamp];
skipper = 2 * nChan * (30000/60); % Skip 30000/60 = 500 samples so that only 60 Hz so skip 2bytes*nchan*samples
% Add appropriate index to samp0 to align with this channel
samp0 = 2*(ch-1);
fid = fopen(fullfile(path, binName), 'rb');
fseek(fid,samp0, 'bof');
dataArray = fread(fid, sizeA, 'int16=>double',skipper);
fclose(fid);
end % ReadBin
