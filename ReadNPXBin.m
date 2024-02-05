
% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [nChan,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.

% Modified by Masa to remove str2double as variables needed from meta file 
% is already in double format
function dataArray = ReadNPXBin(samp0, nSamp, meta, binName, path)

    nChan = meta.nSavedChans;

    nFileSamp = meta.fileSizeBytes / (2 * nChan);
    samp0 = max(samp0, 0);
    nSamp = min(nSamp, nFileSamp - samp0);
    dataArray = [];

%     if isfield(meta,'selected_channels')
%         % in progress not working
%         sizeA = [1, nSamp];
%         start_sample = samp0 * 2 * nChan;
%         fid = fopen(fullfile(path, binName), 'rb');
% 
%         for n = 100
%             fseek(fid, start_sample + (meta.selected_channels(n)-1)*nSamp, 'bof');
%             dataArray(n,:)= fread(fid, nSamp, 'int16=>double');
%         end
%         fclose(fid);
% 
%     else
%     
        sizeA = [nChan, nSamp];
        start_sample = samp0 * 2 * nChan;
        fid = fopen(fullfile(path, binName), 'rb');
        fseek(fid, start_sample, 'bof');
        dataArray = fread(fid, sizeA, 'int16=>double');
        fclose(fid);
   
%     end
      

end % ReadBin
