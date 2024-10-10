classdef NPadmin
    %NPADMIN Neuropixel helper functions
    %   Just a place to hold all the function involved in helping to
    %   perform all the admin for reading in the neuropixel data
    
    properties
    end
    
    methods (Static)
        % Meta functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [meta] = ReadMeta(binName, path)
            
            % Create the matching metafile name
            [dumPath,name,dumExt] = fileparts(binName);
            metaName = strcat(name, '.meta');
            
            % Parse ini file into cell entries C{1}{i} = C{2}{i}
            fid = fopen(fullfile(path, metaName), 'r');
            % -------------------------------------------------------------
            %    Need 'BufSize' adjustment for MATLAB earlier than 2014
            %    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
            C = textscan(fid, '%[^=] = %[^\r\n]');
            % -------------------------------------------------------------
            fclose(fid);
            
            % New empty struct
            meta = struct();
            
            % Convert each cell entry into a struct entry
            for i = 1:length(C{1})
                tag = C{1}{i};
                if tag(1) == '~'
                    % remake tag excluding first character
                    tag = sprintf('%s', tag(2:end));
                end
                meta = setfield(meta, tag, C{2}{i});
            end
            
            % convert strings to doubles for useful params
            meta.nSavedChans   = str2double(meta.nSavedChans);
            meta.fileTimeSecs  = str2double(meta.fileTimeSecs);
            meta.fileSizeBytes = str2double(meta.fileSizeBytes);
            if isfield(meta,'imSampRate')
                meta.imSampRate = str2double(meta.imSampRate);
            end
            if isfield(meta,'niSampRate')
                meta.niSampRate = str2double(meta.niSampRate);
            end
                    
            meta.nFileSamp = meta.fileSizeBytes / (2 * meta.nSavedChans);
        end % ReadMeta
        
        function [tableconfig, type] = getNPChannelConfig(meta)
            
            warning('off', 'MATLAB:table:RowsAddedExistingVars');
            warning('off','MATLAB:table:RowsAddedNewVars');
            
            config = NPadmin.IMROtoElectrodes(meta.imroTbl); % Get the mapping
            type   = config.type;
            % Get saved channels
            [spikeGLXchans,~,recorded_idx] = NPadmin.convert_snsSaveChanSubset_to_channels(meta);
            ks = NPadmin.electrodesToKcoords(config); % get x/y coordinates (used in ks)           
            
            % Make table
            tableParams = {...
                NaN,    'Channel';...
                NaN,    'Bank';...
                NaN,    'SpikeGLXchan0';...
                NaN,    'Shank';...
                NaN,    'Electrode';...
                false,  'Recorded_idx';...
                false,  'Internal_reference';...
                NaN,    'Ks_xcoord';...
                NaN,    'Ks_ycoord';...
                };
            
            tableconfig = table(tableParams{:,1});
            tableconfig.Properties.VariableNames = tableParams(:,2);
                        
            tableconfig.Channel(1:384,1)       = config.channels; % 1-indexed
            tableconfig.Bank(1:384,1)          = config.banks;
            tableconfig.SpikeGLXchan0(1:384,1) = spikeGLXchans;   % 0-indexed
            tableconfig.Shank(1:384,1)         = config.shanks;
            tableconfig.Electrode(1:384,1)     = config.elecs;    % 1-indexed
            tableconfig.Recorded_idx(1:384,1)  = recorded_idx;
            tableconfig.Ks_xcoord(1:384,1)     = ks.xcoords;
            tableconfig.Ks_ycoord(1:384,1)     = ks.ycoords;
                   
            % Internal reference channel 
            switch type
                case 0
                    refsites = [192, 576,960]; % 1-indexed
                    tableconfig.Internal_reference(ismember(tableconfig.Electrode,refsites),1)= true;
                case 24
                    refchannel = 127; % 0-indexed
                    tableconfig.Internal_reference(ismember(tableconfig.SpikeGLXchan0,refchannel),1)= true;
            end
            
            warning('on', 'MATLAB:table:RowsAddedExistingVars');
            warning('on','MATLAB:table:RowsAddedNewVars');
            
        end % getNPChannelConfig
        
        % Mapping functions
        function kcoords = electrodesToKcoords(config)
            if config.type == 0
                % input unique electrode site number for neuropixels 1.0 (aka phase 3B)
                % returns channel number, and x/y coordinates on probe
                %
                % mapping based on kilosort standard:
                %
                %   :   :   :
                %   |7  |8  |   ^
                %   |  5|  6|   |  60
                %   |3  |4  |   |  40
                %   |  1|  2|   |  20
                %    \     /    |
                %     \   /     |
                %      \ /      |
                %               |
                %  <------------
                %   59 43 27 11
                %
                %
                % Soraya Dunn 2019
                
                
                for n = 1:numel(config.elecs)
                    electrode = config.elecs(n);
                    row = ceil(electrode/2);
                    if electrode <= 384
                        bank = 0;
                    elseif and(electrode>384, electrode<=384*2)
                        bank = 1;
                    elseif and(electrode>384*2, electrode<=960)
                        bank = 2;
                    else
                        keyboard
                    end
                    channel = electrode - bank*384;  % starting at 1
                    
                    kcoords.chanMap(n,1) = channel;
                    
                    if isOdd(electrode)
                        if isOdd(row)
                            kcoords.xcoords(n,1) = 43;
                        else
                            kcoords.xcoords(n,1) = 59;
                        end
                    else
                        if isOdd(row)
                            kcoords.xcoords(n,1) = 11;
                        else
                            kcoords.xcoords(n,1) = 27;
                        end
                        
                    end
                    kcoords.ycoords(n,1) = row * 20;
                end
                
            elseif config.type == 24 | config.type == 2013 | config.type == 2014
                % mapping based on NP1.0 kilosort standard and adjust using
                %     Fig 1 of https://www.biorxiv.org/content/10.1101/2020.10.27.358291v1.full.pdf:
                %
                %   :   :   :
                %   |7  |8  |   ^
                %   |5  |6  |   |  45
                %   |3  |4  |   |  30
                %   |1  |2  |   |  15
                %    \     /    |
                %     \   /     |
                %      \ /      |
                %               |
                %  <------------
                %   59   27
                %
                %
                % Katarina Poole 2021
                %     center-to-center spacing between shanks on the four-shank version is 250 µm
                for n = 1:numel(config.elecs)
                    electrode = config.elecs(n);
                    shank = config.shanks(n);
                    kcoords.chanMap(n,1) = config.channels(n);
                    row = ceil(electrode/2);
                    kcoords.ycoords(n,1) = row * 15;
                    if isOdd(electrode)
                        kcoords.xcoords(n,1) = 59 + ((shank-1)*250);
                    else
                        kcoords.xcoords(n,1) = 27 + ((shank-1)*250);
                    end
                end
            else
                warndlg('This type cannot be loaded, please pick another file')
                return
            end
        end
        
        function config = IMROtoElectrodes(IMRO)
            % Gets electrode info from IMRO table
            % Katarina Poole 2021
            elecs = [];
            shanks = [];
            banks = [];
            channels = [];
            % Remove first and last bracket
            IMRO = IMRO(2:end-1);
            IMROcell = split(IMRO,')(');
            % Get the header
            header = split(IMROcell{1,1},',');
            if str2double(header{1,1}) == 24 | str2double(header{1,1}) == 2013 | str2double(header{1,1}) == 2014
                IMROcell = IMROcell(2:end,1);
                for n = 1:size(IMROcell,1)
                    currRow = str2num(IMROcell{n,1})+1; % Index from 1
                    channels(n) = currRow(1);
                    shanks(n) = currRow(2);
                    banks(n) = currRow(3);
                    elecs(n) = currRow(5);
                    refs(n) = currRow(4);
                end
            elseif str2double(header{1,1}) == 0
                IMROcell = IMROcell(2:end,1);
                for n = 1:size(IMROcell,1)
                    currRow = str2num(IMROcell{n,1})+1;
                    channels(n) = currRow(1);
                    shanks(n) = 1;
                    banks(n) = currRow(2)-1;
                    elecs(n) = channels(n) + banks(n)*384;
                    refs(n) = currRow(3);
                end
            else
                warndlg('This type cannot be loaded, please pick another file')
                return
            end
            config.channels = channels;
            config.shanks = shanks;
            config.banks = banks;
            config.elecs = elecs;
            config.type = str2double(header{1,1});
            reftype = unique(refs);
            %             if numel(reftype)==1
            config.reference = reftype;
%             else
%                 keyboard % not quite sure what to do in the case of different channels having different references... but that shouldn't happen the way we usually record!
%             end
        end
        
        function shanks = xcoordsToShank(xcoords)
            shanks = floor(xcoords/250)+1;
        end
        
        % Channel admin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [sGLXchans,channels,recorded_idx] = convert_snsSaveChanSubset_to_channels(meta)
            
            %
            % NP 1.0 spikeGLX saved channel assignments: AP channels 0:383, LF channels 384:767, sync: 768
            % returns:
            % all recorded channels with spikeGLX indices (excl sync), NaN for missing channels
            % recorded channels starting at index 1 (for both LF and AP channels)
            % logical index of recorded channels (vector of length 384)
            %
            % Soraya Dunn 2019
            % Edited by Katarina Poole 2021 to work with NP2.0
            
            sGLXchans = transpose(str2num(meta.snsSaveChanSubset));
                  
            switch meta.imDatPrb_type
                case '0' % For 1.0 probes
                    sGLXchans(sGLXchans == 768) = []; % remove sync, only neural channels left
                    
                    if round(meta.imSampRate) == 2500 % LF channels
                        set1 = 384:767;   % LF channels
                        missing_chans = find(~ismember(set1,sGLXchans));
                        for n = 1:numel(missing_chans)
                            sGLXchans = [sGLXchans(1:missing_chans(n)-1); NaN; sGLXchans(missing_chans(n):end)];
                        end
                        channels = sGLXchans - 383; % channels 1-indexed
                        
                    elseif round(meta.imSampRate) == 30000 % AP channels
                        set1 = 0:383;  % AP channels
                        missing_chans = find(~ismember(set1,sGLXchans));
                        for n = 1:numel(missing_chans)
                            sGLXchans = [sGLXchans(1:missing_chans(n)-1); NaN; sGLXchans(missing_chans(n):end)];
                        end
                        channels = sGLXchans + 1;
                    end
                otherwise % 2.0 probes (do not have separate LF channels)
                    sGLXchans(sGLXchans == 384) = []; % remove sync, only neural channels left
                    set1 = 0:383;  % AP channels
                    missing_chans = find(~ismember(set1,sGLXchans));
                    for n = 1:numel(missing_chans)
                        sGLXchans = [sGLXchans(1:missing_chans(n)-1); NaN; sGLXchans(missing_chans(n):end)];
                    end
                    channels = sGLXchans + 1;
               
            end
            
            recorded_idx = ~isnan(channels);
            
        end
        
        function chans = OriginalChans(meta)
            % Return array of original channel IDs. As an example,
            % suppose we want the imec gain for the ith channel stored
            % in the binary data. A gain array can be obtained using
            % ChanGainsIM() but we need an original channel index to
            % do the look-up. Because you can selectively save channels
            % the ith channel in the file isn't necessarily the ith
            % acquired channel, so use this function to convert from
            % ith stored to original index.
            %
            % Note: In SpikeGLX channels are 0-based, but MATLAB uses
            % 1-based indexing, so we add 1 to the original IDs here.
            if strcmp(meta.snsSaveChanSubset, 'all')
                chans = (1:str2double(meta.nSavedChans));
            else
                chans = str2num(meta.snsSaveChanSubset);
                chans = chans + 1;
            end
        end % OriginalChans
        
        function [AP,LF,SY] = ChannelCountsIM(meta)
            % =========================================================
            % Return counts of each imec channel type that compose
            % the timepoints stored in binary file.
            %
            M = str2num(meta.snsApLfSy);
            AP = M(1);
            LF = M(2);
            SY = M(3);
        end % ChannelCountsIM
        
        function [internalRefs] = getInteralRefElectrodes(type)
            % 1-indexed electrode numbers
            switch type
                case 0
                    internalRefs = [192,576,960];
                case 21
                    internalRefs = [127,507,887,1251];
                case 24
                    internalRefs = [127,511,895,1279];
            end
            
        end
        
        function NPMAP = getNPMAP(type)
            parentPath = fileparts(fileparts(which('NPadmin')));
            mapDir = [parentPath,'\utils\tables\'];
            fileName = 'Neuropixels_Electrode_Mapping.xlsx';
            fileDir = fullfile(mapDir,fileName);
            switch type % Load in the mapping file
                case 24                   
                    NPMAP.Shank1 = readtable(fileDir,'Sheet','Multi-Shank 1 (SR_CHAIN1)');
                    NPMAP.Shank2 = readtable(fileDir,'Sheet','Multi-Shank 2 (SR_CHAIN2)');
                    NPMAP.Shank3 = readtable(fileDir,'Sheet','Multi-Shank 3 (SR_CHAIN3)');
                    NPMAP.Shank4 = readtable(fileDir,'Sheet','Multi-Shank 4 (SR_CHAIN4)');
                case 0
                    NPMAP.Shank1 = readtable(fileDir,'Sheet','1.0 Single Shank (SR_CHAIN1)');
                case 21
                    NPMAP.Shank1 = readtable(fileDir,'Sheet','2.0 Single Shank (SR_CHAIN1)');
            end
        end
        
        function [bHv] = getBhvInfo(binName,binPath)
            
        end
        
        % Gain adjust %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dataArray = GainCorrectIM(dataArray, chanList, meta)
            
            %
            % Having acquired a block of raw imec data using ReadBin(),
            % convert values to gain-corrected voltages. The conversion
            % is only applied to the saved-channel indices in chanList.
            %
            % Output in uV
            %
            
            
            % Look up gain with acquired channel ID
            chans = NPadmin.OriginalChans(meta);  % chans output is 1-indexed
            [APgain,LFgain] = NPadmin.ChanGainsIM(meta);
            nAP = length(APgain);
            nNu = nAP * 2;
            
            % Common conversion factor
            fI2V = NPadmin.Int2Volts(meta);
            
            for i = 1:length(chanList)
%                 j = chanList(i);    % index into timepoint
%                 k = chans(j);       % acquisition index
                k = chanList(i);       % acquisition index
                if k <= nAP
                    conv = fI2V / APgain(k);
                elseif k <= nNu
                    conv = fI2V / LFgain(k - nAP);
                else
                    continue;
                end
                dataArray(i,:) = dataArray(i,:) * conv; % in V % changed index from j to i as getting indexing error when putting channels in that don't start at 1, think it's an error in original spikeGLX code               
            end
            dataArray = dataArray * 10^6; % in uV - useful if you want to save LFP data as int16
        end
            
        function [APgain,LFgain] = ChanGainsIM(meta)
            % Return gain arrays for imec channels (from 06/01/22 spikeGLX_Datafile_tools http://billkarsh.github.io/SpikeGLX/#post-processing-tools).
            %
            % Index into these with original (acquired) channel IDs.
            %      
            
            if isfield(meta,'imDatPrb_type')
                probeType = str2num(meta.imDatPrb_type);
            else
                probeType = 0;
            end
            if (probeType == 21) || (probeType == 24)
                [AP,LF,~] = NPadmin.ChannelCountsIM(meta);
                % NP 2.0; APgain = 80 for all channels
                APgain = zeros(AP,1,'double');
                APgain = APgain + 80;
                % No LF channels, set gain = 0
                LFgain = zeros(LF,1,'double');
            else
                % 3A or 3B data?
                % 3A metadata has field "typeEnabled" which was replaced
                % with "typeImEnabled" and "typeNiEnabled" in 3B.
                % The 3B imro table has an additional field for the
                % high pass filter enabled/disabled
                if isfield(meta,'typeEnabled')
                    % 3A data
                    C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
                        'EndOfLine', ')', 'HeaderLines', 1 );
                else
                    % 3B data
                    C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
                        'EndOfLine', ')', 'HeaderLines', 1 );
                end
                APgain = double(cell2mat(C(1)));
                LFgain = double(cell2mat(C(2)));
            end
        end % ChanGainsIM
        
        function fI2V = Int2Volts(meta)
            % Return a multiplicative factor for converting 16-bit
            % file data to voltage. This does not take gain into
            % account. The full conversion with gain is:
            %
            %   dataVolts = dataInt * fI2V / gain.
            %
            % Note that each channel may have its own gain.
            % Updated 06/01/22
            if strcmp(meta.typeThis, 'imec')
                if isfield(meta,'imMaxInt')
                    maxInt = str2num(meta.imMaxInt);
                else
                    maxInt = 512;
                end
                fI2V = str2double(meta.imAiRangeMax) / maxInt;
            else
                fI2V = str2double(meta.niAiRangeMax) / 32768;
            end
        end % Int2Volts
               
        function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
            % =========================================================
            % Return gain for ith channel stored in the nidq file.
            %
            % ichan is a saved channel index, rather than an original
            % (acquired) index.
            %
            if ichan <= savedMN
                gain = str2double(meta.niMNGain);
            elseif ichan <= savedMN + savedMA
                gain = str2double(meta.niMAGain);
            else
                gain = 1;
            end
        end % ChanGainNI
        
        
        % file paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [p] = getBinPath(rootPath,binName)
            % returns full file path of bin file in root directory (or
            % whatever point in the dir tree above the file)
            % returns empty vector and displays warning if no file found
            binDir = dir(fullfile(rootPath,['**\' binName]));
            if isempty(binDir)
                p = [];
                disp(['Could not find ' binName ' in ' rootPath])
            else
                p      = binDir.folder;
            end
        end
        
        % data extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % =========================================================
        
        function dataArray = ReadBin(samp0, nSamp, meta, binName, path)
            % From spikeGLX offline tools
            % Read nSamp timepoints from the binary file, starting
            % at timepoint offset samp0. The returned array has
            % dimensions [nChan,nSamp]. Note that nSamp returned
            % is the lesser of: {nSamp, timepoints available}.
            %
            % IMPORTANT: samp0 and nSamp must be integers.
            %
            nChan = meta.nSavedChans; 
            
            nFileSamp = meta.fileSizeBytes / (2 * nChan);
            samp0 = max(samp0, 0);
            nSamp = min(nSamp, nFileSamp - samp0);
            
            sizeA = [nChan, nSamp];
            
            fid = fopen(fullfile(path, binName), 'rb');
            fseek(fid, samp0 * 2 * nChan, 'bof');
            dataArray = fread(fid, sizeA, 'int16=>int16');
            fclose(fid);
        end % ReadBin
        
        
        function [tline1k,tline] = makeTimelines(meta)
            tline   = transpose(1/meta.imSampRate : 1/meta.imSampRate : meta.fileTimeSecs);
            tline1k = transpose(tline(1) : 1/1000 : tline(end));
        end
        
        function chanDat = loadOneChanfromBin(filename,numChans,loadChan)
            % based on spikerepo github extractSyncChannelFromFile
            % extracts as int16 (same format as data is saved)
            maxReadSize = 1e5;
            
            fid = fopen(filename, 'r');
            
            % skip over the first samples of the other channels
            try
                q = fread(fid, (loadChan-1), 'int16=>int16');
            catch
                chanDat = []; % retunr empty vector if file could not be read
                return
            end
            
            d = dir(filename);
            nSamp = d.bytes/2/numChans;
            
            chanDat = zeros(1, nSamp, 'int16');
            
            nBatch = floor(nSamp/maxReadSize);
            for b = 1:nBatch
                dat = fread(fid, [1, maxReadSize], 'int16=>int16', (numChans-1)*2); % skipping other channels
                chanDat((b-1)*maxReadSize+1:b*maxReadSize) = dat;
            end
            % all the other samples
            dat = fread(fid, [1, Inf], 'int16=>int16', (numChans-1)*2); % skipping other channels
            
            chanDat(nBatch*maxReadSize+1:end) = dat;
            
            fclose(fid);
        end
        
        function digArray = extract_NI_digital(digChan, dLineList)
            
            %
            % Return an array [lines X timepoints] of uint8 values for
            % a specified set of digital lines.
            
            % dLineList is 0-based (ie will be within 0:7)
            
            [~,nSamp] = size(digChan);
            digArray = zeros(numel(dLineList), nSamp, 'uint8');
            
            for i = 1:numel(dLineList)
                digArray(i,:) = bitget(digChan, dLineList(i)+1, 'int16');
            end
            
        end
        
        function [trialTimes, SyncTimes] = getTrialTimes(binName,binPath) % specific to KP
            %NPTRIALTIMES Gets trial times from sync chan
            %   Detailed explanation goes here
            
            % Set params
            nd = 6; % I think its the 6th channel for the sync and 1 is the imec
            syncimec = 1;
            % Look at nidq.bin
            NIbinName = [binName(1:strfind(binName,'.imec')-1),'.nidq.bin'];
            d = dir(fullfile(binPath, '..'));
            NIbinPath = [d(1).folder,'\'];
            meta     = NPadmin.ReadMeta(NIbinName, NIbinPath);
            digChan = NPadmin.loadOneChanfromBin(fullfile(NIbinPath,NIbinName),meta.nSavedChans,meta.nSavedChans);
            dLine    = nd-1;
            syncLine = syncimec-1;
            digArray = NPadmin.extract_NI_digital(digChan,dLine);
            syncArray = NPadmin.extract_NI_digital(digChan,syncLine);
            
            logIdx = logical(digArray);
            syncPulse = logical(syncArray);
            NStartSamples = strfind(logIdx,[0 1])+1;
            NSyncSamples = strfind(syncPulse,[0 1])+1;
            Nfs = meta.niSampRate;
            trialTimes = (NStartSamples/Nfs)+0.05; % add 0.05 since the pulse starts 0.05 before stim onset
            SyncTimes = (NSyncSamples/Nfs);
        end 
        
    end %methods
end %classdef
