
% Write out coordinates for a Neuropixels 3A, 1.0 or 2.0 metadata file.
% Format selected with the outType variable.
% Jennifer Colonell, Janelia Research Campus
%
function SGLXMetaToCoords_ChannelMap_masa(metaName)


% Output selection: 0 for text coordinate file;
%                   1 for Kilosort or Kilosort2 channel map file;
%                   2 for strings to paste into JRClust .prm file
outType = 1;

% Ask user for metadata file

%metaName = dir(fullfile(path, '*.meta'));
%metaName = metaName(1).name; % edd edit ?
% Shank separation for multishank
shankSep = 250;

% Parse in file to get the metadata structure
meta = ReadMeta(fullfile(metaName.folder,metaName.name));

if isfield(meta, 'imDatPrb_type')
    pType = str2num(meta.imDatPrb_type);
else
    pType = 0; %3A probe
end


if isfield(meta,'snsGeomMap')
    [nShank, shankWidth, shankPitch, shankind, xcoords, ycoords, connected] = geomMapToGeom(meta);
elseif isfield(meta,'snsShankMap')
    [nShank, shankWidth, shankPitch, shankind, xcoords, ycoords, connected] = shankMapToGeom(meta);

elseif pType <= 1
    %Neuropixels 1.0 or 3A probe
    [elecInd, connected] = NP10_ElecInd(meta);

    % Get saved channels
    chans = OriginalChans(meta);
    [AP,~,SY] = ChannelCountsIM(meta);
    chans = chans(1:AP);        %1-based channel numbers

    % Trim elecInd and shankInd to include only saved channels
    elecInd = elecInd(chans);
    shankind = zeros(size(elecInd));

    % Get XY coords for saved channels
    [xcoords, ycoords] = XYCoord10(meta, elecInd);
else
    % Parse imro table for shank and electrode indicies
    [elecInd, shankind, bankMask, connected] = NP20_ElecInd(meta);

    % Get saved channels
    chans = OriginalChans(meta);
    [AP,LF,SY] = ChannelCountsIM(meta);
    chans = chans(1:AP);        %1-based channel numbers

    % Trim elecInd and shankInd to include only saved channels
    elecInd = elecInd(chans);
    shankind = shankind(chans);

    % Get XY coords for saved channels
    [xcoords, ycoords] = XYCoord20(meta, elecInd, bankMask, shankind);

end


% Build output name and write out file
[~,fname,~] = fileparts(metaName.name);

switch outType
    case 0      %tab delimited, chan, x, y, shank
        newName = [fname,'-siteCoords.txt'];
        fid = fopen( newName, 'w');
        for i = 1:numel(elecInd)
            currX = shankind(i)*shankSep + xcoords(i);
            fprintf( fid, '%d\t%d\t%d\t%d\n', chans(i)-1, currX, ycoords(i), shankind(i));
        end
        fclose(fid);

    case 1     %KS2 *.mat
        newName = [fname,'_kilosortChanMap.mat'];
        chanMap = (1:numel(xcoords))';
        chanMap0ind = chanMap - 1;
        connected = logical(connected);
        xcoords = shankind*shankSep + xcoords;   %KS2 not yet using kcoords, so x coord includes shank sep
        kcoords = shankind + 1;     %KS1 uses kcoords to force templates to be on one shank
        name = fname;

        save( fullfile(metaName.folder, newName), 'chanMap', 'chanMap0ind', 'connected', 'name', 'xcoords', 'ycoords', 'kcoords' );

    case 2  %strings to copy into JRC prm file
        newName = [fname,'_forJRCprm.txt'];
        nchan = numel(chans);
        fid = fopen( newName, 'w' );
        fprintf( fid, 'shankMap = [' );
        for i = 1:nchan-1
            fprintf( fid, '%d,', shankind(i) + 1 ); % switch to 1-based for MATLAB
        end
        fprintf( fid, '%d];\n',shankind(nchan) + 1 );

        xcoords = shankind*shankSep + xcoords;

        fprintf( fid, 'siteLoc = [' );
        for i = 1:nchan-1
            fprintf(fid, '%d,%d;', xcoords(i), ycoords(i));
        end
        fprintf( fid, '%d,%d];\n', xcoords(nchan), ycoords(nchan) );

        fprintf( fid, 'siteMap = [' );
        for i = 1:nchan-1
            fprintf( fid, '%d,', chans(i) );
        end
        fprintf( fid, '%d];\n', chans(nchan) );
        fclose(fid);
end

end


function [meta] = ReadMeta(metaName)

% Parse ini file into cell entries C{1}{i} = C{2}{i}
fid = fopen(metaName, 'r');
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
end % ReadMeta


% =========================================================
% Return shank and electrode number for NP2.0.
%
% Index into these with original (acquired) channel IDs.
%
function [elecInd, shankInd, bankMask, connected] = NP20_ElecInd(meta)
pType = str2num(meta.imDatPrb_type);
if pType == 21
    % Single shank probe
    % imro table entries: (channel, bank, refType, electrode #)
    C = textscan(meta.imroTbl, '(%*s %d %*s %d', ...
        'EndOfLine', ')', 'HeaderLines', 1 );
    elecInd = int32(cell2mat(C(2)));
    bankMask = int32(cell2mat(C(1)));
    shankInd = zeros(size(elecInd));
    connected = ones(size(elecInd));
    exChan = findDisabled(meta);
    for i = 1:numel(exChan)
        connected(elecInd == exChan(i)) = 0;
    end

else
    % 4 shank probe
    % imro table entries: (channel, shank, bank, refType, electrode #)
    C = textscan(meta.imroTbl, '(%d %d %d %*s %d', ...
        'EndOfLine', ')', 'HeaderLines', 1 );
    chan = double(cell2mat(C(1)));
    elecInd = int32(cell2mat(C(4)));
    bankMask = int32(cell2mat(C(3)));
    shankInd = double(cell2mat(C(2)));
    connected = ones(size(chan));
    exChan = findDisabled(meta);
    %exChan = [127];
    for i = 1:numel(exChan)
        connected(chan == exChan(i)) = 0;
    end
end
end % NP20_ElecInd


% =========================================================
% Return shank and electrode number for NP1.0.
%
% Index into these with original (acquired) channel IDs.
%
function [elecInd, connected] = NP10_ElecInd(meta)

% 3A or 3B data?
% 3A metadata has field "typeEnabled" which was replaced
% with "typeImEnabled" and "typeNiEnabled" in 3B.
% The 3B imro table has an additional field for the
% high pass filter enabled/disabled
% Note that the textscan funtion places line breaks at each
% instance of the 'EndofLine' character -- here, ')'
% 'HeaderLines' = 1 skips the initial entry in the table with
% the probe type and number of entries.
if isfield(meta,'typeEnabled')
    % 3A data
    C = textscan(meta.imroTbl, '(%d %d %*s %*s %*s', ...
        'EndOfLine', ')', 'HeaderLines', 1 );
    exChan = findDisabled(meta);
    %exChan = [36, 75, 112, 151, 188, 227, 264, 303, 340, 373];
else
    % 3B data
    C = textscan(meta.imroTbl, '(%d %d %*s %*s %*s %*s', ...
        'EndOfLine', ')', 'HeaderLines', 1 );
    exChan = findDisabled(meta);
    %exChan = [191];

end
chan = double(cell2mat(C(1)));
bank = double(cell2mat(C(2)));
elecInd = bank*384 + chan;
connected = ones(size(chan));
for i = 1:numel(exChan)
    connected(chan == exChan(i)) = 0;
end

end % NP10_ElecInd

% =========================================================
% Read shank map for any probe type and return list
% of channels that are disabled. This will include the
% reference channels
%
% Note that the textscan funtion places line breaks at each
% instance of the 'EndofLine' character -- here, ')'
% 'HeaderLines' = 1 skips the initial entry in the table with
% the number of shanks, columns, and rows
function [exChan] = findDisabled(meta)
% read in the shank map
C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
    'EndOfLine', ')', 'HeaderLines', 1 );
enabled = double(cell2mat(C(4)));
% There's an entry in the shank map for each saved channel.
% Get the array of saved channels:
chan = OriginalChans(meta);
% Find out how many are non-SY chans
[AP,~,~] = ChannelCountsIM(meta);
exChan = [];
for i = 1:AP
    if enabled(i) == 0
        exChan = [exChan, chan(i)];
    end
end
end % findDisabled

% =========================================================
% Return x y coords for electrode index for 2.0 probes
%
%
function [xCoord, yCoord] = XYCoord20(meta, elecInd, bankMask, shankind)

pType = str2num(meta.imDatPrb_type);

nElec = 1280;   %per shank; pattern repeats for the four shanks
vSep = 15;   % in um
hSep = 32;

elecPos = zeros(nElec, 2);

elecPos(1:2:end,1) = 0;         %sites 0,2,4...
elecPos(2:2:end,1) = hSep;      %sites 1,3,5...

% fill in y values
viHalf = (0:(nElec/2-1))';                %row numbers
elecPos(1:2:end,2) = viHalf * vSep;       %sites 0,2,4...
elecPos(2:2:end,2) = elecPos(1:2:end,2);  %sites 1,3,5...


xCoord = elecPos(elecInd+1,1);
yCoord = elecPos(elecInd+1,2);

if pType == 21
    % single shank probe. Plot only lowest selected electrode
    figure(1)
    % plot all positions
    scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
    scatter( xCoord, yCoord, 100, 'b', 'square', 'filled' );hold on;
    xlim([-16,64]);
    ylim([-10,10000]);
    title('NP 2.0 single shank view');
    hold off;
else
    % four shank probe, no multiple connections
    figure(1)
    shankSep = 250;
    for sI = 0:3
        cc = find(shankind == sI);
        scatter( shankSep*sI + elecPos(:,1), elecPos(:,2), 30, 'k', 'square' ); hold on;
        scatter( shankSep*sI + xCoord(cc), yCoord(cc), 20, 'b', 'square', 'filled' ); hold on;
    end
    xlim([-16,3*shankSep+64]);
    ylim([-10,10000]);
    title('NP2.0 MS shank view');
    hold off;
end


end % XY20Coord

% =========================================================
% Return x y coords for electrode index for 1.0 probes
%
%
function [xCoord, yCoord] = XYCoord10(meta, elecInd)

nElec = 960;   %per shank; pattern repeats for the four shanks
vSep = 20;   % in um
hSep = 32;

elecPos = zeros(nElec, 2);

elecPos(1:4:end,1) = hSep/2;            %sites 0,4,8...
elecPos(2:4:end,1) =  (3/2)*hSep;       %sites 1,5,9...
elecPos(3:4:end,1) = 0;                 %sites 2,6,10...
elecPos(4:4:end,1) =  hSep;             %sites 3,7,11...
elecPos(:,1) = elecPos(:,1) + 11;       %x offset on the shank

% fill in y values
viHalf = (0:(nElec/2-1))';                %row numbers
elecPos(1:2:end,2) = viHalf * vSep;       %sites 0,2,4...
elecPos(2:2:end,2) = elecPos(1:2:end,2);  %sites 1,3,5...

xCoord = elecPos(elecInd+1,1);
yCoord = elecPos(elecInd+1,2);

% single shank probe. Plot only lowest selected electrode
figure(1)
% plot all positions
scatter( elecPos(:,1), elecPos(:,2), 150, 'k', 'square' ); hold on;
scatter( xCoord, yCoord, 100, 'b', 'square', 'filled' );hold on;
xlim([0,70]);
ylim([-10,8000]);
title('NP 1.0 single shank view');
hold off;

end % XY10Coord


% =========================================================
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
%
function chans = OriginalChans(meta)
if strcmp(meta.snsSaveChanSubset, 'all')
    chans = (1:str2double(meta.nSavedChans));
else
    chans = str2num(meta.snsSaveChanSubset);
    chans = chans + 1;
end
end % OriginalChans


% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
function [AP,LF,SY] = ChannelCountsIM(meta)
M = str2num(meta.snsApLfSy);
AP = M(1);
LF = M(2);
SY = M(3);
end % ChannelCountsIM


% =========================================================
% Parse snsGeomMap for XY coordinates
%
function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = geomMapToGeom(meta)

C = textscan(meta.snsGeomMap, '(%d:%d:%d:%d', ...
    'EndOfLine', ')', 'HeaderLines', 1 );
shankInd = double(cell2mat(C(1)));
xCoord = double(cell2mat(C(2)));
yCoord = double(cell2mat(C(3)));
connected = double(cell2mat(C(4)));

% parse header for number of shanks
geomStr = meta.snsGeomMap;
headStr = extractBefore(geomStr,')(');
headParts = split(headStr,',');
nShank = str2double(headParts{2});
shankWidth = str2double(headParts{4});
shankPitch = str2double(headParts{3});
end % geomMapToGeom


% =========================================================
% Get XY coordinates from snsShankMap plus hard coded geom values
%
function [nShank, shankWidth, shankPitch, shankInd, xCoord, yCoord, connected] = shankMapToGeom(meta)
% get number of saved AP channels (some early metadata files have a
% SYNC entry in the snsChanMap
[nchan,~,~] = ChannelCountsIM(meta);

C = textscan(meta.snsShankMap, '(%d:%d:%d:%d', ...
    'EndOfLine', ')', 'HeaderLines', 1 );
shankInd = double(cell2mat(C(1)));
colInd = double(cell2mat(C(2)));
rowInd = double(cell2mat(C(3)));
connected = double(cell2mat(C(4)));

% trim these to the number of saved channels
shankInd = shankInd(1:nchan);
colInd = colInd(1:nchan);
rowInd = rowInd(1:nchan);
connected = connected(1:nchan);

geom = getGeomParams(meta);

oddRows = logical(mod(rowInd,2));
evenRows = ~oddRows;
xCoord = colInd*geom.horzPitch;
xCoord(evenRows) = xCoord(evenRows) + geom.even_xOff ;
xCoord(oddRows) = xCoord(oddRows) + geom.odd_xOff;
yCoord = rowInd*geom.vertPitch;

nShank = geom.nShank;
shankWidth = geom.shankWidth;
shankPitch = geom.shankPitch;
end % shankMapToGeom


% =========================================================
% Return geometry paramters for supported probe types
% These are used to calculate positions from metadata
% that includes only ~snsShankMap
%
function geom = getGeomParams(meta)
% create map
geomTypeMap = makeTypeMap();

% get probe part number; if absent, this is a 3A
if isfield(meta,'imDatPrb_pn')
    pn = meta.imDatPrb_pn;
else
    pn = '3A';
end

if geomTypeMap.isKey(pn)
    geomType = geomTypeMap(pn);
else
    fprintf('unsupported probe part number\n');
    return;
end

switch geomType
    case 'np1_stag_70um'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'nhp_lin_70um'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'nhp_stag_125um_med'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 87;
        geom.vertPitch = 20;
        geom.rowsPerShank = 1368;
        geom.elecPerShank = 2496;
    case 'nhp_stag_125um_long'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 87;
        geom.vertPitch = 20;
        geom.rowsPerShank = 2208;
        geom.elecPerShank = 4416;
    case 'nhp_lin_125um_med'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 103;
        geom.vertPitch = 20;
        geom.rowsPerShank = 1368;
        geom.elecPerShank = 2496;
    case 'nhp_lin_125um_long'
        geom.nShank = 1;
        geom.shankWidth = 125;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 103;
        geom.vertPitch = 20;
        geom.rowsPerShank = 2208;
        geom.elecPerShank = 4416;
    case 'uhd_8col_1bank'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 14;
        geom.odd_xOff = 14;
        geom.horzPitch = 6;
        geom.vertPitch = 6;
        geom.rowsPerShank = 48;
        geom.elecPerShank = 384;
   case 'uhd_8col_16bank'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 14;
        geom.odd_xOff = 14;
        geom.horzPitch = 6;
        geom.vertPitch = 6;
        geom.rowsPerShank = 768;
        geom.elecPerShank = 6144;
    case 'np2_ss'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 15;
        geom.rowsPerShank = 640;
        geom.elecPerShank = 1280;
     case 'np2_4s'
        geom.nShank = 4;
        geom.shankWidth = 70;
        geom.shankPitch = 250;
        geom.even_xOff = 27;
        geom.odd_xOff = 27;
        geom.horzPitch = 32;
        geom.vertPitch = 15;
        geom.rowsPerShank = 640;
        geom.elecPerShank = 1280;
    case 'NP1120'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 6.75;
        geom.odd_xOff = 6.75;
        geom.horzPitch = 4.5;
        geom.vertPitch = 4.5;
        geom.rowsPerShank = 192;
        geom.elecPerShank = 384;
    case 'NP1121'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 6.25;
        geom.odd_xOff = 6.25;
        geom.horzPitch = 3;
        geom.vertPitch = 3;
        geom.rowsPerShank = 384;
        geom.elecPerShank = 384;
    case 'NP1122'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 12.5;
        geom.odd_xOff = 12.5;
        geom.horzPitch = 3;
        geom.vertPitch = 3;
        geom.rowsPerShank = 24;
        geom.elecPerShank = 384;
    case 'NP1123'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 10.25;
        geom.odd_xOff = 10.25;
        geom.horzPitch = 4.5;
        geom.vertPitch = 4.5;
        geom.rowsPerShank = 32;
        geom.elecPerShank = 384;
    case 'NP1300'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 11;
        geom.odd_xOff = 11;
        geom.horzPitch = 48;
        geom.vertPitch = 20;
        geom.rowsPerShank = 480;
        geom.elecPerShank = 960;
    case 'NP1200'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 27;
        geom.odd_xOff = 11;
        geom.horzPitch = 32;
        geom.vertPitch = 20;
        geom.rowsPerShank = 64;
        geom.elecPerShank = 128;
    case 'NXT3000'
        geom.nShank = 1;
        geom.shankWidth = 70;
        geom.shankPitch = 0;
        geom.even_xOff = 53;
        geom.odd_xOff = 53;
        geom.horzPitch = 0;
        geom.vertPitch = 15;
        geom.rowsPerShank = 128;
        geom.elecPerShank = 128;
    otherwise
        % shouldn't see this case
        fprintf('unsupported probe part number\n');
        return;
end
end % getGeomParams




% =========================================================
% Return geometry paramters for supported probe types
% Note that geom only contains enough info to calculate
% positions for the electrodes listed in snsShankMap
%
function M = makeTypeMap()
% many part numbers have the same geometry parameters ;
% make a map that pairs geometry type (value) with probe part number (key)
M = containers.Map('KeyType','char','ValueType','char');

M('3A') = 'np1_stag_70um';
M('PRB_1_4_0480_1') = 'np1_stag_70um';
M('PRB_1_4_0480_1_C') = 'np1_stag_70um';
M('NP1010') = 'np1_stag_70um'; 
M('NP1011') = 'np1_stag_70um';
M('NP1012') = 'np1_stag_70um';
M('NP1013') = 'np1_stag_70um';

M('NP1015') = 'nhp_lin_70um';
M('NP1015') = 'nhp_lin_70um';
M('NP1016') = 'nhp_lin_70um';
M('NP1017') = 'nhp_lin_70um';
   
M('NP1020') = 'nhp_stag_125um_med';
M('NP1021') = 'nhp_stag_125um_med';
M('NP1030') = 'nhp_stag_125um_long';
M('NP1031') = 'nhp_stag_125um_long';

M('NP1022') = 'nhp_lin_125um_med';
M('NP1032') = 'nhp_lin_125um_long';

M('NP1100') = 'uhd_8col_1bank';
M('NP1110') = 'uhd_8col_16bank';

M('PRB2_1_2_0640_0') = 'np2_ss';
M('PRB2_1_4_0480_1') = 'np2_ss';
M('NP2000') = 'np2_ss';
M('NP2003') = 'np2_ss';
M('NP2004') = 'np2_ss';

M('PRB2_4_2_0640_0') = 'np2_4s';
M('PRB2_4_4_0480_1') = 'np2_4s';
M('NP2010') = 'np2_4s';
M('NP2013') = 'np2_4s';
M('NP2014') = 'np2_4s';

M('NP1120') = 'NP1120';
M('NP1121') = 'NP1121';
M('NP1122') = 'NP1122';
M('NP1123') = 'NP1123';
M('NP1300') = 'NP1300';

M('NP1200') = 'NP1200';
M('NXT3000') = 'NXT3000';
end % makeTypeMap
