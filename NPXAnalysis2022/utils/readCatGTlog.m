function catGT_table = readCatGTlog(filename)

%% Edward Horrocks 2021
% This function reads the CatGT.log file produced from concatenating spikeGLX
% recordings using CatGT and produces a matrix with the required
% information to find the start points of each recording

% Output: 
% Matrix of size (nRecordings, 4) where each row has the format
% [g_index, out_start_smp, inp_gap_smp, out_zeros_smp]
% The values have the following meaning:
% StimFirstSample(g_idx) = out_start_smp(g_idx)+out_zeros_smp(g_idx)
% StimLastSample(g_idx) = out_start_smp(g_idx+1);
%
% filename = 'CatGT_20210817.log';
%
% [Thd 140052433868608 CPU 0 4/28/22 17:32:37.186] Cmdline: CatGT -dir=/home/edd/temp_spikesorting/M22009/20220414 -run=M22009_20220414 -g=0,4 -t=0 -prb_fld -prb=0 -ap -zerofillmax=500 -aphipass=300 -gbldmx -gfix=0,0.1,0.02 -SY=0,-1,6,500
% [Thd 140052433868608 CPU 0 4/28/22 17:35:39.867] Gap before file 'M22009_20220414_g1_t0.imec0.ap.bin' out_start_smp=10116430 inp_gap_smp=3764962 out_zeros_smp=15000
% [Thd 140052433868608 CPU 0 4/28/22 17:42:05.231] Gap before file 'M22009_20220414_g2_t0.imec0.ap.bin' out_start_smp=32436698 inp_gap_smp=5887428 out_zeros_smp=15000
% [Thd 140052433868608 CPU 0 4/28/22 17:53:46.924] Gap before file 'M22009_20220414_g3_t0.imec0.ap.bin' out_start_smp=73167201 inp_gap_smp=14226813 out_zeros_smp=15000
% [Thd 140052433868608 CPU 0 4/28/22 18:03:51.475] Gap before file 'M22009_20220414_g4_t0.imec0.ap.bin' out_start_smp=107585600 inp_gap_smp=24829881 out_zeros_smp=15000
% [Thd 140052433868608 CPU 0 4/28/22 18:06:12.241] Run M22009_20220414_g0 Gfix prb 0 edits/sec 0.041
% [Thd 140052433868608 CPU 0 4/28/22 18:06:13.395] 


%% load data using textscan
delimiter = {''};
    startRow = 1;
    endRow = inf;
formatSpec = '%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end
fclose(fileID);
logText = dataArray{:, 1};

%% get relevant lines of text file
% filter .log file to get only the lines with required info.
text2find = 'Gap before file';
rel_lines = [];

for irow = 1:size(logText,1)
    k = strfind(logText(irow,:),text2find);
    if ~isempty(k)
        rel_lines(end+1) = irow;
    end
end

filt_text = logText(rel_lines);

%% get required values
% each line repsents a gated recoring and has format: 
% [g_index, out_start_smp, inp_gap_smp, out_zeros_smp]

logVals = [];

% get gVals
for irow = 1:size(filt_text)
    str = char(filt_text(irow,:));
    
    % get gX idx (based
    temp_val = extractBetween(str,'_g','_t');
    logVals(irow,1) = str2num(temp_val{1});
    
    str(strfind(str, '=')) = [];
    Key = 'out_start_smp';
    Index = strfind(str, Key);
    logVals(irow,2) = sscanf(str(Index(1) + length(Key):end), '%g', 1);
    Key = 'inp_gap_smp';
    Index = strfind(str, Key);
    logVals(irow,3) = sscanf(str(Index(1) + length(Key):end), '%g', 1);
    Key = 'out_zeros_smp';
    Index = strfind(str, Key);
    logVals(irow,4) = sscanf(str(Index(1) + length(Key):end), '%g', 1);
    
    
end

catGT_table = array2table(logVals, 'VariableNames', {'g_idx', 'out_start_smp', 'inp_gap_smp', 'out_zeros_smp'});

end

