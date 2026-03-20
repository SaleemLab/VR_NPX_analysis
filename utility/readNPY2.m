function data = readNPY2(filename)
    fid = fopen(filename, 'rb');
    assert(fid > 0, 'Could not open file');
    
    % Read the magic string
    magic = fread(fid, 6, 'uint8=>char')';
    assert(strcmp(magic, char([147 'NUMPY'])), 'Not a valid NPY file');
    
    % Read version number
    version = fread(fid, 2, 'uint8');
    
    % Read header length
    if version(1) == 1
        headerLen = fread(fid, 1, 'uint16');
    elseif version(1) == 2
        headerLen = fread(fid, 1, 'uint32');
    else
        error('Unsupported NPY file version');
    end
    
    % Read header
    header = fread(fid, headerLen, 'uint8=>char')';
    header = strtrim(header); % Remove padding spaces
    
    % Extract dtype and shape from header
    dtypeMatch = regexp(header, '''descr'':\s*''([^'']+)''', 'tokens');
    shapeMatch = regexp(header, '''shape'':\s*\(([^)]+)\)', 'tokens');
    assert(~isempty(dtypeMatch), 'Data type not found in header');
    assert(~isempty(shapeMatch), 'Shape not found in header');
    
    dtype = dtypeMatch{1}{1};
    shape = str2num(strrep(shapeMatch{1}{1}, ',', ' ')); %#ok<ST2NM>
    
    % Map numpy dtype to MATLAB type
    switch dtype
        case '<f8', matlabType = 'double';
        case '<f4', matlabType = 'single';
        case '<i4', matlabType = 'int32';
        case '<i8', matlabType = 'int64';
        case '|u1', matlabType = 'uint8';
        otherwise, error('Unsupported dtype: %s', dtype);
    end
    
    % Read the data
    data = fread(fid, prod(shape), [matlabType '=>' matlabType]);
    data = reshape(data, shape);
    
    fclose(fid);
end