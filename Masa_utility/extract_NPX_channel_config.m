function [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column)

% load chan_config and sorted config (from one of the four columns) based on ImecMeta files

% Find config files for ephys data information
[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
% imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,LF_FILE));

switch options.importMode
    case 'LF'
        file_to_use = LF_FILE;
    case 'KS'
        file_to_use = AP_FILE;
end

imecMeta = NPadmin.ReadMeta(file_to_use,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);
switch str2double(imecMeta.imDatPrb_type)
    case 0 
    electrode_spacing_um = 20;
    otherwise
    electrode_spacing_um = 15;
end

% Set probe columns (4 columns for NPX1)
% (x coordinate = 11, 27, 43 or 59 micron)
shanks_available = unique(chan_config.Shank);
for n=1:size(shanks_available)
    cols_available = unique(chan_config.Ks_xcoord(chan_config.Shank==shanks_available(n)));
end

% Get channels from one column (x coordinate = 11, 27, 43 or 59 micron)
col_ID = cols_available(column); % (1) is 11
sorted_config = sortrows(chan_config,'Ks_ycoord','descend');
col_idx = sorted_config.Ks_xcoord == col_ID;
sorted_config = sorted_config(col_idx,:);

end