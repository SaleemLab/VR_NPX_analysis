function [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column)

% load chan_config and sorted config (from one of the four columns) based on ImecMeta files

% Find config files for ephys data information
% [AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);
% imecMeta = ReadMeta(fullfile(options.EPHYS_DATAPATH,LF_FILE));
[AP_FILE,LF_FILE] = findImecMetaFile(options.ANALYSIS_DATAPATH,options.probe_id);

switch options.importMode
    case 'LF'
        file_to_use = LF_FILE;
    case 'KS'
        file_to_use = AP_FILE;
end

file_to_use = strrep(file_to_use,'meta','bin'); % replace meta to bin.

imecMeta = NPadmin.ReadMeta(file_to_use,options.ANALYSIS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);
% switch str2double(imecMeta.imDatPrb_type)
%     case 0 
%     electrode_spacing_um = 20;
%     otherwise
%     electrode_spacing_um = 15;
% end

sorted_config = [];

% Get channels from one column (x coordinate = 11, 27, 43 or 59 micron)
if ~isempty(column)
    % Set probe columns (4 columns for NPX1)
    % (x coordinate = 11, 27, 43 or 59 micron)
    cols_available = [];

    shank_id = []; % the corresponding shank_id
    shanks_available = unique(chan_config.Shank);
    for n=1:size(shanks_available)
        cols_available = [cols_available unique(chan_config.Ks_xcoord(chan_config.Shank==shanks_available(n)))];
        shank_id = [shank_id shanks_available(n)];
    end

    col_ID = cols_available(column); % (1) is 11
    sorted_config = sortrows(chan_config,'Ks_ycoord','descend');
    col_idx = sorted_config.Ks_xcoord == col_ID;
    sorted_config = sorted_config(col_idx,:);
end
end