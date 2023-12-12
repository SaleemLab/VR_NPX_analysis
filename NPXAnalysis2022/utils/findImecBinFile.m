% Programme to return the full path of bin files in target folder
% SFS 14th April 2022
function [ap_file,lf_file] = findImecBinFile(EPHYS_DATAPATH)
td = dir(fullfile(EPHYS_DATAPATH,'*.bin'));
fnames = {};
[fnames{1:length(td)}] = deal(td.name);
ap_file = find(contains(fnames,'ap.bin'));
if ~isempty(ap_file)
    ap_file = fnames{ap_file};
else ap_file = [];
end
lf_file = find(contains(fnames,'lf.bin'));
if ~isempty(lf_file)
    lf_file = fnames{lf_file};
else lf_file = [];
end
