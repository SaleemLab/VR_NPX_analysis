% Programme to return the full path of bin files in target folder
% Masa 21/02/2024
function [ap_file,lf_file] = findImecMetaFile(EPHYS_DATAPATH,probe_id,gFileNum)
td = dir(fullfile(EPHYS_DATAPATH,sprintf('*g%i*imec%i*.meta',gFileNum,probe_id)));
fnames = {};
[fnames{1:length(td)}] = deal(td.name);
ap_file = find(contains(fnames,'ap.meta'));
if ~isempty(ap_file)
    ap_file = fnames{ap_file};
else ap_file = [];
end
lf_file = find(contains(fnames,'lf.meta'));
if ~isempty(lf_file)
    lf_file = fnames{lf_file};
else lf_file = [];
end
