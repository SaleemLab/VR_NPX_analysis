% Find g folder
function folderName = findGFolder(fpath,targetG)
% Get the folders in the ephys folder
tb = dir(fpath);
% Find the g folders in that
tc = {};
[tc{1:length(tb)}] = deal(tb.name);
td = tc(contains(tc,'_g'));
% Get the G nums out of the folder name
te = regexpi(td,'_g(\w*)','tokens');
tf = str2double(cellfun(@(x) x{1}{1}, te,'UniformOutput',false));
% Match to target
tp = find(tf == targetG);
if ~isempty(tp)
    folderName = td{tp};
else
    error('No match found')
end
