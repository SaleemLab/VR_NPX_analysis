% Function to load spikes from a kilosort file
%
% History
% Tomaso Muzzu - 21/06/2019 - adapted from Jonas code
% SGS 21st April 2020 Adjusted from get_ks_Spikes to make uigetfile_n_dir
% call non-compulsory
function SpikeData = get_ks_Spikes_new(directory)




% load spikes from kilosort. Only loads good and mua.
DIRname = uigetfile_n_dir(directory,'Select the folder containing the clustered recording of interest');
sp = loadKSdir(DIRname{1});

% filename_clusters = fullfile(DIRname,'cluster_groups.csv');
%[cids, cgs] = readClusterGroupsCSV(filename_clusters);
% 0 = noise
% 1 = mua
% 2 = good
% 3 = unsorted
% noise_list = cids(cgs==0);

noise_list = [];
clear SpikeInfo
SpikeInfo{1,1} = 'Spiketimes';
SpikeInfo{2,1} = 'Unit type';
SpikeInfo{3,1} = 'Shank';
SpikeInfo{4,1} = 'Channel';
SpikeInfo{5,1} = 'ID';
SpikeInfo{6,1} = 'Frame rate';

for icell = 1:numel(sp.cids)
    
    SpikeInfo{1,icell+1} = sp.st(sp.clu==sp.cids(icell));
    SpikeInfo{2,icell+1} = sp.cgs(icell);
    SpikeInfo{3,icell+1} = 99;
    SpikeInfo{4,icell+1} = 99;
    SpikeInfo{5,icell+1} = sp.cids(icell);
    SpikeInfo{6,icell+1} = sp.sample_rate; 
    
end %for icell

% load metadata info
MetaDataFilePath = dir([DIRname{1} filesep 'M*dat_meta.mat']);
MetaData = load([DIRname{1} filesep MetaDataFilePath.name]);

% save SpikeData
SpikeData.SpikeInfo = SpikeInfo;
SpikeData.MetaData = MetaData;


end %function