function [clusters chan_config sorted_config] = select_clusters_NPX(options,varargin)

%%% Function to extract clusters based on different sorters
% Was modified/simplified from old function: 
% [clusters chan_config sorted_config] = load_KS_NPX1(options,column,varargin)
% Modified from original version of load_KS_NPX1 function


% Extract NPX channel configuration
options.importMode = 'KS';
[file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

% Default values
p = inputParser;
addParameter(p,'selected_channels',[1 size(chan_config,1)],@isnumeric) % Select channels for analysis (default is all the channles)
addParameter(p,'group','good clusters',@isstr) % spike data grouped 'good clusters'
addParameter(p,'plot_option',0,@isnumeric) % option to plot
addParameter(p,'SR',60,@isnumeric) % SR for cluster spike count time course (normally 60Hz)
addParameter(p,'tvec',[],@isnumeric) % Time vector

% addParameter(p,'cell_exporer','off',@isstr) % cell explorer option
% removed (replaced by few lines of code that can classify the cell based on )
addParameter(p,'waveform','off',@isstr) % extract individual spike waveform for each cluster (not mean waveform)
addParameter(p,'sorter','off',@isstr) % What sorters to use. If off, it is the default kilosort 3.
% If 'KS3' or 'KS2', it is spike interface KS3 and KS2.


% addParameter(p,'stdev',[],@isnumeric)
% addParameter(p,'show','on',@isstr)
% addParameter(p,'noise',[],@ismatrix)
% addParameter(p,'passband',[125 300],@isnumeric)
% addParameter(p,'EMGThresh',.9,@isnumeric);
% addParameter(p,'saveMat',false,@islogical);
% addParameter(p,'minDuration',20,@isnumeric)
% addParameter(p,'plotType',1,@isnumeric)
% addParameter(p,'savepath',[],@isstr)

% assign parameters (either defaults or given)
parse(p,varargin{:});
selected_channels = p.Results.selected_channels;
group = p.Results.group;
tvec = p.Results.tvec;
SR = p.Results.SR;
sorter = p.Results.sorter;

amplitude_cutoff = cluster_metrics.amplitude_cutoff; % <0.1 (just in case there are cells that are selectively active)
isi_violation = cluster_metrics.isi_viol; % < 0.5
presence_ratio = cluster_metrics.presence_ratio; % > 0.5 %Should be at least active more than 50% of the time
troughToPeak = cluster_metrics.duration; % Duration (in seconds) between peak and trough
amplitude = cluster_metrics.amplitude;

good_unit = cluster_metrics.amplitude_cutoff <= 0.1...
    &cluster_metrics.isi_viol <= 0.1...
    &cluster_metrics.amplitude >=50;

end