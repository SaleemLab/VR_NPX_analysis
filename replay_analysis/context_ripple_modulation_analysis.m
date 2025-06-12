%% Context-selective ripple modulation in V1 and HPC

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('C:\Users\masah\OneDrive\Documents\corticohippocampal_replay')
    analysis_folder = 'C:\Users\masah\OneDrive\Documents\corticohippocampal_replay';
elseif exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end
load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))
% load(fullfile(analysis_folder,'slow_waves_all_markov_POST.mat'))
load(fullfile(analysis_folder,'ripples_all_POST.mat'))
load(fullfile(analysis_folder,'spindles_all_POST.mat'))
load(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'))
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripples_event_info.mat'),'event_info');
load(fullfile(analysis_folder,'V1-HPC sleep interaction','merged_UP_DOWN_ripples_event_info.mat'),'merged_event_info');

load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA.mat'),'UP_DOWN_ripple_PSTH_MUA');
PSTH_MUA = UP_DOWN_ripple_PSTH_MUA;
load(fullfile(analysis_folder,'V1-HPC sleep interaction','UP_DOWN_ripple_PSTH_MUA_baseline.mat'));
PSTH_MUA_baseline = UP_DOWN_ripple_PSTH_MUA;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));
probability_psth_whole = probability;

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'bayesian_reactivation_V1_all_POST.mat'))
load(fullfile(analysis_folder,'bayesian_reactivation_all_POST.mat'))

sessions_to_process = 1:max(slow_waves_all(1).UP_session_count);


%% Context selective ripple modulation 

for nsession = 1:length(sessions_to_process)



end














%% Extract key information for DOWN UP, ripples and spindles
 
%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end

%% Extract key information for DOWN UP, ripples and spindles
 
%%%%%%% Find reference channel/shank
cortex_ref_shank = [];
HPC_ref_shank = [];

for nsession = 1:max(ripples_all(1).session_count)
    for probe_no = 1:2
        cortex_ref_shank(nsession,probe_no) = find(slow_waves_all(probe_no).shank_id{nsession} == slow_waves_all(probe_no).shank{nsession}(slow_waves_all(probe_no).channel{nsession} == slow_waves_all(probe_no).best_channel(nsession))...
            &slow_waves_all(probe_no).probe_hemisphere{nsession} == probe_no);
        % [~,idx] = min(abs(ripples_all(probe_no).SWR_peaktimes{nsession}' - ripples_all(probe_no).peaktimes(ripples_all(probe_no).session_count==nsession))');
        % ripple_counts = histcounts(idx,length(ripples_all(probe_no).shank_id{nsession}));
        % [~,HPC_ref_shank(nsession,probe_no)] = max(ripple_counts);

        shank_id = find(ripples_all(probe_no).probe_hemisphere{nsession} == probe_no);
        HPC_ref_shank(nsession,probe_no) = shank_id(ripples_all(probe_no).best_channel(nsession));

    end
end

%%%%%%%%%%%%%%%%%% Ripple info
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','ripples_spindles_probability_whole.mat'));
session_count = [ripples_all(1).session_count(ripples_all(1).SWS_index==1); ripples_all(2).session_count(ripples_all(2).SWS_index==1)];
subject_id = str2double(cellstr(ripples_all(1).subject(session_count,end-1:end)));
[~, ~, subject_id] = unique(subject_id);

%%% ripple power
ripple_info.ripple_power = [ripples_all(1).peak_zscore(ripples_all(1).SWS_index==1); ripples_all(2).peak_zscore(ripples_all(2).SWS_index==1)];

%%% spindle co-occurance
[~,spindle_index,~,index] =RestrictInts(merged_event_info.ripples_ints,merged_event_info.spindles_ints);
ripple_info.spindle_presence = spindle_index;
ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.spindle_presence_hemi(find(spindle_index)) = merged_event_info.spindles_hemisphere_id(index);

%%% spindle phase
spindle_phase=[];
for probe_no = 1:2
    spindle_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        spindle_phase{probe_no} = [spindle_phase{probe_no} ripples_all(probe_no).spindle_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

spindle_phase = [spindle_phase{1}(:,ripples_all(1).SWS_index==1) spindle_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.spindle_phase = spindle_phase';

%%% SO phase
SO_phase=[];
for probe_no = 1:2
    SO_phase{probe_no}=[];
    for nsession = 1:length(sessions_to_process)
        SO_phase{probe_no} = [SO_phase{probe_no} ripples_all(probe_no).SO_phase_ripple_peaktime{nsession}(cortex_ref_shank(nsession,:),:)];
    end
end

SO_phase = [SO_phase{1}(:,ripples_all(1).SWS_index==1) SO_phase{2}(:,ripples_all(2).SWS_index==1)];

ripple_info.SO_phase = SO_phase';

%%% early UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,1)+0.1]);
ripple_info.early_UP_index = UP_index;
ripple_info.early_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.early_UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);


%%% late UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1)-0.1 merged_event_info.DOWN_ints(:,1)]);
ripple_info.late_UP_index = UP_index;
ripple_info.late_UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.late_UP_index_hemi(find(UP_index)) = merged_event_info.DOWN_hemisphere_id(index);

%%% UP transition co-occurance
[~,UP_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.UP_ints(:,1) merged_event_info.UP_ints(:,2)]);
ripple_info.UP_index = UP_index;
ripple_info.UP_index_hemi = zeros(size(UP_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.UP_index_hemi(find(UP_index)) = merged_event_info.UP_hemisphere_id(index);

%%% DOWN transition co-occurance
[~,DOWN_index,~,index] = RestrictInts(merged_event_info.ripples_ints,[merged_event_info.DOWN_ints(:,1) merged_event_info.DOWN_ints(:,2)]);
ripple_info.DOWN_index = DOWN_index;
ripple_info.DOWN_index_hemi = zeros(size(DOWN_index));
% ripple_info.spindle_presence_hemi = zeros(size(spindle_index));
ripple_info.DOWN_index_hemi(find(DOWN_index)) = merged_event_info.DOWN_hemisphere_id(index);