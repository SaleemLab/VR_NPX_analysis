%% Loading extracted SO ripples spindles info for plotting and analysis 
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\OneDrive\Documents\GitHub\VR_NPX_analysis'))


if exist('D:\corticohippocampal_replay')>0
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

load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_whole.mat'));

load(fullfile(analysis_folder,'KDE_reactivation_DOWN_all_POST.mat'),'KDE_reactivation_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_UP_all_POST.mat'),'KDE_reactivation_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_all_POST.mat'),'KDE_reactivation_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_DOWN_all_POST.mat'),'KDE_reactivation_V1_DOWN_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_UP_all_POST.mat'),'KDE_reactivation_V1_UP_all')
load(fullfile(analysis_folder,'KDE_reactivation_V1_all_POST.mat'),'KDE_reactivation_V1_all')


all_sessions = max(slow_waves_all(1).DOWN_session_count);
sessions_to_process = 1:all_sessions;
for nsession = 1:all_sessions
    
    for nprobe = 1:2
        event_index = intersect(find(slow_waves_all(nprobe).UP_session_count==nsession),probability(nprobe).UP_all_index);
        for nevent = 1:length(event_index)


        end

        KDE_reactivation_V1_UP_all(nprobe).event_id{nsession}

        probability(nprobe).DOWN_all_index

        
        probability(nprobe).UP_all_index
    end

end


% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov_normalised.mat'));
% load(fullfile(analysis_folder,'V1-HPC sleep interaction','SO_ripples_probability_markov.mat'));