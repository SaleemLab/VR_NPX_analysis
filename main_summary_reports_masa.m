if ismac
    %%merge all sessions
    addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))
    % Load clusters_all_all from base_folder
    base_folder = fullfile('/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data');
else
    addpath(genpath('C:\Users\adam.tong\Documents\GitHub\V1_MEC_acute_MAT'))
    base_folder = fullfile('C:\Users\adam.tong\OneDrive - University College London\data');
end

load(fullfile(base_folder, 'clusters_all_ks3.mat'));
load(fullfile(base_folder,'spatial_responses_V1_HVA_ks3.mat'))

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
V1_index = (clusters_all.region == "V1");
HVA_index =  clusters_all.region =='HVA';
spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = (V1_index| HVA_index) & spatially_tuned_neurons;
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);

session_count = 0;
save_path = fullfile(base_folder,'summary_figures','population_summary_V1');
mkdir(save_path)
for iS = 5:19
    session_count= session_count + 1;
    clusters_this_session = clusters_all.session_count(overall_cluster_index) == iS;
    cluster_id_this_session = cluster_id(clusters_this_session);
    first_cluster_id = num2str(cluster_id_this_session(1));
    session_name = ['M',first_cluster_id(1:7)];
    ids_index = zeros(size(cluster_id));
    ids_index(clusters_this_session) = 1;
    ids_index = logical(ids_index);
    cluster_summary(clusters_all,cluster_id(ids_index),...
        'SMI',spatial_response(ids_index,:), ...
        spatial_response_extended(ids_index,:),...
        spatial_response_all(ids_index,:), ...
        spatial_response_extended_all(ids_index,:),'plot_population',1,'break_loop',1,'plot_cluster',0);
    saveas(gcf,fullfile(save_path,[session_name,'.jpg']),'jpg');
    close all
end
save_path = fullfile(base_folder,'summary_figures','per_neuron_report_V1');
mkdir(save_path)
session_count = 0;
for iS = 1:19
    session_count= session_count + 1;
    clusters_this_session = clusters_all.session_count(overall_cluster_index) == iS;
    cluster_id_this_session = cluster_id(clusters_this_session);
    first_cluster_id = num2str(cluster_id_this_session(1));
    session_name = ['M',first_cluster_id(1:7)];
    ids_index_this_session = zeros(size(cluster_id));
    ids_index_this_session(clusters_this_session) = 1;
    ids_index_this_session = logical(ids_index_this_session);

    mkdir(fullfile(save_path,session_name));
    for iC = 1:length(cluster_id_this_session)

        ids_index = zeros(size(cluster_id));
        ids_index(ismember(cluster_id,cluster_id_this_session(iC))) = 1;
        ids_index = logical(ids_index);
        cluster_summary(clusters_all,cluster_id(ids_index),...
            'SMI',spatial_response(ids_index,:), ...
            spatial_response_extended(ids_index,:),...
            spatial_response_all(ids_index,:), ...
            spatial_response_extended_all(ids_index,:),'plot_population',0,'break_loop',1,'plot_cluster',1);
        saveas(gcf,fullfile(save_path,session_name,[num2str(cluster_id_this_session(iC)),'.jpg']),'jpg');
        close all


    end
end


%% MEC clusters

load(fullfile(base_folder, 'clusters_all_ks3.mat'));
load(fullfile(base_folder,'spatial_responses_MEC_ks3.mat'))

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
MEC_index = (clusters_all.region == "MEC");


spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = MEC_index & spatially_tuned_neurons;
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);

session_count = 0;
save_path = fullfile(base_folder,'summary_figures','population_summary_MEC');
mkdir(save_path)
for iS = 5:19
    session_count= session_count + 1;
    clusters_this_session = clusters_all.session_count(overall_cluster_index) == iS;
    cluster_id_this_session = cluster_id(clusters_this_session);
    first_cluster_id = num2str(cluster_id_this_session(1));
    session_name = ['M',first_cluster_id(1:7)];
    ids_index = zeros(size(cluster_id));
    ids_index(clusters_this_session) = 1;
    ids_index = logical(ids_index);
    cluster_summary(clusters_all,cluster_id(ids_index),...
        'MEC',spatial_response(ids_index,:), ...
        spatial_response_extended(ids_index,:),...
        spatial_response_all(ids_index,:), ...
        spatial_response_extended_all(ids_index,:),'plot_population',1,'break_loop',1,'plot_cluster',0);
    saveas(gcf,fullfile(save_path,[session_name,'.jpg']),'jpg');
    close all
end
save_path = fullfile(base_folder,'summary_figures','per_neuron_report_MEC');
mkdir(save_path)
session_count = 0;
for iS = 5:19
    session_count= session_count + 1;
    clusters_this_session = clusters_all.session_count(overall_cluster_index) == iS;
    cluster_id_this_session = cluster_id(clusters_this_session);
    first_cluster_id = num2str(cluster_id_this_session(1));
    session_name = ['M',first_cluster_id(1:7)];
    ids_index_this_session = zeros(size(cluster_id));
    ids_index_this_session(clusters_this_session) = 1;
    ids_index_this_session = logical(ids_index_this_session);

    mkdir(fullfile(save_path,session_name));
    for iC = 1:length(cluster_id_this_session)

        ids_index = zeros(size(cluster_id));
        ids_index(ismember(cluster_id,cluster_id_this_session(iC))) = 1;
        ids_index = logical(ids_index);
        cluster_summary(clusters_all,cluster_id(ids_index),...
            'MEC',spatial_response(ids_index,:), ...
            spatial_response_extended(ids_index,:),...
            spatial_response_all(ids_index,:), ...
            spatial_response_extended_all(ids_index,:),'plot_population',0,'break_loop',1,'plot_cluster',1);
        saveas(gcf,fullfile(save_path,session_name,[num2str(cluster_id_this_session(iC)),'.jpg']),'jpg');
        close all


    end
end

%% MEC clusters

load(fullfile(base_folder, 'clusters_all_ks3.mat'));
load(fullfile(base_folder,'spatial_responses_HPC_ks3.mat'))

% select clusters you want to plot for population analysis e.g. only plot neurons from V1
HPC_index = (clusters_all.region == "HPC");


spatially_tuned_neurons = (clusters_all.odd_even_stability(:,1) >=0.95 ...
    |clusters_all.odd_even_stability(:,2) >= 0.95) ...
    & (clusters_all.peak_percentile(:,1) >=0.95 | clusters_all.peak_percentile(:,2) >=0.95)...
    &(clusters_all.first_second_stability(:,1)>0.95 | clusters_all.first_second_stability(:,2)>0.95)...
    &(clusters_all.reliability(:,1)>0.005 | clusters_all.reliability(:,2)>0.005);
overall_cluster_index = HPC_index & spatially_tuned_neurons;
cluster_id = clusters_all.cluster_id(overall_cluster_index);
% Call population_spatial_summary function with cluster_id
% population_spatial_summary(clusters_all, cluster_id, spatial_response,spatial_response_extended);

session_count = 0;
save_path = fullfile(base_folder,'summary_figures','population_summary_HPC');
mkdir(save_path)
for iS = 5:19
    session_count= session_count + 1;
    clusters_this_session = clusters_all.session_count(overall_cluster_index) == iS;
    cluster_id_this_session = cluster_id(clusters_this_session);
    first_cluster_id = num2str(cluster_id_this_session(1));
    session_name = ['M',first_cluster_id(1:7)];
    ids_index = zeros(size(cluster_id));
    ids_index(clusters_this_session) = 1;
    ids_index = logical(ids_index);
    cluster_summary(clusters_all,cluster_id(ids_index),...
        'MEC',spatial_response(ids_index,:), ...
        spatial_response_extended(ids_index,:),...
        spatial_response_all(ids_index,:), ...
        spatial_response_extended_all(ids_index,:),'plot_population',1,'break_loop',1,'plot_cluster',0);
    saveas(gcf,fullfile(save_path,[session_name,'.jpg']),'jpg');
    close all
end
save_path = fullfile(base_folder,'summary_figures','per_neuron_report_HPC');
mkdir(save_path)
session_count = 0;
for iS = 5:19
    session_count= session_count + 1;
    clusters_this_session = clusters_all.session_count(overall_cluster_index) == iS;
    cluster_id_this_session = cluster_id(clusters_this_session);
    first_cluster_id = num2str(cluster_id_this_session(1));
    session_name = ['M',first_cluster_id(1:7)];
    ids_index_this_session = zeros(size(cluster_id));
    ids_index_this_session(clusters_this_session) = 1;
    ids_index_this_session = logical(ids_index_this_session);

    mkdir(fullfile(save_path,session_name));
    for iC = 1:length(cluster_id_this_session)

        ids_index = zeros(size(cluster_id));
        ids_index(ismember(cluster_id,cluster_id_this_session(iC))) = 1;
        ids_index = logical(ids_index);
        cluster_summary(clusters_all,cluster_id(ids_index),...
            'MEC',spatial_response(ids_index,:), ...
            spatial_response_extended(ids_index,:),...
            spatial_response_all(ids_index,:), ...
            spatial_response_extended_all(ids_index,:),'plot_population',0,'break_loop',1,'plot_cluster',1);
        saveas(gcf,fullfile(save_path,session_name,[num2str(cluster_id_this_session(iC)),'.jpg']),'jpg');
        close all


    end
end

