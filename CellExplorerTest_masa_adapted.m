%%Cell Explorer On Kilosort Data
%add cell explorer to your path
SUBJECT = {'M23087'};
DATES = ['20231207';'20231208'];

SUBJECT = {'M23029'};
DATES = ['20230707'];

% Delete old cell explorer files (RUN if needed)
for nsubject = 1:length(SUBJECT)
    for iDate = 1:size(DATES,1)
        for nprobe = 1:2
            DATE = DATES(iDate,:);
            if isunix
                addpath(genpath('/research/USERS/Masa/CellExplorer'))
                basepath = ['/research/DATA/SUBJECTS/',SUBJECT{nsubject},'/ephys/',DATE,'/kilosort'];
            elseif ispc
%                             addpath(genpath('Z:\ibn-vision\USERS\Masa\code\CellExplorer'))
                 addpath(genpath('P:\corticohippocampal_replay\Masa\CellExplorer')) % Make sure only one cellexplorer folder is in the path
                
                if exist(['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort_probe_',num2str(nprobe)]) == 7
                    basepath = ['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort_probe_',num2str(nprobe)];
                elseif exist(['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort']) == 7
                    basepath = ['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort'];
                end
            end

            cd(basepath)
            delete([SUBJECT{nsubject},'_',DATE,'.spikes.cellinfo.mat'])
            delete([SUBJECT{nsubject},'_',DATE,'.session.mat'])
            delete([SUBJECT{nsubject},'_',DATE,'.cell_metrics.cellinfo.mat'])
            delete([SUBJECT{nsubject},'_',DATE,'.noiseLevel.channelInfo.mat'])

        end
    end
end


for nsubject = 1:length(SUBJECT)
    for iDate = 1:size(DATES,1)
        for nprobe = 1:2
            DATE = DATES(iDate,:);
            if isunix
                addpath(genpath('/research/USERS/Masa/CellExplorer'))
                basepath = ['/research/DATA/SUBJECTS/',SUBJECT{nsubject},'/ephys/',DATE,'/kilosort'];
            elseif ispc
%                             addpath(genpath('Z:\ibn-vision\USERS\Masa\code\CellExplorer'))
                 addpath(genpath('P:\corticohippocampal_replay\Masa\CellExplorer')) % Make sure only one cellexplorer folder is in the path
                
                if exist(['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort_probe_',num2str(nprobe)]) == 7
                    basepath = ['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort_probe_',num2str(nprobe)];
                elseif exist(['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort']) == 7
                    basepath = ['Z:\ibn-vision\DATA\SUBJECTS\',SUBJECT{nsubject},'\ephys\',DATE,'\kilosort'];
                end
            end

            cd(basepath)

            % M22069_20221201_g0_tcat.imec0.ap.bin       : raw data
            % rez.mat                   : metadata from KiloSort
            % PP02_2020-07-10.xml       : metadata from NeuroSuite ??
            % *.npy and *.tsv files 	: Spike data from Phy

            % Code for making 'all good' cluster_group.tsv for cell
            % explorer (it was put after validateSessionStruct() before 22/09/2023, meaning it was not correct)
            if ~exist('cluster_group_phy.tsv')
                phy_cluster_group = tdfread(fullfile(basepath,'cluster_group.tsv'));
                cluster_KSLabel = tdfread(fullfile(basepath,'cluster_KSLabel.tsv'));
                copyfile('cluster_group.tsv','cluster_group_phy.tsv')

                if length(cluster_KSLabel.cluster_id) == length(phy_cluster_group.cluster_id)
                    phy_cluster_group = tdfread(fullfile(basepath,'cluster_group_phy.tsv'));
                else % grab kslabel id
                    phy_cluster_group = [];
                    phy_cluster_group.cluster_id = cluster_KSLabel.cluster_id;
                    phy_cluster_group.group = cluster_KSLabel.KSLabel;
                end
                %                 phy_cluster_group = tdfread('cluster_group.tsv')
%                 count = 0;
                new_cluster_group = [];
                new_cluster_group.cluster_id = phy_cluster_group.cluster_id;
                for n = 1:length(phy_cluster_group.group)
                    new_cluster_group.group(n,:) = 'good';
                end

                cd(basepath)
                writetable(struct2table(new_cluster_group),'cluster_group.tsv','filetype','text', 'delimiter','\t')

            end

            % Run a template script to generate and import relevant session-level metadata
            session = sessionTemplate(basepath);

            session.general.name = [SUBJECT{nsubject},'_',DATE];
            session.general.date = DATE;
            session.general.sessionType = 'Acute';

            session.animal.name = SUBJECT{nsubject};
            session.animal.species = 'Mouse';
            session.animal.sex = 'Female';
            session.animal.strain = 'C57B1/6';

            session.extracellular.srLfp = 2500;
            %session.extracellular.fileName = ['/research/DATA/SUBJECTS/',SUBJECT,'/ephys/',DATE,'/kilosort/',SUBJECT,'_',DATE,'_g0_tcat.imec0.ap.bin'];
            session.extracellular.fileName = ['temp_wh.dat'];
            % session = gui_session(session);
            % You can validate that all required and optional fields for CellExplorer has been entered
            validateSessionStruct(session);


            % cluster_group = tdfread(fullfile(basepath,'cluster_group.tsv'));

            %% 2 Run the cell metrics pipeline 'ProcessCellMetrics' using the session struct as input
            tic
            cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'monoSynaptic_connections'},'showWaveforms',false,'sessionSummaryFigure',false);
            toc
            % Several files are generated here
            %
            % basename.cell_metrics.cellinfo.mat    : cell_metrics
            % basename.session.mat                  : session-level metadata
            % basename.spikes.cellinfo.mat          : spikes struct
            %
            % All files and structs are documented on the CellExplorer website:
            % https://cellexplorer.org/datastructure/data-structure-and-format/
            %
            % Once created the files can be loaded with dedicated scripts
            % session = loadSession;
            % cell_metrics = loadCellMetrics;
            % spikes = loadSpikes;

        end
    end
end
%% 3.1 Visualize the cell metrics in CellExplorer

cell_metrics = CellExplorer('metrics',cell_metrics);

%% get cell types (can be used without cell explorer)

for iunit = 1:numel(goodUnits)
    [ccg, t] = CCG(goodUnits(iunit).spike_times, ones(size(goodUnits(iunit).spike_times)),...
        'binSize', 0.0005, 'duration', 0.1,'norm', 'rate');
    goodUnits(iunit).acg = ccg;
    fit_params_out = fit_ACG(ccg,false);

    goodUnits(iunit).tau_rise = fit_params_out.acg_tau_rise;
end

narrow_idx = find([goodUnits.duration]<=0.45);
wide_idx = find([goodUnits.duration]>0.45 & [goodUnits.tau_rise]>6);
pyr_idx = find(~ismember(1:numel(goodUnits), [narrow_idx,wide_idx]));

[goodUnits.cellType] = deal(nan);
[goodUnits(pyr_idx).cellType] = deal(1);
[goodUnits(narrow_idx).cellType] = deal(2);
[goodUnits(wide_idx).cellType] = deal(3);

%% 3.2 Batch of sessions from Mice and rats (3686 cells from 111 sessions)

% load('/Volumes/Peter_SSD_4/cell_metrics/cell_metrics_peter_viktor.mat');
cd(basepath)
cd('Z:\ibn-vision\DATA\SUBJECTS\M23029\ephys\20230706\kilosort_probe_1')
load cluster_table
load M23029_20230706.cell_metrics.cellinfo.mat

cd('Z:\ibn-vision\DATA\SUBJECTS\M23028\ephys\20230706\analysis')
load best_channels

% cell_metrics = CellExplorer('metrics',cell_metrics);
cell_metrics.KS_cluster_id = double(cluster_table.cluster_id)';
cell_metrics.good_unit = double(cluster_table.cluster_id)';

cell_type = [];
SUA_count = 1;
V1_cell_type = [];
CA1_cell_type = [];
cell_type_all = [];

% Putatively CA1 from 190 to 220
% V1 from 234 to 334
for cell = 1:length(cell_metrics.putativeCellType)
    % Percentage of spikes with refractory period violation should be less
    % than 1% (Fraction of ISIs less than 2 ms)
    cell_metrics.good_unit = 0;
    cell_metrics.ks_good_unit = 0;
    if strcmp(cluster_table.label{cell},'good')
        cell_metrics.ks_good_unit = 1;
%     if cell_metrics.refractoryPeriodViolation(cell) < 0.5
        cell_metrics.good_unit = 1;

        if strcmp(cell_metrics.putativeCellType{cell},'Narrow Interneuron')
            cell_type_all(cell) = 1;
            cell_metrics.cell_type(cell) = 1;
        elseif strcmp(cell_metrics.putativeCellType{cell},'Wide Interneuron')
            cell_type_all(cell) = 2;
            cell_metrics.cell_type(cell) = 2;
        elseif strcmp(cell_metrics.putativeCellType{cell},'Pyramidal Cell')
            cell_type_all(cell) = 3;
            cell_metrics.cell_type(cell) = 3;
        else
            cell_type_all(cell) = 0;
            cell_metrics.cell_type(cell) = 0;
        end
        SUA_count = SUA_count + 1;

        if cell_metrics.maxWaveformCh(cell) < best_channels{1}.CA1_channel+10 & cell_metrics.maxWaveformCh(cell) > best_channels{1}.CA1_channel-10
            CA1_cell_type(cell) = cell_type_all(cell);
        else
            CA1_cell_type(cell) = 0;
        end

        if cell_metrics.maxWaveformCh(cell) < best_channels{1}.first_in_brain_channel & cell_metrics.maxWaveformCh(cell) > best_channels{1}.L5_channel-10
            V1_cell_type(cell) = cell_type_all(cell);
        else
            V1_cell_type(cell) = 0;
        end
    end
end
   

cell_type_all = cell_type_all
brain_region_cells = [cell_type_all; V1_cell_type; CA1_cell_type];
title_texts = {'All single units','V1 single units','CA1 single units'};
% sgtitle('All single units')

for region = 1:3
    cell_type = brain_region_cells(region,:);
    figure
    subplot(2,2,1)
    bar(1,sum(cell_type==1),'k')
    hold on
    bar(2,sum(cell_type==2),'b')
    bar(3,sum(cell_type==3),'r')
    xticks = [1 2 3];
    xticklabels({'Narrow Interneuron','Wide Interneuron','Pyramidal Neuron'})
    ylabel('Number of neurons detected')
    

    subplot(2,2,2)

        scatter3(cell_metrics.troughToPeak(find(cell_type==1)),...
        cell_metrics.firingRate(find(cell_type==1)),cell_metrics.acg_tau_rise(find(cell_type==1)),'k')
    hold on
    scatter3(cell_metrics.troughToPeak(find(cell_type==2)),...
        cell_metrics.firingRate(find(cell_type==2)),cell_metrics.acg_tau_rise(find(cell_type==2)),'b')
    hold on
    scatter3(cell_metrics.troughToPeak(find(cell_type==3)),...
        cell_metrics.firingRate(find(cell_type==3)),cell_metrics.acg_tau_rise(find(cell_type==3)),'r')
    hold on
    xlabel('Trough to Peak (ms)')
    ylabel('firingRate')
    zlabel('ACG Tau rise')

    % Visualise FR difference (which is not used as the criteria for cell classification)
    subplot(2,2,3)
    histogram(cell_metrics.firingRate(find(cell_type==1)),...
        'EdgeColor','none','FaceColor','k','Normalization','cdf','BinEdges',0:0.1:30)
    hold on
    histogram(cell_metrics.firingRate(find(cell_type==2)),...
        'EdgeColor','none','FaceColor','b','Normalization','cdf','BinEdges',0:0.1:30)
    histogram(cell_metrics.firingRate(find(cell_type==3)),...
        'EdgeColor','none','FaceColor','r','Normalization','cdf','BinEdges',0:0.1:30)
    legend('Narrow Interneruon','Wide Interneuron','Pyramidal Neuron')
    xlabel('Firing Rate')

    % Visualise FR difference (which is not used as the criteria for cell classification)
    subplot(2,2,4)
    histogram(cell_metrics.burstIndex_Mizuseki2012(find(cell_type==1)),...
        'EdgeColor','none','FaceColor','k','Normalization','cdf','BinEdges',0:0.01:1)
    hold on
    histogram(cell_metrics.burstIndex_Mizuseki2012(find(cell_type==2)),...
        'EdgeColor','none','FaceColor','b','Normalization','cdf','BinEdges',0:0.01:1)
    histogram(cell_metrics.burstIndex_Mizuseki2012(find(cell_type==3)),...
        'EdgeColor','none','FaceColor','r','Normalization','cdf','BinEdges',0:0.01:1)
    xlabel('Fraction of burstiness (less than 6 ms)')

    
    sgtitle(title_texts{region})
end


session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
session_info(1).probe(1).importMode = 'KS';
column = 1

shanks_available = unique(chan_config.Shank);
for n=1:size(shanks_available)
    cols_available = unique(chan_config.Ks_xcoord(chan_config.Shank==shanks_available(n)));
end

 [file_to_use,imecMeta, chan_config, sorted_config] = extract_NPX_channel_config(session_info(1).probe(1),column);
for column = 1:4
    col_ID = cols_available(column); % (1) is 11
    sorted_config_spikes = sortrows(chan_config,'Ks_ycoord','descend');
    col_idx = sorted_config_spikes.Ks_xcoord == col_ID;
    sorted_config_spikes = sorted_config_spikes(col_idx,:);

    for nchannel = 1:size(sorted_config_spikes,1)
        narrow_int(nchannel,column) = sum(SUA(sorted_config_spikes.Channel(nchannel)).cell_type == 1);
        wide_int(nchannel,column) = sum(SUA(sorted_config_spikes.Channel(nchannel)).cell_type == 2);
        pyramidal(nchannel,column) = sum(SUA(sorted_config_spikes.Channel(nchannel)).cell_type == 3);
    end
end

figure
hold on
plot(sum(narrow_int,2),sorted_config_spikes.Ks_ycoord-10,'k')
plot(sum(wide_int,2),sorted_config_spikes.Ks_ycoord-10,'b')
plot(sum(pyramidal,2),sorted_config_spikes.Ks_ycoord-10,'r')
% barh(narrow_int,'FaceColor','k','FaceAlpha',0.3)
% barh(wide_int,'FaceColor','b','FaceAlpha',0.3)
% barh(pyramidal,'FaceColor','r','FaceAlpha',0.3)
xlim([0 12])
legend('Narrow interneuon','Wide interneuron','Pyramidal neuron')
% yticklabels(sorted_config_spikes.Ks_ycoord(1:2:end))



%% 3.3 Work with several sessions (batch-mode)

basepaths = {'/your/data/path/basename1/','/your/data/path/basename2/'};
basenames = {'basename1','basename2'};

cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

cell_metrics = CellExplorer('metrics',cell_metrics);


%% 4.1 NeuroScope2: Two 6-shank silicon probes implanted bilaterally in CA1 (128 channels; 150 cells)

cd('/Volumes/Peter_SSD_4/CellExplorerTutorial/MS22/Peter_MS22_180629_110319_concat');
NeuroScope2


%% 4.2 NeuroScope2: Inspect Neuropixels data

cd('/Volumes/Peter_SSD_4/NeuropixelsData/PP01/PP01_2020-06-29_13-15-57');
NeuroScope2