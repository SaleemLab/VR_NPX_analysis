%% Kilosort + cell explorer basic visualisation of cluster spiking
% Set the data folders and processing parameters
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))

if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
%     ROOTPATH = 'X:\ibn-vision';
    ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive
%     ROOTPATH = '/research';
end

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'RUN';

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info(n).probe(1).SUBJECT,'ephys',session_info(n).probe(1).SESSION,'analysis'))
        load best_channels
        load extracted_PSD

        raw_LFP = [];
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            column = 1;
            session_info(n).probe(nprobe).task_type = stimulus_name{n};
            options = session_info(n).probe(nprobe);
            %             options.importMode = 'LF';
            options.probe_no = nprobe; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            options.importMode = 'KS';
            options.gFileNum = gFileNum(n);
            % Load all spike data sorted according to the channel position

            [SUA.probe{nprobe} chan_config sorted_config] = load_KS_NPX1(options,column,'LFP_tvec',[],'group','by channel','cell_exporer','on');

            load([options.KS_DATAPATH,'\',options.SUBJECT,'_',options.SESSION,'.cell_metrics.cellinfo.mat'])
            load([options.KS_DATAPATH,'\','cluster_table.mat'])


            cell_metrics.KS_cluster_id = double(cluster_table.cluster_id)';
            cell_metrics.good_unit = double(cluster_table.cluster_id)';

            cols_available = unique(chan_config.Ks_xcoord);

            for column = 1:length(cols_available)
                col_ID = cols_available(column); % (1) is 11
                sorted_config_spikes = sortrows(chan_config,'Ks_ycoord','descend');
                col_idx = sorted_config_spikes.Ks_xcoord == col_ID;
                sorted_config_spikes = sorted_config_spikes(col_idx,:);

                for nchannel = 1:size(sorted_config_spikes,1)
                    narrow_int(nchannel,column) = sum(SUA.probe{nprobe}(sorted_config_spikes.Channel(nchannel)).cell_type == 1);
                    wide_int(nchannel,column) = sum(SUA.probe{nprobe}(sorted_config_spikes.Channel(nchannel)).cell_type == 2);
                    pyramidal(nchannel,column) = sum(SUA.probe{nprobe}(sorted_config_spikes.Channel(nchannel)).cell_type == 3);
                end
            end

            fig = figure
            fig.Position = [334.7143 102.7143 1100 830]

            subplot(1,3,1);
            colour_line= {[177,0,38]/256,[213,62,79]/256,[252,141,89]/256,...
                [153,213,148]/256,[26,152,80]/256,[66,146,198]/256,[8,69,148]/256};
            freq_legends = {'0.5 -3 Hz','4-12 Hz','9 - 17 Hz','30-60 Hz','60-100 Hz','125-300 Hz','300-600 Hz'};
            selected_frequency = [1 2 6 7];
            %     subplot(1,4,2)
            hold on
            plot([0 1],[chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel) chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel)],'--k','LineWidth',2)

            if ~isempty(best_channels{nprobe}.L4_channel)
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-5)],'--b','LineWidth',0.5)
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)+5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)+5)],'--b','LineWidth',0.5)
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel))],'--b','LineWidth',2)
            end

            plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel)-10)],'--c','LineWidth',0.5)
            if isempty(best_channels{nprobe}.L4_channel) | find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10 < find(chan_config.Channel ==best_channels{nprobe}.L4_channel)-6
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10)],'--c','LineWidth',0.5)
            else
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L4_channel)-6) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-6)],'--c','LineWidth',0.5)
            end
            plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel))],'--c','LineWidth',2)


            if ~isempty(best_channels{nprobe}.CA1_channel)
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)-10)],'--r','LineWidth',0.5)
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)+10)],'--r','LineWidth',0.5)
                plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel))],'--r','LineWidth',2)
            end

            for i = selected_frequency
                pl(i) = plot(power{nprobe}(:,i)./max(power{nprobe}(:,i)),sorted_config.Ks_ycoord','Color',colour_line{i})

                %     p(n) = plot(power_differnece,1:96,colour_line{n})
                hold on
            end
            % legend('1-3','4-12','30-60','60-100','125-300')
            ylabel('probe depth (um)')
            legend([pl(selected_frequency)],{freq_legends{selected_frequency}},'Location','southeast','Color','none');
            ylim([0 4000])
            set(gca,"TickDir","out",'box', 'off','Color','none')

            subplot(1,3,2)
            hold on
%             plot(sum(narrow_int,2),sorted_config_spikes.Ks_ycoord,'k')
%              plot(sum(wide_int,2),sorted_config_spikes.Ks_ycoord,'b')
%               plot(sum(pyramidal,2),sorted_config_spikes.Ks_ycoord,'r')
%             barh(sorted_config_spikes.Ks_ycoord,sum(narrow_int,2),'k','FaceAlpha',0.5,'EdgeColor','none')
%             barh(sorted_config_spikes.Ks_ycoord,sum(wide_int,2),'b','FaceAlpha',0.5,'EdgeColor','none')
%             barh(sorted_config_spikes.Ks_ycoord,sum(pyramidal,2),'r','FaceAlpha',0.5,'EdgeColor','none')
            b = barh(sorted_config_spikes.Ks_ycoord,[sum(narrow_int,2) sum(wide_int,2) sum(pyramidal,2)],'stacked','FaceAlpha',0.5,'EdgeColor','none')
            xmax = max(sum([sum(narrow_int,2) sum(wide_int,2) sum(pyramidal,2)]'));
            b(1).FaceColor = 'k';
            b(2).FaceColor = 'b';
            b(3).FaceColor = 'r';

            % barh(narrow_int,'FaceColor','k','FaceAlpha',0.3)
            % barh(wide_int,'FaceColor','b','FaceAlpha',0.3)
            % barh(pyramidal,'FaceColor','r','FaceAlpha',0.3)

            hold on
            plot([0 12],[chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel) chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel)],'--k','LineWidth',2)

            if ~isempty(best_channels{nprobe}.L4_channel)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-5)],'--b','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)+5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)+5)],'--b','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel))],'--b','LineWidth',2)
            end

            plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel)-10)],'--c','LineWidth',0.5)
            if isempty(best_channels{nprobe}.L4_channel) | find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10 < find(chan_config.Channel ==best_channels{nprobe}.L4_channel)-6
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10)],'--c','LineWidth',0.5)
            else
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L4_channel)-6) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-6)],'--c','LineWidth',0.5)
            end
            plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel))],'--c','LineWidth',2)


            if ~isempty(best_channels{nprobe}.CA1_channel)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)-10)],'--r','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)+10)],'--r','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel))],'--r','LineWidth',2)
            end

            xlim([0 xmax])
            xticks(0:2:xmax)
            legend('Narrow interneuon','Wide interneuron','Pyramidal neuron','Color','none')
            % yticklabels(sorted_config_spikes.Ks_ycoord(1:2:end))
            xlabel('No of units')
            ylabel('probe depth')
            yyaxis right
            ylim([-chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel) 4000-chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel)])
            set(gca,"TickDir","out",'box', 'off','Color','none')


            subplot(1,3,3)
            hold on
            plot(sum(narrow_int,2),sorted_config_spikes.Ks_ycoord,'k')
             plot(sum(wide_int,2),sorted_config_spikes.Ks_ycoord,'b')
              plot(sum(pyramidal,2),sorted_config_spikes.Ks_ycoord,'r')
%             barh(sorted_config_spikes.Ks_ycoord,sum(narrow_int,2),'k','FaceAlpha',0.5,'EdgeColor','none')
%             barh(sorted_config_spikes.Ks_ycoord,sum(wide_int,2),'b','FaceAlpha',0.5,'EdgeColor','none')
%             barh(sorted_config_spikes.Ks_ycoord,sum(pyramidal,2),'r','FaceAlpha',0.5,'EdgeColor','none')
%             b = barh(sorted_config_spikes.Ks_ycoord,[sum(narrow_int,2) sum(wide_int,2) sum(pyramidal,2)],'stacked','FaceAlpha',0.5,'EdgeColor','none')
            xmax = max(sum([sum(narrow_int,2) sum(wide_int,2) sum(pyramidal,2)]'));
            b(1).FaceColor = 'k';
            b(2).FaceColor = 'b';
            b(3).FaceColor = 'r';

            % barh(narrow_int,'FaceColor','k','FaceAlpha',0.3)
            % barh(wide_int,'FaceColor','b','FaceAlpha',0.3)
            % barh(pyramidal,'FaceColor','r','FaceAlpha',0.3)

            hold on
            plot([0 12],[chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel) chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel)],'--k','LineWidth',2)

            if ~isempty(best_channels{nprobe}.L4_channel)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-5)],'--b','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)+5) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)+5)],'--b','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel))],'--b','LineWidth',2)
            end

            plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel)-10)],'--c','LineWidth',0.5)
            if isempty(best_channels{nprobe}.L4_channel) | find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10 < find(chan_config.Channel ==best_channels{nprobe}.L4_channel)-6
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)+10)],'--c','LineWidth',0.5)
            else
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L4_channel)-6) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)-6)],'--c','LineWidth',0.5)
            end
            plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel))],'--c','LineWidth',2)


            if ~isempty(best_channels{nprobe}.CA1_channel)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)-10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)-10)],'--r','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)+10) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)+10)],'--r','LineWidth',0.5)
                plot([0 xmax],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel))],'--r','LineWidth',2)
            end

            xlim([0 xmax])
            legend('Narrow interneuon','Wide interneuron','Pyramidal neuron','Color','none')
            % yticklabels(sorted_config_spikes.Ks_ycoord(1:2:end))
            xticks(0:2:xmax)
            xlabel('No of units')
            ylabel('probe depth')
            yyaxis right
            ylim([-chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel) 4000-chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel)])
            set(gca,"TickDir","out",'box', 'off','Color','none')



            sgtitle(sprintf('%s %s cell explorer unit classification probe %i',options.SUBJECT,options.SESSION,nprobe))
            filename = sprintf('%s %s cell explorer unit classification probe %i.pdf',options.SUBJECT,options.SESSION,nprobe)
            saveas(gcf,filename)
            filename = sprintf('%s %s cell explorer unit classification probe %i.fig',options.SUBJECT,options.SESSION,nprobe)
            saveas(gcf,filename)
            close all
        end

    end
end

%%
cd(options.KS_DATAPATH)
load M22069_20221201.cell_metrics.cellinfo.mat
load cluster_table

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

        if cell_metrics.maxWaveformCh(cell) < best_channels{nprobe}.CA1_channel+10 & cell_metrics.maxWaveformCh(cell) > best_channels{nprobe}.CA1_channel-10
            CA1_cell_type(cell) = cell_type_all(cell);
        else
            CA1_cell_type(cell) = 0;
        end

        if cell_metrics.maxWaveformCh(cell) < best_channels{nprobe}.first_in_brain_channel & cell_metrics.maxWaveformCh(cell) > best_channels{nprobe}.first_in_brain_channel-100
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

save cell_metrics cell_metrics
%%
[AP_FILE,LF_FILE] = findImecBinFile(options.EPHYS_DATAPATH);

imecMeta = NPadmin.ReadMeta(AP_FILE,options.EPHYS_DATAPATH);
chan_config = NPadmin.getNPChannelConfig(imecMeta);

switch str2double(imecMeta.imDatPrb_type)
    case 0 
    electrode_spacing_um = 20;
    otherwise
    electrode_spacing_um = 15;
end

% Set probe columns (4 columns for NPX1)
% (x coordinate = 11, 27, 43 or 59 micron)
shanks_available = unique(chan_config.Shank);
for n=1:size(shanks_available)
    cols_available = unique(chan_config.Ks_xcoord(chan_config.Shank==shanks_available(n)));
end

% Get channels from one column (x coordinate = 11, 27, 43 or 59 micron)
col_ID = cols_available(1); % (1) is 11
KS_metrics = readtable(fullfile(options.KS_DATAPATH,'metrics.csv'));
[these_spike_times,nominal_KSLabel,cluster_id,peakChannel,maxSpkTime] = import_ks_spiketimes(options.KS_DATAPATH,options.gFileNum,options.KS_CATGT_FNAME,imecMeta.imSampRate)
mean_waveforms = readNPY([options.KS_DATAPATH,'/mean_waveforms.npy']); % Information about mean waveform


MUA = [];
SUA = [];
tic
for nchannel = 1:size(chan_config,1)
    clusters_this_channel = find(peakChannel == nchannel)-1;
    [~,index,]= intersect(cluster_id,clusters_this_channel);
    good_units_index = index(find(nominal_KSLabel(index)=='good'));
    %     good_units_this_channel = clusters_this_channel(find(nominal_KSLabel(index)=='good'));
    SUA(nchannel).spike_times= [];
    SUA(nchannel).spike_id = [];
    SUA(nchannel).cell_type  = [];

    if ~isempty(index) % If any clusters
        if ~isempty(good_units_index) % If any SUA
            SUA(nchannel).spike_times= [];
            SUA(nchannel).spike_id = [];
            SUA(nchannel).cell_type  = [];
            for unit = 1:length(good_units_index)

                SUA(nchannel).spike_times = [SUA(nchannel).spike_times; these_spike_times{good_units_index(unit)}];
                SUA(nchannel).spike_id = [SUA(nchannel).spike_id; good_units_index(unit)*ones(length(these_spike_times{good_units_index(unit)}),1)];
                SUA(nchannel).cell_type = [SUA(nchannel).cell_type  cell_metrics.cell_type(good_units_index(unit))]; % Putative cell types
            end
        end
    end
end
toc


%% Distribution of cell types
% Sort channel for spike time data

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

%% Latency


events = slow_waves.ints.UP(:,1);
events = ripples.onset;
events = V1_ripples.onset;
% events = MousePos.stimuli_onset(MousePos.stimuli_track == 1)
% events = MousePos.stimuli_onset(MousePos.stimuli_track == 2)

cell_id = unique(V1_spike_times(:,1));
% cell_id = unique(V1_SWR_spike_times(:,1));

% events = ripples.onset;
windows = [-0.5 1];
SUA_channel_psth= [];
SUA_channel_latency = [];
bin_size = 0.001;

MUA_filter_length = 50;
SD_alpha = 5; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w=w./sum(w);  %gaussian kernel 10 ms STD, 2 std width

MUA_filter_length = 10;
SD_alpha = 2; %2 std width
MUA_filter_alpha = (MUA_filter_length-1)/SD_alpha;
w1=gausswin(MUA_filter_length,MUA_filter_alpha); %41,4
w1=w1./sum(w1);  %gaussian kernel 10 ms STD, 2 std width

events = slow_waves.ints.UP(:,1);
events = slow_waves.ints.UP(:,1);
events = ripples.SWS_onset;
events = V1_ripples.onset;
events = spindles.onset;
figure
SUA_channel_latency = [];
SUA_channel_psth= [];
for column = 1:4

%     SUA_channel_psth = [];
    col_ID = cols_available(column); % (1) is 11
    sorted_config_spikes = sortrows(chan_config,'Ks_ycoord','descend');
    col_idx = sorted_config_spikes.Ks_xcoord == col_ID;
    sorted_config_spikes = sorted_config_spikes(col_idx,:);

    for nchannel = 1:size(sorted_config_spikes,1)
       spikes_this_channel = SUA( sorted_config_spikes.Channel(nchannel)).spike_times;
%         spikes_this_channel = SUA(nchannel).spike_times;
        %     mean_FR(cell) = length(spikes_this_cell)/(spikes_this_cell(end)-spikes_this_cell(1));

        windows = [0 0.04];
        [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_channel, events, windows, bin_size);
        for event = 1:length(events)
            if ~isempty(find(binnedArray(event,:) == 1))
                SUA_channel_latency(nchannel,event) = median(find(binnedArray(event,:) == 1));
            else
                SUA_channel_latency(nchannel,event) = nan;
            end
        end
        SUA_channel_psth(nchannel,:) = sum(binnedArray)./max(sum(binnedArray));
    end

%     subplot(2,2,1)
%     imagesc(bins,1:1:size(SUA_channel_psth,1),SUA_channel_psth)
%     % imagesc(bins,1:1:length(index),zscore(V1_UP_psth(index,:),0,2))
%     ylabel('Depth')
%     xlabel('Time relative to event onset (s)')
%     cbar = colorbar
%     cbar.Label.String = 'Firing Rate (z score)';
%     hold on
%     plot([0 0],get(gca,'ylim'),'r','LineWidth',5)
% %     clim([0 4])
%     yticks(1:4:size(SUA_channel_psth,1))
%     yticklabels(sorted_config_spikes.Ks_ycoord(1:4:end))
%         ylim([1 40])
% 
    subplot(2,2,column)
    hold on
    scatter(nanmean(SUA_channel_latency,2),sorted_config_spikes.Ks_ycoord,'k','filled','o','MarkerFaceAlpha',0.3)
    ylim([2300 3500])
    xlim([0 40])
    title(sprintf('Column %i/4',column))
end
sgtitle('Down-UP transition median MUA latency')
sgtitle('Spindle median MUA latency')
