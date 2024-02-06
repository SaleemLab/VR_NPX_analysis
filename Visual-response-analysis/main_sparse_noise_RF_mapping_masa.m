% Programme to estiamte receptive fields from sparse noise in NP1 / NP2
% data
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise
% cd('/research/USERS/Masa/code')



%% Visual tuning based on SparseNoise
clear all
SUBJECTS = {'M23017','M23028','M23029'};
% SUBJECTS = {'M23087'};

SUBJECTS = {'M23153'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_type = 'Masa2tracks';
ROOTPATH = 'Z:\ibn-vision';
Stimulus_type = 'SparseNoise_fullscreen';
initMap = [];
probe_hemisphere_text = {'Left','Right'};

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    %     cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info.probe(1).SUBJECT,'ephys',session_info.probe(1).SESSION,'analysis'))
    load('best_channels.mat')
    for n = 1:length(session_info)
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            options = session_info(n).probe(nprobe);
            bin_size = 1/60;
            AnalysisTimeWindow = [0 1/60*7];


            options.importMode = 'KS'; % LF or MUA or KS
            % options.importMode = 'LF'; % LF or MUA or KS
            options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
            options.stim_dur = 0.1;
            options.AnalysisTimeWindow = [0 1/60*7];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
            options.ks_unitType = 'good'; % 'mua', 'good' or ''
            options.PD_FLAG = 1;
            options.paradigm = 'masa';
            %         options.gFileNum = gFileNum;
            options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

            load(fullfile(options.ANALYSIS_DATAPATH,"extracted_behaviour.mat"))
            load(fullfile(options.ANALYSIS_DATAPATH,"extracted_task_info.mat"))
            load(fullfile(options.ANALYSIS_DATAPATH,"extracted_clusters.mat"))

            % currently hard-coded
            good_unit_index = (clusters.amplitude_cutoff <= 0.1...
                &clusters.isi_viol <= 0.1...
                &clusters.amplitude >=50);
            good_unit = clusters.cluster_id(good_unit_index);

            all_spike_data = [clusters.spike_id(find(ismember(clusters.spike_id, good_unit))) clusters.spike_times(find(ismember(clusters.spike_id, good_unit)))];

            %         num_cell = length(unique(spike_data(:,1)));
            for unit_id = 1:length(good_unit)
                spikes_this_cell = all_spike_data(find(all_spike_data(:,1) == good_unit(unit_id)),2);

                [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spikes_this_cell, Task_info.stim_onset, AnalysisTimeWindow, bin_size);
                resps(unit_id,:,:) = binnedArray';
            end

            % Calculate sparsenoise RFs
            stim_matrix = cat(3,Task_info.stim_matrix{:}); % N rows x M cols x nFrames

            sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
            sn_options.mapSampleRate = 60; % Hz
            sn_options.mapsToShow = {'linear','black','white','contrast'};
            sn_options.mapMethod = 'fitlm'; % fitlm mean
            sn_options.framesToShow = 1:7;  % at 60 Hz would be 8 ms, 40, 72 etc
            sn_options.plotflag = 0; % 1 if you want to see details of the receptive fields

            horizontalAngles = -30:(30+90)/8:90;
            verticalAngles = -120:240/16:120;

            SUA_data = [];
            On_response = [];
            Off_response = [];
            average_On = [];
            average_Off = [];
            average_contrast = [];
            contrast_response = [];
            %         On_grid_distribution = zeros(size(stim_matrix,[1 2]));
            %         Off_grid_distribution = zeros(size(stim_matrix,[1 2]));
            %
            %         % x 240 degree [-120 120]
            %         % y 120 degree [-30 90]
            %         for x = 1:size(stim_matrix,1)
            %                             for y = 1:size(stim_matrix,2)
            %                                 this_location = squeeze(stim_matrix(x,y,:));
            %                                 On_trials = find(this_location==1); % find white quad
            %                                 On_grid_distribution(x,y) = sum(this_location == 1);
            %                                 Off_grid_distribution(x,y) = sum(this_location == -1);
            %
            %                             end
            %         end
            %         Overall_grid_distribution = On_grid_distribution + Off_grid_distribution;
            %         % check on distribution of stimuli on the grid
            %         figure;
            %         subplot(3,1,1);imagesc(On_grid_distribution');colorbar;caxis([0 max(max(On_grid_distribution))])
            %         subplot(3,1,2);imagesc(Off_grid_distribution');colorbar;caxis([0 max(max(Off_grid_distribution))])
            %         subplot(3,1,3);imagesc(Overall_grid_distribution');colorbar;caxis([0 max(max(Overall_grid_distribution))])

            unit_no = 0;
            for unit_id = 1:size(resps,1)

                %                 if unit_id> 300 || unit_id <600
                %                     sn_options.plotflag = 1;
                %                 end
                SUA_data = squeeze(resps(unit_id,:,:))'; % returns 1 x N time bins x nFrames;
                initMap_temp = sparseNoiseAnalysis(stim_matrix,SUA_data,[],[],sn_options);
                initMap{unit_id,1} = initMap_temp;

            end

            RF.probe(options.probe_no).cluster_id = good_unit;
            RF.probe(options.probe_no).peak_channel = clusters(options.probe_no).peak_channel(good_unit_index);
            RF.probe(options.probe_no).peak_location = chan_config.Ks_ycoord(RF.probe(options.probe_no).peak_channel);
            RF.probe(options.probe_no).shank = chan_config.Shank(RF.probe(options.probe_no).peak_channel);

            RF.probe(options.probe_no).RF_map = initMap;
            initMap = [];
        end

        save(fullfile(options.ANALYSIS_DATAPATH,'receptiveFields.mat'),'RF')
%         save(fullfile(options.ANALYSIS_DATAPATH,'receptiveFields_without_pd.mat'),'RF')


        fig = figure
        fig.Position = [300 150 1000 620]
        fig.Name = sprintf('%s %s Averaged V1 Receptive field',options.SUBJECT,options.SESSION);


        for nprobe = 1:length(session_info.probe) % For each session, how many probes
            options = session_info.probe(nprobe);
            options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            options.importMode = 'KS';

            [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

            for unit_id = 1:length(RF.probe(options.probe_no).cluster_id)

                initMap_temp = RF.probe(options.probe_no).RF_map{unit_id}(:,:,:,4);

                for x = 1:size(initMap_temp,1)
                    for y = 1:size(initMap_temp,2)
                        RF_map(unit_id,x,y) = mean(initMap_temp(x,y,4:7));
                    end
                end
                temp = RF_map(unit_id,:,:)/max(max(RF_map(unit_id,:,:)));
                temp(isinf(temp)) = nan;
                RF_map(unit_id,:,:) = temp;
                %         imagesc(initMap{options.probe_no}{100}(:,:,4,4))
            end



            %         c = 1;
            %         for unit = 200:239
            %
            %             subplot(5,8,c)
            %             imagesc(flip(squeeze(RF_map(unit,:,:))))
            %             c = c + 1;
            %         end
            subplot(2,2,nprobe)
            scal_f = 10; % scale image by this before...
            sigma = 3; % ...filtering by this
            V1_cell_id = find(RF.probe(options.probe_no).peak_location <= chan_config.Ks_ycoord(best_channels{options.probe_no}.first_in_brain_channel) & RF.probe(options.probe_no).peak_location >= chan_config.Ks_ycoord(best_channels{options.probe_no}.first_in_brain_channel)-800);
            
            
            thisMap_s = imresize(squeeze(nanmean(RF_map(V1_cell_id,:,:))),scal_f);
            thisMap_s = imgaussfilt(thisMap_s,sigma);
            imagesc(flip(thisMap_s/max(max(thisMap_s))))
            xlabel('Azimuth')
            ylabel('Elevation')
            xticks(linspace(1,size(thisMap_s,2),13))
            xticklabels(-120:20:120)
            yticks(linspace(1,size(thisMap_s,1),7))
            yticklabels(flip(-30:20:90))
            title(sprintf('Averaged %s V1 Receptive field',probe_hemisphere_text{options.probe_hemisphere}))
            colorbar
            set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
        end
        mkdir('visual')
        save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','visual'),[])
    end
end

% load(fullfile(EPHYS_DATAPATH,'analysis','receptiveFields.mat'))

%% plotting SparseNoise
fig = figure
fig.Position = [300 150 1000 620]
fig.Name = sprintf('%s %s Averaged V1 Receptive field',options.SUBJECT,options.SESSION);
RF_map = [];
options.probe_no = 1
for unit_id = 1:length(RF.probe(options.probe_no).cluster_id)

    initMap_temp = RF.probe(options.probe_no).RF_map{unit_id}(:,:,:,4);

    for x = 1:size(initMap_temp,1)
        for y = 1:size(initMap_temp,2)
            RF_map(unit_id,x,y) = mean(initMap_temp(x,y,4:7));
        end
    end
%     temp = zscore(RF_map(unit_id,:,:),0,"all");
     temp = RF_map(unit_id,:,:);
    temp(isinf(temp)) = nan;
    RF_map(unit_id,:,:) = temp;
    %         imagesc(initMap{options.probe_no}{100}(:,:,4,4))
end


for nshank = 1:4
    subplot(2,2,nshank)
    scal_f = 10; % scale image by this before...
    sigma = 3; % ...filtering by this
    V1_cell_id = find(RF.probe(options.probe_no).shank == nshank);

    thisMap_s = squeeze(nanmean(RF_map(V1_cell_id,:,:)));

    thisMap_s = imresize(squeeze(nanmean(RF_map(V1_cell_id,:,:))),scal_f);
    thisMap_s = imgaussfilt(thisMap_s,sigma);
%     thisMap_s = zscore(thisMap_s,0,'all');
    imagesc(flip(thisMap_s))


    xlabel('Azimuth')
    ylabel('Elevation')
    xticks(linspace(1,size(thisMap_s,2),13))
        xticklabels(-120:20:120)
    yticks(linspace(1,size(thisMap_s,1),7))
        yticklabels(flip(-30:20:90))
    title(sprintf('Averaged Shank %i V1 Receptive field',nshank))
    colorbar
    set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
end



for nshank = 1:4
    %     subplot(2,2,nshank)
    scal_f = 10; % scale image by this before...
    sigma = 3; % ...filtering by this
    V1_cell_id = find(RF.probe(options.probe_no).shank == nshank);
    fig = figure
    fig.Position = [300 150 1000 620]
    fig.Name = sprintf('%s %s Shank %i Averaged V1 Receptive field',options.SUBJECT,options.SESSION,nshank);

    for unit = 1:16

        subplot(4,4,unit)
        thisMap_s = squeeze(RF_map(V1_cell_id(unit),:,:));
        thisMap_s = imresize(squeeze(thisMap_s),scal_f);
        thisMap_s = imgaussfilt(thisMap_s,sigma);
        %     thisMap_s = zscore(thisMap_s,0,'all');
        imagesc(flip(thisMap_s))


        xlabel('Azimuth')
        ylabel('Elevation')
        xticks(linspace(1,size(thisMap_s,2),13))
        xticklabels(-120:20:120)
        yticks(linspace(1,size(thisMap_s,1),7))
        yticklabels(flip(-30:20:90))
        title(sprintf('Unit %i Shank %i V1 Receptive field',V1_cell_id(unit),nshank))
        colorbar
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
end