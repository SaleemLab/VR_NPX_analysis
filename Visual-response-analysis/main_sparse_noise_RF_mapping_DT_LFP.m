% Programme to estiamte receptive fields from sparse noise in NP1 / NP2
% data
% Dependencies: NPXAnalysis2022;
%               Visual-response-analysis/Stimuli/sparsenoise
% cd('/research/USERS/Masa/code')



%% Visual tuning based on SparseNoise
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

clear all
% SUBJECTS = {'M23017','M23028','M23029'};
% SUBJECTS = {'M23087'};
SUBJECTS = {'M24073'};

% SUBJECTS = {'M23153'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS,'V1-MEC');
% Stimulus_type = 'Masa2tracks';
ROOTPATH = 'Z:\ibn-vision';
% Stimulus_type = 'SparseNoise_fullscreen';
Stimulus_type = 'SparseNoise';
initMap = [];
probe_hemisphere_text = {'MEC','V1'};

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    %     cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info.probe(1).SUBJECT,'ephys',session_info.probe(1).SESSION,'analysis'))
    %     load(fullfile(options.ANALYSIS_DATAPATH,'..',"best_channels.mat"))

    for n = 1:length(session_info)
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            options = session_info(n).probe(nprobe);

            load(fullfile(options.ANALYSIS_DATAPATH,"extracted_behaviour.mat"))
            load(fullfile(options.ANALYSIS_DATAPATH,"extracted_task_info.mat"))


            bin_size = 1/1250;
            AnalysisTimeWindow = [-0.02 0.1];
            options.importMode = 'KS'; % LF or MUA or KS
            % options.importMode = 'LF'; % LF or MUA or KS
            options.BinWidth = bin_size; % resolution (in s) of output resps (e.g. 1/60)
            options.stim_dur = 0.1;
            options.AnalysisTimeWindow = AnalysisTimeWindow;% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
            options.ks_unitType = 'good'; % 'mua', 'good' or ''
            options.PD_FLAG = 1;
            options.paradigm = 'masa';
            %         options.gFileNum = gFileNum;
            options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

            kClusters = kmeans(chan_config.Ks_ycoord,2);
            if mean(chan_config.Ks_ycoord(kClusters==1)) > mean(chan_config.Ks_ycoord(kClusters==2))
                V1_channel_ids = kClusters==1;

            elseif mean(chan_config.Ks_ycoord(kClusters==2)) > mean(chan_config.Ks_ycoord(kClusters==1))
                V1_channel_ids = kClusters==2;

            end

            stimTimes = Task_info.stim_onset;
            [resps] = perievent_LFP_response(stimTimes,AnalysisTimeWindow,options);
            tvec = linspace(AnalysisTimeWindow(1),AnalysisTimeWindow(2),size(resps,2));

            % Calculate sparsenoise RFs
            stim_matrix = cat(3,Task_info.stim_matrix{:}); % N rows x M cols x nFrames
            RF_map = zeros(size(stim_matrix,1),size(stim_matrix,2),size(resps,1),size(resps,2));
            for x = 1:size(stim_matrix,1)
                for y = 1:size(stim_matrix,2)
                    this_location = squeeze(stim_matrix(x,y,:));
                    On_trials = find(this_location==1); % find white quad
                    for nchannel = 1:size(resps,1)
                        RF_map(x,y,nchannel,:) = squeeze(mean(resps(nchannel,:,On_trials),3));
                        peak_map(x,y,nchannel) = min(squeeze(mean(resps(nchannel,:,On_trials),3)));
                    end
                end
            end

            % RF(nprobe).resps = resps;
            RF(nprobe).RF_map = RF_map;
            RF(nprobe).peak_map = peak_map;
            RF(nprobe).V1_channel_ids = V1_channel_ids;
            RF(nprobe).shanks = chan_config.Shank;
            RF(nprobe).xcoord = chan_config.Ks_xcoord;
            RF(nprobe).ycoord = chan_config.Ks_ycoord;

            %         save(fullfile(options.ANALYSIS_DATAPATH,'receptiveFields_ks3.mat'),'RF')
            %         save(fullfile(options.ANALYSIS_DATAPATH,'receptiveFields_without_pd.mat'),'RF')
        end
        save(fullfile(options.ANALYSIS_DATAPATH,'receptiveFields_LFP.mat'),'RF')
    end
end
% load(fullfile(EPHYS_DATAPATH,'analysis','receptiveFields.mat'))

imagesc(mean(peak_map(:,:,V1_channel_ids),3))

%% plotting SparseNoise
clear all
SUBJECTS = {'M24073'};
Dates = {'20250124'};
nsession = 1;
% experiment_info = subject_session_stimuli_mapping(SUBJECTS,'bilateral');
% session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,'SparseNoise'));
% options = session_info(nsession).probe(1);
load(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',SUBJECTS{1},'analysis',Dates{nsession},'SparseNoise_4','receptiveFields_LFP.mat'),'RF')
load(fullfile('Z:\ibn-vision\DATA\SUBJECTS\',SUBJECTS{1},'analysis',Dates{nsession},'SparseNoise_4','session_info.mat'))
% hemisphere_texts = {'Left','Right'}

for nprobe = 1:2
    options = session_info.probe(nprobe);
    fig = figure
    fig.Position = [34 60 1850 920]
    fig.Name = sprintf('%s %s Averaged V1-MEC-paraSub Receptive field',options.SUBJECT,options.SESSION);
    sgtitle(sprintf('%s %s Averaged V1 Receptive field',options.SUBJECT,options.SESSION))
    RF_map = RF(nprobe).peak_map;
    % RF(nprobe)
    options.importMode = 'KS';
     [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);

    % V1_channel_ids=RF(nprobe).V1_channel_ids;
    % imagesc(mean(RF(nprobe).peak_map(:,:,V1_channel_ids),3))
    temp = round(linspace(max(chan_config.Ks_ycoord),min(chan_config.Ks_ycoord),37));
    channel_depths_range=[];
    channel_depths_range(:,1) = temp(1:end-1);
    channel_depths_range(:,2) = temp(2:end);

    % channel_depths_range = [3200 2700;2700 2200;2200 1700;1700 1300];

    for nregion = 1:36
        subplot(6,6,nregion)
        % scal_f = 10; % scale image by this before...
        % sigma = 3; % ...filtering by this
        scal_f = 2; % scale image by this before...
        sigma = 1/3; % ...filtering by this
        % V1_channel_ids = find(RF(nprobe).shanks == nshank&RF(nprobe).V1_channel_ids==1);
        V1_channel_ids=find(chan_config.Ks_ycoord<channel_depths_range(nregion,1) & chan_config.Ks_ycoord>channel_depths_range(nregion,2));

        thisMap_s = imresize(squeeze(nanmean(RF_map(:,:,V1_channel_ids),3)),scal_f);
        thisMap_s = imgaussfilt(thisMap_s,sigma);
        %     thisMap_s = zscore(thisMap_s,0,'all');
        imagesc(flip(thisMap_s))


        xlabel('Azimuth')
        ylabel('Elevation')
        xticks(linspace(1,size(thisMap_s,2),13))
        xticklabels(-120:20:120)
        yticks(linspace(1,size(thisMap_s,1),7))
        yticklabels(flip(-30:20:90))
        title(sprintf('Averaged RF for depths %i micron to %i micron',channel_depths_range(nregion,1),channel_depths_range(nregion,2)))
        colorbar
        colormap((gray))
        set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
    end
end



figure
n = 1;
RF_map= RF(nprobe).RF_map;
peak_map = RF(nprobe).peak_map;

for x = 1:8
    for y = 1:16
        subplot(8,16,n)
        if peak_map(x,y,233) == min(min(squeeze(peak_map(:,:,233))))
            plot(squeeze(RF_map(x,y,233,:)),'r')
        else
            plot(squeeze(RF_map(x,y,233,:)),'k')
        end
        ylim([-0.0001 0.0001])
        n = n+1;
    end
end

%% 
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