%% Putative pipeline for processing NPX1 data recorded from hippocampus and V1
% This is a higher-order multi-purpose pipeline that calls dependent functions for
% different analysis piepline

% For mapping of visual receptive field, Please refer to
% Sparse_noise_RF_mapping_masa.mat (subject to change)

% main_NPX_data_preprocessing

% Go to SUA_analysis_masa for kilosort + cell explorer cell classification
% Go to CellExplorerTest_masa_adapted for cell exploerer pipeline


%% Set the data folders and processing parameters
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
% addpath('Z:\ibn-vision\USERS\Masa\code\Masa_utility')
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\NPXAnalysis\NPXAnalysis2022'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\visual_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\LFP_analysis'));
% addpath(genpath('Z:\ibn-vision\USERS\Masa\code\spikes'));

if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
%     ROOTPATH = 'X:\ibn-vision';
    ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive
%     ROOTPATH = '/research';
end

SUBJECTS = {'M23017','M23028','M23029'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'Masa2tracks';

%% import and align and store Bonsai data
Stimulus_type = 'Masa2tracks';
for nsession =4:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
            session_info(n).probe(nprobe).task_type = stimulus_name{n};
        end

        if length(session_info(n).probe) >1
            align_probes_NX1(session_info(n).probe);
        end

        options = session_info(n).probe(1);
        [Behaviour,position] = import_and_align_Masa_VR_Bonsai(stimulus_name{n},session_info(n).probe);
        
        if exist(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis')) == 0
            mkdir(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
        end
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
   
        save(sprintf('bonsai_behaviour%s.mat',erase(stimulus_name{n},Stimulus_type)),'Behaviour')
        save(sprintf('extracted_position%s.mat',erase(stimulus_name{n},Stimulus_type)),'position')
        
        if contains(stimulus_name{n},'RUN')
            lap_times = extract_laps_masa(1,Behaviour,position)
            save extracted_laps lap_times
        end


        % figure
        % hold on
        % plot(MousePos.sglxTime,MousePos.pos)
        % scatter(MousePos.stimuli_onset(MousePos.stimuli_track == 1),1000*MousePos.stimuli_track(MousePos.stimuli_track == 1),'r')
        % scatter(MousePos.stimuli_onset(MousePos.stimuli_track == 2),100*MousePos.stimuli_track(MousePos.stimuli_track == 2),'b')
    end

end

% 
% Stimulus_type = 'Masa2tracks_POST';
% for nsession =5:length(experiment_info)
%     session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
%     for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
%         for nprobe = 1:length(session_info(n).probe) % For each session, how many probes
%             session_info(n).probe(nprobe).task_type = stimulus_name{n};
%         end
% 
%         if length(session_info(n).probe) >1
%             align_probes_NX1(session_info(n).probe);
%         end
%     end
% end

%% PSD analysis and LFP profile
clear all
Stimulus_type = 'RUN'; % extract LFP during RUN
ROOTPATH = 'Z:\ibn-vision'
SUBJECTS = {'M23087'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
   
    PSD = [];
    power = [];
    best_channels = [];
    for nprobe = 1:length(session_info.probe)
        session_info.probe(nprobe).importMode = 'LF';
        options = session_info.probe(nprobe);

        column = 1;
        if nprobe ~= 1
            % Information about AP band probe 1 sample size to align probe
            % 2 LFP traces.
            session_info.probe(1).importMode = 'KS';
            [~, imecMeta, ~, ~] = extract_NPX_channel_config(session_info.probe(1),column);
            [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column,'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp);
        else
            [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column);
        end

        [PSD{nprobe} power{nprobe} best_channels{nprobe}] = calculate_channel_PSD(raw_LFP,SR,sorted_config,options,'plot_option',1)
        % [gamma_coherence gamma_phase_coherence] = gamma_coherence_analysis(raw_LFP,tvec,SR,best_channels,sorted_config) % This function
        if exist(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis')) == 0
            mkdir(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
        end
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
%         title(sprintf('%s %s PSD profile probe %i',options.SUBJECT,options.SESSION,nprobe))
%         filename = sprintf('%s %s PSD profile probe %i.pdf',options.SUBJECT,options.SESSION,nprobe)
%         saveas(gcf,filename)
%         filename = sprintf('%s %s PSD profile probe %i.fig',options.SUBJECT,options.SESSION,nprobe)
%         saveas(gcf,filename)
    end
    save extracted_PSD PSD power
    save best_channels best_channels
end


%% L4 based on checkerboard
addpath(genpath('Z:\ibn-vision\USERS\Masa\code'))
if ismac
    ROOTPATH = '/Users/s.solomon/Filestore/Research2/ibn-vision';
else
%     ROOTPATH = 'X:\ibn-vision';
    ROOTPATH = 'Z:\ibn-vision'; % New server mapped to z drive
%     ROOTPATH = '/research';
end

% SUBJECTS = {'M23017','M23028','M23029'};
SUBJECTS = {'M23087'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
Stimulus_type = 'Checkerboard';


for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        for nprobe = 1:length(session_info.probe) % For each session, how many probes
            options= session_info(n).probe(nprobe);
            options.ROOTPATH = ROOTPATH;
            options.importMode = 'LF';
            options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

            [lfpAvg.nprobe(options.probe_no),csd.nprobe(options.probe_no),power,best_channels] = determine_best_channels(options,options.probe_no);
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'),[])
            close

            column = 1;
            [LF_FILE imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);% Since it is LF
            [best_channels_updated{options.probe_no}] = update_best_channels(options,sorted_config);
            
            % If updating all
            best_channels{options.probe_no}.first_in_brain_channel = best_channels_updated{options.probe_no}(1);
            best_channels{options.probe_no}.L4_channel = best_channels_updated{options.probe_no}(2);
%             best_channels{nprobe}.L4_channel = [];
            best_channels{options.probe_no}.L5_channel = best_channels_updated{options.probe_no}(3);
            best_channels{options.probe_no}.CA1_channel = best_channels_updated{options.probe_no}(4);
            
            close
            
            save(fullfile(options.ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis',"best_channels.mat"),"best_channels")
            % Replot based on updated channels
            plot_perievent_CSD_LFP(lfpAvg.nprobe(options.probe_no),csd.nprobe(options.probe_no),power{options.probe_no},chan_config,sorted_config,best_channels{options.probe_no},options)
            save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'),[])
        end

    end
end


%% ripple power during immobility and during running (not complete as of 22/09/23)
Stimulus_type = 'RUN'; % extract LFP during RUN
for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,'RUN'));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
    load best_channels
    load extracted_PSD
    load(sprintf('extracted_position%s.mat',erase(stimulus_name{n},'Masa2tracks')))


    for nprobe = 1:length(session_info.probe)
        options = session_info.probe(nprobe);
        options.importMode = 'LF';
        column = 1;
        if nprobe ~= 1
            session_info.probe(1).importMode = 'KS';
            [~, imecMeta, ~, ~] = extract_NPX_channel_config(session_info.probe(1),column);
            [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column,'probe_no',nprobe,'probe_1_total_sample',imecMeta.nFileSamp);
        else
            [raw_LFP tvec SR chan_config sorted_config] = load_LFP_NPX1(options,column);
        end

        if ~isempty(tvec)
            LFP_tvec = tvec;
        else
            LFP_tvec = [];
        end
        tic
        
        % Ripple band
        parameters = list_of_parameters;
        filter_type  = 'bandpass';
        filter_width = [125 300];                 % range of frequencies in Hz you want to filter between
        filter_order = round(6*SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
        norm_freq_range = filter_width/(SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_gamma = fir1(filter_order, norm_freq_range,filter_type);

        for nchannel = 1:size(sorted_config,1)
            LFP.probe(nprobe).ripple(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
            %     LFP(nchannel).gamma_zscore = zscore(abs(hilbert(LFP(nchannel).gamma)));
            LFP.probe(nprobe).ripple_hilbert(nchannel,:) = abs(hilbert(filtfilt(b_gamma,1,raw_LFP(nchannel,:))));
        end
        toc
        tic
        % Slow wave band
        parameters = list_of_parameters;
        filter_type  = 'bandpass';
        filter_width = [0.5 3];                 % range of frequencies in Hz you want to filter between
        filter_order = round(6*SR/(max(filter_width)-min(filter_width)));  % creates filter for ripple
        norm_freq_range = filter_width/(SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
        b_gamma = fir1(filter_order, norm_freq_range,filter_type);

        for nchannel = 1:size(sorted_config,1)
            LFP.probe(nprobe).SO(nchannel,:) = filtfilt(b_gamma,1,raw_LFP(nchannel,:));
            %     LFP(nchannel).gamma_zscore = zscore(abs(hilbert(LFP(nchannel).gamma)));
            LFP.probe(nprobe).SO_hilbert(nchannel,:) = abs(hilbert(filtfilt(b_gamma,1,raw_LFP(nchannel,:))));
        end
        toc

        speed_during_LFP = interp1(position.t,position.v,LFP_tvec,'nearest');
        
%         immobility_period = speed_during_LFP(speed_during_LFP <= 5);
%         running_period = speed_during_LFP(speed_during_LFP > 5);
        
        for nchannel = 1:size(sorted_config,1)
            
            immobility_channel_power.SO{nprobe}(nchannel) = mean(LFP.probe(nprobe).SO_hilbert(nchannel,speed_during_LFP <= 1));
            running_channel_power.SO{nprobe}(nchannel) = mean(LFP.probe(nprobe).SO_hilbert(nchannel,speed_during_LFP > 20));

            immobility_channel_power.ripple{nprobe}(nchannel) = mean(LFP.probe(nprobe).ripple_hilbert(nchannel,speed_during_LFP <= 1));
            running_channel_power.ripple{nprobe}(nchannel) = mean(LFP.probe(nprobe).ripple_hilbert(nchannel,speed_during_LFP > 20));
        end
        figure
        plot(immobility_channel_power.SO{nprobe},sorted_config.Ks_ycoord);hold on
        plot(running_channel_power.SO{nprobe},sorted_config.Ks_ycoord);

        plot([0 1],[chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel) chan_config.Ks_ycoord(best_channels{nprobe}.first_in_brain_channel)],'--k','LineWidth',2)
        plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L4_channel))],'--b','LineWidth',2)
        plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel ==best_channels{nprobe}.L5_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.L5_channel))],'--c','LineWidth',2)

        if ~isempty(best_channels{nprobe}.CA1_channel)
            plot([0 1],[chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel)) chan_config.Ks_ycoord(find(chan_config.Channel == best_channels{nprobe}.CA1_channel))],'--r','LineWidth',2)
        end
        xlim([0 max(immobility_channel_power.SO{nprobe})])
    end
end

%% Manual update best channels (this is manual for now...)
Stimulus_type = 'Checkerboard';
for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
    load best_channels 
    load extracted_PSD
%     openfig('Checkerboard event (all filtered).fig')
    for nprobe = 1:length(session_info.probe)
        options = session_info.probe(nprobe);
        options.importMode = 'LF';
%         options.importMode = 'KS';
        options.probe_no = options.probe_id + 1;
        column = 1;
        [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,column);
       
        [best_channels_updated{nprobe}] = update_best_channels(options,PSD{nprobe},power{nprobe},best_channels{nprobe},sorted_config);
        
        % If updating all 
        best_channels{nprobe}.first_in_brain_channel = best_channels_updated{nprobe}(1);
        best_channels{nprobe}.L4_channel = best_channels_updated{nprobe}(2);
        best_channels{nprobe}.L5_channel = best_channels_updated{nprobe}(3);
        best_channels{nprobe}.CA1_channel = best_channels_updated{nprobe}(4);

        save best_channels best_channels

    end
    

    save best_channels best_channels
end


%% Visual tuning based on SparseNoise
clear all
SUBJECTS = {'M23017','M23028','M23029'};
SUBJECTS = {'M23087'};
experiment_info = subject_session_stimuli_mapping(SUBJECTS);
% Stimulus_type = 'Masa2tracks';
ROOTPATH = 'Z:\ibn-vision';
Stimulus_type = 'SparseNoise_fullscreen';
initMap = [];
probe_hemisphere_text = {'Left','Right'};

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        cd(fullfile(ROOTPATH,'DATA','SUBJECTS',session_info.probe(1).SUBJECT,'ephys',session_info.probe(1).SESSION,'analysis'))
        load('best_channels.mat')

    for nprobe = 1:length(session_info.probe) % For each session, how many probes
        options = session_info.probe(nprobe);
        options.BinWidth = 1/60;
        options.importMode = 'KS'; % LF or MUA or KS
        % options.importMode = 'LF'; % LF or MUA or KS
        options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
        options.stim_dur = 0.1;
        options.AnalysisTimeWindow = [0 1/60*7];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
        options.ks_unitType = 'good'; % 'mua', 'good' or ''
        options.PD_FLAG = 1;
        options.paradigm = 'masa';
        options.gFileNum = gFileNum;
        options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

        [file_to_use imecMeta chan_config sorted_config] = extract_NPX_channel_config(options,1);
        [resps,otherData,stimData,~,wheelData,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
        onDuration = photodiodeData.stim_off.sglxTime-photodiodeData.stim_on.sglxTime;
        offDuration = photodiodeData.stim_on.sglxTime(2:end)-photodiodeData.stim_off.sglxTime(1:end-1);
        % Calculate sparsenoise RFs
        stim_matrix = cat(3,stimData.stim_matrix{:}); % N rows x M cols x nFrames
        
        sn_options.grid_size = [size(stim_matrix,1) size(stim_matrix,2)];
        sn_options.mapSampleRate = 60; % Hz
        sn_options.mapsToShow = {'linear','black','white','contrast'};
        sn_options.mapMethod = 'fitlm'; % fitlm mean
        sn_options.framesToShow = 1:7;  % at 60 Hz would be 8 ms, 40, 72 etc
        sn_options.plotflag = 0; % 1 if you want to see details of the receptive fields
        
        horizontalAngles = -30:(30+90)/8:90;
        verticalAngles = -120:240/16:120;
        % Extract wheel and eye data
        wheel_data = squeeze(nanmean(otherData(1,:,:),2))';
        eye_data = squeeze(nanmean(otherData(2,:,:),2))';
        
        

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

        RF.probe(options.probe_no).cluster_id = options.cluster_id;
        RF.probe(options.probe_no).peak_channel = options.peakChannel;
        RF.probe(options.probe_no).peak_location = chan_config.Ks_ycoord(options.peakChannel);

        RF.probe(options.probe_no).RF_map = initMap;
        initMap = [];

    end
    save(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','receptiveFields.mat'),'RF')


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

% load(fullfile(EPHYS_DATAPATH,'analysis','receptiveFields.mat'))


%% Visual tuning based on Static Gratings
Stimulus_type = 'StaticGratings';
for nsession =10:length(experiment_info)
    session_info = experiment_info(nsession).stimuli_type(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
        gFileNum = experiment_info(nsession).gFileNum(contains(experiment_info(nsession).StimulusName,Stimulus_type));

    for nprobe = 1:length(session_info.probe) % For each session, how many probes
        options = session_info.probe(nprobe);
        options.BinWidth = 1/60;
        options.importMode = 'KS'; % LF or MUA or KS
        % options.importMode = 'LF'; % LF or MUA or KS
        options.BinWidth = 1/60; % resolution (in s) of output resps (e.g. 1/60)
        options.stim_dur = 0.1;
        options.AnalysisTimeWindow = [0 2];% two-element vector specifying time window around stim-on (e.g. [-0.25 1.25])
        options.ks_unitType = 'good'; % 'mua', 'good' or ''
        options.PD_FLAG = 1;
        options.paradigm = 'SG';
        options.gFileNum = gFileNum;
        options.probe_no = options.probe_id+1; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)

        [resps,otherData,stimData,~,wheelData,photodiodeData,timeVector,options] = extractAndCollateNPData(options);
        StaticGratings.probe(options.probe_no).resps=resps;
        StaticGratings.probe(options.probe_no).otherData=otherData;
        StaticGratings.probe(options.probe_no).stimData=stimData;
        StaticGratings.probe(options.probe_no).wheelData=wheelData;
        StaticGratings.probe(options.probe_no).photodiodeData=photodiodeData;
        StaticGratings.probe(options.probe_no).timeVector=timeVector;
        StaticGratings.probe(options.probe_no).options=options;
        switch(stimulus_name{1,1})
            case {'StaticGratings_short','StaticGratings'}
                stim_orientation = readmatrix('X:\ibn-vision\CODE\DEV\BONSAI\Diao\dome_dual_DT\Grating_trials_short.CSV');
            case {'StaticGratings_long'}
                stim_orientation = readmatrix('X:\ibn-vision\CODE\DEV\BONSAI\Diao\dome_dual_DT\Grating_trials.CSV');
        end
        orientation_angles = unique(stim_orientation);
        stim_orientation_tmp = stim_orientation;
        stim_orientation_tmp(wheelData.staticgrating_idx_error) = -1;
        stim_orientation_wo_error = find(stim_orientation_tmp>-1);
        stim_index_wo_error = zeros(200,1);
        stim_index_wo_error(stim_orientation_wo_error) = 1;
        grating_response = cell(length(orientation_angles),1);
        avg_grating_response = zeros([size(resps,[1 2]),length(orientation_angles)]);
        for iAngle = 1:length(orientation_angles)
            orientation_idx = find(stim_orientation == orientation_angles(iAngle));
            grating_response{iAngle,1} = resps(:,:,orientation_idx);
            avg_grating_response(:,:,iAngle) = mean(grating_response{iAngle,1},3);
        end
        grating_response_wo_error = cell(length(orientation_angles),1);
        avg_grating_response_wo_error = zeros([size(resps,[1 2]),length(orientation_angles)]);
        for iAngle = 1:length(orientation_angles)
            orientation_idx = find(stim_orientation == orientation_angles(iAngle));
            iangle_idx = zeros(size(stim_orientation,1),1);
            iangle_idx(orientation_idx) = 1;
            iangle_idx_wo_error = iangle_idx & stim_index_wo_error;
            grating_response_wo_error{iAngle,1} = resps(:,:,iangle_idx_wo_error);
            avg_grating_response_wo_error(:,:,iAngle) = mean(grating_response_wo_error{iAngle,1},3);
        end
        StaticGratings.probe(options.probe_no).stim_orientation=stim_orientation;
        StaticGratings.probe(options.probe_no).stim_index_wo_error=stim_index_wo_error;
        StaticGratings.probe(options.probe_no).grating_response=grating_response;
        StaticGratings.probe(options.probe_no).grating_response_wo_error=grating_response_wo_error;
        tuning_curve = plot_tuning_curve(grating_response,0,0);
        StaticGratings.probe(options.probe_no).tuning_curve=tuning_curve;



    end
    cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
    save('StaticGrating.mat','StaticGratings')
end

%% Static Gratings Tuning Curve Plotting
cd(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis'))
load('StaticGratings.mat')
tuning_curve = plot_tuning_curve(StaticGratings.probe(1).grating_response,1,0);