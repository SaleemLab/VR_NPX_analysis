%% Main DLAG analysis

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\mDLAG\mDLAG'))
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\mDLAG\DLAG'))
% addpath core_mdlag
% addpath plotting
% addpath simulation
% addpath util
% addpath variable_transformations
% addpath descriptive_statistics
% addpath demo
cd('C:\Users\masahiro.takigawa\Documents\GitHub\mDLAG\mDLAG\demo')

%% mDLAG Populational level V1-HPC interaction during context-selective ripples 

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]
% load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'),'place_fields_all_L','place_fields_all_R','place_fields_all_combined')

psthBinSize = 0.02; %create 1ms bin
for nsession = [1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'));

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = merged_clusters;
            sorting_option = 'spikeinterface';
        elseif exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = clusters_ks3;
            sorting_option = 'spikeinterface';
        else
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            sorting_option = 'old';
        end

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_candidate_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         tic
%         load(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'))
%         toc

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        ripple_modulation_L = [];
        ripple_modulation_R = [];
        ripple_modulation_combined = [];
        load(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'))

        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;
        %             metric_param.region = @(x) contains(x,'V1');
        metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id);
        %             metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id(ia));
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        [C,ia,ic] = unique(selected_clusters.merged_cluster_id);

        if length(merged_clusters)>1 % combine ripples (keep ripple hemisphere and event context ID)
            event_id = [ones(1,length(ripples(1).T1_onset)) ones(1,length(ripples(2).T1_onset)) 2*ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
            [event_times,index] = sort([ripples(1).T1_onset ripples(2).T1_onset ripples(1).T2_onset ripples(2).T2_onset]);
            probe_hemisphere = [session_info(n).probe.probe_hemisphere];
            temp = [ones(1,length(ripples(1).T1_onset)) 2*ones(1,length(ripples(2).T1_onset)) ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
        else
            event_id = [ones(1,length(ripples(1).T1_onset)) 2*ones(1,length(ripples(1).T2_onset))];
            [event_times,index] = sort([ripples(1).T1_onset ripples(2).T1_onset ripples(1).T2_onset]);
            probe_hemisphere = [session_info(n).probe.probe_hemisphere];
            temp = ones(1,length(event_id));
        end

        probe_id = temp;
        for nprobe = 1:length(merged_clusters)
            probe_id(temp==nprobe) = probe_hemisphere(nprobe);
        end
        clear temp

        %         [event_times,index] = sort([Task_info.start_time_all]);
        %         event_id = Task_info.track_ID_all;

        event_id = event_id(index); % context id
        probe_id = probe_id(index); % probe id 

        % whole session timevec variables
        timevec = Behaviour.tvec';
        timevec_edge = (timevec(1)-(psthBinSize)/2....
            :psthBinSize:...
            timevec(end)+(psthBinSize)/2)';

        all_spike_counts = [];
        
        %         event_id = [ones(1,size(ripple_modulation_combined(1).spike_count{ncell},1))  2*ones(1,size(ripple_modulation_combined(2).spike_count{ncell},1))];
        %                 ripple_cells = unique([find(ripple_modulation_combined(1).ripple_modulation_percentile>=0.95) find(ripple_modulation_combined(2).ripple_modulation_percentile>=0.95)]);
        ripple_cells = 1:length(ripple_modulation_combined(1).ripple_modulation_percentile);

        for ncell = 1:length(ripple_cells)

            spike_times = clusters_combined.spike_times(clusters_combined.merged_spike_id == ripple_modulation_combined(1).cluster_id(ripple_cells(ncell)));

            %             % Convolve spike count time series with Gaussian window
            %             y = histcounts(spike_times, timevec_edge)';
            %             y = conv(y, gaussianWindow, 'same')/psthBinSize; % all spikes

            spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');
            spike_times = spike_times(spike_speed<5);

            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times,...
                event_times, [-0.1 0.5], psthBinSize);
            all_spike_counts(ncell,:,:) = binnedArray'; % ncell X time x events
        end

        eventSpikes = [];
        yDims = [];
        if length(session_info(n).probe)>1
            yDims(1) = sum(contains(cell_regions,'V1_L'));
            yDims(2) = sum(contains(cell_regions,'V1_R'));
            yDims(3) = sum(contains(cell_regions,'HPC'));
            selected_cell_id = [find(contains(cell_regions,'V1_L')) find(contains(cell_regions,'V1_R'))...
                find(contains(cell_regions,'HPC'))];
        elseif session_info(n).probe(1).probe_hemisphere==1
            yDims(1) = sum(contains(cell_regions,'V1_L'));
%             yDims(2) = sum(contains(cell_regions,'V1_R'));
            yDims(2) = sum(contains(cell_regions,'HPC'));
            selected_cell_id = [find(contains(cell_regions,'V1_L'))...
                find(contains(cell_regions,'HPC'))];
        elseif session_info(n).probe(1).probe_hemisphere==2
%             yDims(1) = sum(contains(cell_regions,'V1_L'));
            yDims(1) = sum(contains(cell_regions,'V1_R'));
            yDims(2) = sum(contains(cell_regions,'HPC'));
            selected_cell_id = [find(contains(cell_regions,'V1_R'))...
                find(contains(cell_regions,'HPC'))];
        end

        cell_regions = ripple_modulation_combined(1).region(ripple_cells);
        for nevent = 1:size(all_spike_counts,3)
            eventSpikes(nevent).trialId = nevent;% event id
            eventSpikes(nevent).probeID = probe_id(nevent);% event originated from probe 1 or 2 (1 = L and 2 = R)
            eventSpikes(nevent).trackID = event_id(nevent);% Track L and Track R event  (1 = L and 2 = R)
            eventSpikes(nevent).reactivationID = event_id(nevent);% Track L and Track R event  (1 = L and 2 = R)
            eventSpikes(nevent).yDims = yDims; % just to save region info. 
            % Number of neurons per region: e.g. yDims(1)= 20 and yDims(2)=
            % 30 mean region 1 there are 20 neruons and region 2 there are
            % 30 neruons.
            eventSpikes(nevent).T = size(all_spike_counts,2); % Number of timesteps
            eventSpikes(nevent).y= [] ;
            eventSpikes(nevent).y = all_spike_counts(selected_cell_id,:,nevent);% selected spike counts this event
        end

        % all ripple event id
        options.load_fitted_model=[];
        [estParams,trackedParams,flags] = mDLAG_analysis_pipeline(eventSpikes,psthBinSize,yDims,options);
    end
end














%% CCA Populational level V1-HPC interaction during context-selective ripples 

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]
load(fullfile('D:\corticohippocampal_replay\summary','place_fields_all.mat'),'place_fields_all_L','place_fields_all_R','place_fields_all_combined')

psthBinSize = 0.001; %create 1ms bin
% Define Gaussian window for smoothing
gaussianWindow = gausswin(0.2*1/psthBinSize);
% Normalize to have an area of 1 (i.e., to be a probability distribution)
gaussianWindow = gaussianWindow / sum(gaussianWindow);


for nsession = [1 2 3 4 6 7 8 9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'));

        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP')
        end

        if exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('merged_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = merged_clusters;
            sorting_option = 'spikeinterface';
        elseif exist(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))))
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks3%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters = clusters_ks3;
            sorting_option = 'spikeinterface';
        else
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            sorting_option = 'old';
        end

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_place_fields.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         tic
%         load(fullfile(options.ANALYSIS_DATAPATH,'ripple_modulation.mat'))
%         toc

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        ripple_modulation_L = [];
        ripple_modulation_R = [];
        ripple_modulation_combined = [];

        metric_param = create_cluster_selection_params('sorting_option',sorting_option);
        metric_param.unstable_ids = @(x) x==0;
        %             metric_param.region = @(x) contains(x,'V1');
        metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id);
        %             metric_param.merged_cluster_id = @(x) ismember(x,clusters_combined.merged_cluster_id(ia));
        [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
        [C,ia,ic] = unique(selected_clusters.merged_cluster_id);

        event_id = [ones(1,length(ripples(1).T1_onset)) ones(1,length(ripples(2).T1_onset)) 2*ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
        
        if length(merged_clusters)>1
            probe_hemisphere = [session_info(n).probe.probe_hemisphere];
            temp = [ones(1,length(ripples(1).T1_onset)) 2*ones(1,length(ripples(2).T1_onset)) ones(1,length(ripples(1).T2_onset)) 2*ones(1,length(ripples(2).T2_onset))];
        else
            probe_hemisphere = [session_info(n).probe.probe_hemisphere];
            temp = ones(1,length(event_id));
        end

        % whole session timevec variables
        timevec = Behaviour.tvec';
        timevec_edge = (timevec(1)-(psthBinSize)/2....
            :psthBinSize:...
            timevec(end)+(psthBinSize)/2)';

        probe_id = temp;
        for nprobe = 1:length(merged_clusters)
            probe_id(temp==nprobe) = probe_hemisphere(nprobe);
        end

%         [event_times,index] = sort([ripples(1).T1_onset ripples(2).T1_onset ripples(1).T2_onset ripples(2).T2_onset]);
        [event_times,index] = sort([Task_info.start_time_all]);
        event_id = Task_info.track_ID_all;

        event_id = event_id(index);
        probe_id = probe_id(index);

        all_spike_counts = [];
        %         event_id = [ones(1,size(ripple_modulation_combined(1).spike_count{ncell},1))  2*ones(1,size(ripple_modulation_combined(2).spike_count{ncell},1))];
%                 ripple_cells = unique([find(ripple_modulation_combined(1).ripple_modulation_percentile>=0.95) find(ripple_modulation_combined(2).ripple_modulation_percentile>=0.95)]);
        ripple_cells = 1:length(ripple_modulation_combined(1).ripple_modulation_percentile);

        for ncell = 1:length(ripple_cells)

            spike_times = clusters_combined.spike_times(clusters_combined.merged_spike_id == ripple_modulation_combined(1).cluster_id(ripple_cells(ncell)));

            %             % Convolve spike count time series with Gaussian window
            %             y = histcounts(spike_times, timevec_edge)';
            %             y = conv(y, gaussianWindow, 'same')/psthBinSize; % all spikes

            spike_speed = interp1(Behaviour.tvec,Behaviour.speed,spike_times,'nearest');
            spike_times = spike_times(spike_speed<5);

            [psth, bins, rasterX, rasterY, spikeCounts, binnedArray] = psthAndBA(spike_times,...
                event_times, [-2 2], psthBinSize);
            all_spike_counts(ncell,:,:) = binnedArray(:,bins<0.5&bins>-0.5)'; % ncell X time x events

%             if size(binnedArray,1)==1
%                 psth = conv(binnedArray/psthBinSize,gaussianWindow,'same');
%             else
%                 for nevent = 1:size(binnedArray,1)
%                     psth(nevent,:) = conv(binnedArray(nevent,:)/psthBinSize,gaussianWindow,'same');
%                 end
                %                 psth = (psth-mean(y))/std(y);
                %                 psth_track1 = mean(psth_track1,'omitnan');
%             end

%             all_spike_counts(ncell,:,:) = psth(:,bins<0.5&bins>-0.5)'; % ncell X time x events
            %                 all_spike_counts{track_id}(ncell,:,:) = [ripple_modulation_combined(1).spike_count{ncell}(:,bins<0.5&bins>-0.5);...
            %                     ripple_modulation_combined(2).spike_count{ncell}(:,bins<0.5&bins>-0.5);]; % ncell X time x events
        end

        % V1 L and R
        ripple_spike_count = [];
        cell_regions = ripple_modulation_combined(1).region(ripple_cells);
        for nprobe = 1:length(probe_hemisphere)

            if session_info(n).probe(nprobe).probe_hemisphere == 1
                ripple_spike_count{1} = all_spike_counts(contains(cell_regions,'V1_L'),:,:);
            elseif session_info(n).probe(nprobe).probe_hemisphere == 2
                ripple_spike_count{2} = all_spike_counts(contains(cell_regions,'V1_R'),:,:);
            end

        end

        % HPC
        ripple_spike_count{3} = all_spike_counts(contains(cell_regions,'HPC'),:,:);


        %%% Example computation of a cross-correlation map
        % Figs. 3 and 4

        % The units of the arguments are with respect to the binning window used
        % to bin spikes.
        pairs = [1 2;3 1;3 2];
        pair_name = {'V1 L - V1 R','V1 L - HPC','V1 R - HPC'};
        argOut_all = [];
        for npair = 1:size(pairs,1)
            argIn=[];
            argIn.BinWidth = 2;     % 20ms
            argIn.MaxDelay = 50;    % 500ms
            argIn.TimeStep = 4;     % 40ms
            argIn.WindowLength = 8; % 80ms

            %argIn.NumWorkers = Inf; % Requires Parallel Processing Toolbox

            disp(argIn)
                        argOut = ComputeCorrMap({ripple_spike_count{pairs(npair,1)} ripple_spike_count{pairs(npair,2)}}, event_id, argIn);
%             argOut = ComputeCorrMap({ripple_spike_count{2}(2:2:size(ripple_spike_count{2},1),:,:) ripple_spike_count{2}(1:2:size(ripple_spike_count{2},1),:,:)}, event_id, argIn);

            %%%
            CANONICAL_PAIR_IDX = 1;
            mapDim = size(argOut.CorrMap, 2);
            delays = (-argIn.MaxDelay:argIn.MaxDelay)*10; % Convert to ms
            t = (-0.5*100:argIn.TimeStep:0.5*100)*10; % Convert to ms

            figure(npair);
            subplot(2,2,1)

            imagesc( delays, t, argOut.CorrMap(:,:,CANONICAL_PAIR_IDX)' )

            ax = gca;
            ax.YDir = 'Normal';

            xlabel('Delay')
            ylabel('Time')
            colorbar
            clim([0.1 0.5])
            xlim([-200 200])
            ylim([-500 500])
%             colormap(flip(gray))
            title(pair_name{npair})

            subplot(2,2,2)

            imagesc( delays, t, argOut.CorrMap(:,:,2)' )

            ax = gca;
            ax.YDir = 'Normal';

            xlabel('Delay')
            ylabel('Time')
            colorbar
            clim([0.1 0.5])
            xlim([-200 200])
            ylim([-500 500])
%             colormap(flip(gray))
            title('2nd axis')

            subplot(2,2,3)

            imagesc( delays, t, argOut.FrMap(:,:)' )

            ax = gca;
            ax.YDir = 'Normal';

            xlabel('Delay')
            ylabel('Time')
            colorbar
%             clim([0.2 0.4])
            xlim([-200 200])
            ylim([-500 500])
            %             colormap(flip(gray))
            title('geomertic mean FR')

            subplot(2,2,4)

            plot( delays,mean(argOut.CorrMap(:,t>=-40 & t<=40,CANONICAL_PAIR_IDX),2),'k')
            hold on
            plot( delays,mean(argOut.CorrMap(:,t<=150 & t>=50,CANONICAL_PAIR_IDX),2) ,'r')
            plot( delays,mean(argOut.CorrMap(:,t>=-150 & t<=-50,CANONICAL_PAIR_IDX),2),'b' )
            argOut_all{npair} = argOut;
        end
        %%% Example computation of the interaction structure analysis
        % Fig. 6

        clear argIn

        argIn.TimePeriods = [...
            (  0:20:40)' ( 20:20:60)'; ...
            (128:20:168)' (148:20:188)'] + 5;

        % Can take up to 15min due to the 10-fold cross-validation
        argOut = CovStabilityAcrossTimeAnalysis(spikes, expCond, argIn);

        %%%
        RANK_TO_PLOT = 2;

        normFactor = diag(argOut.CvR(:,:,RANK_TO_PLOT));
        numTimePeriods = size(argIn.TimePeriods, 1);

        figure(2);

        imagesc(argOut.CvR(:,:,RANK_TO_PLOT)./repmat(normFactor', numTimePeriods, 1))

        ax = gca;
        ax.YDir = 'Normal';

        axis square

        xlabel('Time Used For Correlation')
        ylabel('Time Used For Fitting')


    end
end
