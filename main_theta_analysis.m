%% Main theta analysis codes

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\adam.tong\Documents\GitHub\VR_NPX_analysis'))

%% Theta modulation and phase percession

clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = 10
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
        %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'power');
        load(fullfile(options.ANALYSIS_DATAPATH,'..','extracted_PSD.mat'),'power');

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

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end


        % Phase change
        %         PhDiff = phdiffmeasure(x, y, fs, method)

        theta_modulation_L = [];
        theta_modulation_R = [];

        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);
            probe_no = session_info(n).probe(nprobe).probe_id + 1;
            options.probe_no = probe_no; % probe_no is [1,2] it is redundant as we have options.probe_id (0 and 1)
            %                 Behavioural state detection
%             speed_interp = interp1(Behaviour.tvec,Behaviour.speed,LFP(nprobe).tvec','linear');
%             speedTreshold = 1;

            if isfield(LFP(nprobe),'L4')
                cortex_LFP = LFP(nprobe).L4;
            elseif isfield(LFP(nprobe),'L5')
                cortex_LFP = LFP(nprobe).L5;
            elseif isfield(LFP(nprobe),'MEC')

            end

            if isfield(LFP(nprobe),'CA1')
                CA1_LFP = LFP(nprobe).CA1;
            end


            if options(nprobe).probe_hemisphere == 1
                [theta_modulation_L]= phase_precession_absolute_location(LFP(nprobe).tvec,CA1_LFP,place_fields,clusters_combined,Task_info,Behaviour,options);
            elseif options(nprobe).probe_hemisphere == 2
                [theta_modulation_R]= phase_precession_absolute_location(LFP(nprobe).tvec,CA1_LFP,place_fields,clusters_combined,Task_info,Behaviour,options);
            end
        end

%       plot_perievent_spiketimes  
        plot_theta_modulation(theta_modulation_L,theta_modulation_R);
            
    end
end



%% Decoding error and log odds for V1 and HPC
lap_times = [];
clear all
SUBJECTS = {'M23017','M23028','M23029','M23087','M23153'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
Stimulus_type = 'RUN';
% [1 2 3 4 9 10 12 14]

for nsession = [9 10 12 14]
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

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
        load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

        if length(clusters) > 1
            clusters_combined = combine_clusters_from_multiple_probes(merged_clusters(1),merged_clusters(2));
        else
            clusters_combined = merged_clusters;
        end

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_HPC.mat'),'estimated_position_lap_CV_HPC_combined','estimated_position_lap_CV_HPC','estimated_position_lap_CV_shuffled_HPC','estimated_position_lap_CV_shuffled_HPC_combined')
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_HPC.mat'),'probability_ratio_RUN_lap_HPC_combined','probability_ratio_RUN_lap_HPC','estimated_position_lap_CV_shuffled_V1','estimated_position_lap_CV_shuffled_V1_combined')

        load(fullfile(options.ANALYSIS_DATAPATH,'estimated_position_lap_CV_V1.mat'),'estimated_position_lap_CV_V1',"estimated_position_lap_CV_V1_combined")
        load(fullfile(options.ANALYSIS_DATAPATH,'probability_ratio_RUN_lap_V1.mat'),'probability_ratio_RUN_lap_V1','probability_ratio_RUN_lap_V1_combined')
        

        % plotting decoded run trajectory
        for nprobe = 1:length(clusters)
            options = session_info(n).probe(nprobe);

            if exist(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))== 0
                mkdir(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'))
            end

            if clusters(nprobe).probe_hemisphere == 1
                options.region = 'HPC Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                options.region = 'HPC Right';
            end

            plot_decoding_RUN_trajectory(estimated_position_lap_CV_HPC(nprobe).track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])

            if clusters(nprobe).probe_hemisphere == 1
                metric_param.region = @(x) contains(x,'V1_L');
                options.region = 'V1 Left';
            elseif clusters(nprobe).probe_hemisphere == 2
                metric_param.region = @(x) contains(x,'V1_R');
                options.region = 'V1 Right';
            end

            plot_decoding_RUN_trajectory(estimated_position_lap_CV_V1(nprobe).track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
        end

        if length(session_info(n).probe) > 1
            options.region = 'HPC Combined';
            plot_decoding_RUN_trajectory(estimated_position_lap_CV_HPC_combined.track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])

            options.region = 'V1 Combined';
            plot_decoding_RUN_trajectory(estimated_position_lap_CV_V1_combined.track,Task_info,options)
            save_all_figures(fullfile(options.ANALYSIS_DATAPATH,'..','figures','RUN Bayesian decoding'),[])
        end


        % initialise all variables
        for track_id = 1:length(place_fields)
            for temp_track = 1:length(place_fields)
                actual_position{nsession}{track_id} = [];
                actual_speed{nsession}{track_id} = [];
                VR_speed{nsession}{track_id} = [];
                decoded_position_lap_id{nsession}{track_id}  = [];


                decoded_position_HPC_combined{nsession}{track_id} = [];
                decoded_error_HPC_combined{nsession}{track_id}{temp_track}  = [];

                decoded_position_HPC_combined_shuffled{nsession}{track_id} = [];
                decoded_error_HPC_combined_shuffled{nsession}{track_id}{temp_track}  = [];

                decoded_position_V1_combined{nsession}{track_id} = [];
                decoded_error_V1_combined{nsession}{track_id}{temp_track}  = [];

                for nprobe = 1:length(session_info(n).probe)
                    probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;
                    decoded_position_V1{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                    decoded_position_HPC{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} = [];

                    decoded_position_V1_shuffled{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                    decoded_position_HPC_shuffled{probe_hemisphere}{nsession}{track_id} = [];
                    decoded_error_HPC_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                end
            end
        end

        % Position bins from both tracks
        position_bin_across_tracks = estimated_position_lap_CV_V1(1).track(1).lap(1).track(1).position_bin_centres;
        position_bin_across_tracks = [position_bin_across_tracks position_bin_across_tracks+1000];

        for track_id = 1:length(place_fields)
            for temp_track = 1:length(place_fields)
                for lap_id = 1:length(estimated_position_lap_CV_V1(nprobe).track(track_id).lap)

                    if temp_track == 1 %Do not need to loop twice
                        decoded_position_lap_id{nsession}{track_id} = [decoded_position_lap_id{nsession}{track_id} ...
                            lap_id*ones(1,length(estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).run_actual_position))];
                        actual_position{nsession}{track_id} = [actual_position{nsession}{track_id} ...
                            estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).run_actual_position];
                        actual_speed{nsession}{track_id} = [actual_speed{nsession}{track_id} ...
                            estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).actual_run_speed];
                        VR_speed{nsession}{track_id} = [VR_speed{nsession}{track_id} ...
                            estimated_position_lap_CV_V1(1).track(track_id).lap(lap_id).track(track_id).run_speed];
                    end

%                     if length(session_info(n).probe)==1
%                         decoded_error_HPC{nsession}{track_id}{temp_track} = [decoded_error_HPC{nsession}{track_id}{temp_track} ...
%                             estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(temp_track).peak_position...
%                             - estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(track_id).run_actual_position];
%                         if temp_track == 1 %Do not need to loop twice
%                             [~,index] = max([estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(1).run; ...
%                                 estimated_position_lap_CV_HPC.track(track_id).lap(lap_id).track(2).run]);
%                             decoded_position_HPC{nsession}{track_id} = [decoded_position_HPC{nsession}{track_id} index];
%                         end
%                     end

                    for nprobe = 1:length(session_info(n).probe)
                        probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

                        decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        [~,index] = max(estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(temp_track).run);
                        decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_HPC{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_V1(1).track(1).lap(1).track(1).position_bin_centres(index)...
                            - estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        [~,index] = max(estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(temp_track).run);
                        decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_shuffled_V1(1).track(1).lap(1).track(1).position_bin_centres(index)...
                            - estimated_position_lap_CV_shuffled_V1(nprobe).track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        decoded_position_V1_shuffled{probe_hemisphere}{nsession}{track_id} = [];
                        decoded_error_V1_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];
                        decoded_position_HPC_shuffled{probe_hemisphere}{nsession}{track_id} = [];
                        decoded_error_HPC_shuffled{probe_hemisphere}{nsession}{track_id}{temp_track} = [];

                        if temp_track == 1 %Do not need to loop twice
                            [~,index] = max([estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_V1(nprobe).track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_V1{probe_hemisphere}{nsession}{track_id} = [decoded_position_V1{probe_hemisphere}{nsession}{track_id} index];

                            [~,index] = max([estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(1).run; ...
                                estimated_position_lap_CV_HPC(nprobe).track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_HPC{probe_hemisphere}{nsession}{track_id} = [decoded_position_HPC{probe_hemisphere}{nsession}{track_id} index];


                            [~,index] = max([estimated_position_lap_CV_shuffled_V1(nprobe).track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_shuffled_V1(nprobe).track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_V1{probe_hemisphere}{nsession}{track_id} = [decoded_position_V1{probe_hemisphere}{nsession}{track_id} index];

                            [~,index] = max([estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(1).run; ...
                                estimated_position_lap_CV_shuffled_HPC(nprobe).track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_HPC{probe_hemisphere}{nsession}{track_id} = [decoded_position_HPC{probe_hemisphere}{nsession}{track_id} index];
                        end
                    end

                    if length(session_info(n).probe)>1
                        decoded_error_HPC_combined{nsession}{track_id}{temp_track} = [decoded_error_HPC_combined{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        decoded_error_V1_combined{nsession}{track_id}{temp_track} = [decoded_error_V1_combined{nsession}{track_id}{temp_track} ...
                            estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(temp_track).peak_position...
                            - estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(track_id).run_actual_position];

                        if temp_track == 1 %Do not need to loop twice
                            [~,index] = max([estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_HPC_combined.track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_HPC_combined{nsession}{track_id} = [decoded_position_HPC_combined{nsession}{track_id} index];

                            [~,index] = max([estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(1).run;...
                                estimated_position_lap_CV_V1_combined.track(track_id).lap(lap_id).track(2).run]);
                            decoded_position_V1_combined{nsession}{track_id} = [decoded_position_V1_combined{nsession}{track_id} index];
                        end

                    end
                end
            end
        end


        % Plotting decoding performance
        if length(session_info(n).probe)>1
            decoding_performance = plot_within_session_decoded_error_HPC_V1(decoded_position_lap_id,VR_speed,actual_position,estimated_position_lap_CV_V1,decoded_position_V1,decoded_position_HPC,decoded_position_HPC_combined,decoded_position_V1_combined,...
                decoded_error_V1,decoded_error_HPC,decoded_error_HPC_combined,decoded_error_V1_combined,place_fields,session_info(n),nsession)   
        else
          decoding_performance = plot_within_session_decoded_error_HPC_V1(VR_speed,actual_position,estimated_position_lap_CV_V1,decoded_position_V1,decoded_position_HPC,[],[],...
                decoded_error_V1,decoded_error_HPC,[],[],place_fields,session_info(n),nsession)   
        end
       
        save_all_figures(fullfile(ROOTPATH,'DATA','SUBJECTS',options.SUBJECT,'ephys',options.SESSION,'analysis','spatial cells'),[])

    end
end




    decoded_position_V1{probe_hemisphere}{nsession}{track_id};

    scatter(decoded_position_HPC{nsession}{2}(VR_speed{nsession}{2}>5)...
        ,decoded_position_V1{1}{nsession}{2}(VR_speed{nsession}{2}>5),'blue','filled','MarkerFaceAlpha',0.05)

    [N,Xedges,Yedges,binX,binY] = histcounts2(decoded_position_HPC{nsession}{2}(VR_speed{nsession}{2}>5),decoded_position_V1{1}{nsession}{2}(VR_speed{nsession}{2}>5),14);
    imagesc(flip(N)/max(max(N)))
   
