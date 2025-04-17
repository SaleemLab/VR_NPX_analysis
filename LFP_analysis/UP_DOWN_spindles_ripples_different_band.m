%% UP DOWN spindles and ripples with different bands
% pyversion('C:\Users\masahiro.takigawa\.conda\envs\fooof\python')
% pyversion('C:\Users\masah\anaconda3\envs\fooof\python')

addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))
pyenv("ExecutionMode","OutOfProcess")

clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar
% experiment_info=experiment_info([4 5 6 ]);
% experiment_info=experiment_info([4 5 6 18 19 21 34 35 44 45 58 59 60 71]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'Sleep';

for nsession =1:15
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));

    for n = 1:length(session_info) % How many recording sessions for spatial tasks (PRE, RUN and POST)
        options = session_info(n).probe(1);

        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
        if isempty(DIR)
            continue
        end

        DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

        if ~isempty(DIR)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
            session_clusters_RUN=session_clusters;
            clear session_clusters
        end

        if ~isempty(DIR1)
            load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
            session_clusters_RUN=session_clusters;
            clear session_clusters
        end
        tic
        if contains(stimulus_name{n},'Masa2tracks')
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));

            load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
            clusters=clusters_ks4;
        elseif contains(stimulus_name{n},'Sleep')
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
            %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
        else
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
            % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

            load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
            clusters=clusters_ks4;
        end
        toc


        for nprobe = 1:length(slow_waves)
            if isfield(slow_waves(nprobe),'ints')
                slow_waves(nprobe).UP_ints = slow_waves(nprobe).ints.UP;
                slow_waves(nprobe).DOWN_ints = slow_waves(nprobe).ints.DOWN;
            end
        end

        if isfield(slow_waves(1),'ints')
            slow_waves = rmfield(slow_waves,'ints');
        end

        params = create_cluster_selection_params('sorting_option','masa');
        % PSD slope quantification using fooof ()
        tvec = LFP(1).tvec;
        SR = round(1/mean(diff(tvec)));
        nfft_seconds= 2;
        nfft = 2^(nextpow2(SR*nfft_seconds));
        win  = hanning(nfft);

        clipDur = 5; % seconds
        timebin_edges = tvec(1):clipDur:tvec(end); % 5 seconds timebin edges for PSD slope
        nClipSamps = round(SR*clipDur);


        disp('cortical wave analysis started')
        tic
        %%%%%%%%%%%% Cortical wave direction during DOWN state peak
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        % https://www.nature.com/articles/s41467-019-10327-5 levenstein et
        % al used 0.5 to 8 Hz for slow wave UP DOWN detection
        filterparms.deltafilter = [9 17];%spindle band
        filterparms.spindlesfilter = [20 50];%heuristically defined.  room for improvement here.
        % filterparms.gammafilter = [100 400];
        % filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
        % filterparms.gammanormwin = 20; %window for gamma normalization (s)

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
                time_idx = 1:length(LFP(1).best_HPC);
            elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

                time_idx = 1:length(LFP(2).best_HPC);
            else
                time_idx = 1:length(LFP(1).best_HPC) ;
            end

            lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
            lfp.timestamps = LFP(probe_no).tvec(time_idx);
            tvec = LFP(probe_no).tvec(time_idx);
            lfp.data=[];
            DOWN_peaks_shank = [];
            peaks_latency = [];
            DOWN_traveling = [];
            DOWN_peaks_zscore= [];

            slow_waves(probe_no).DOWN_peaks_zscore = [];
            slow_waves(probe_no).DOWN_peaktimes = [];
            % slow_waves(probe_no).DOWN_peaks_latency = [];
            % slow_waves(probe_no).DOWN_traveling = [];
            slow_waves(probe_no).probe_hemisphere = [];
            slow_waves(probe_no).shank_id = [];

            if isempty(slow_waves(probe_no).timestamps)
                continue
            end

            if size(slow_waves(probe_no).timestamps,1) < 10
                continue
            end

            if isfield(LFP(probe_no),'average_V1_xcoord') & ~isempty(behavioural_state_merged.SWS)% if exist best V1 channel for sleep
                probe_id = [];

                if length(LFP)==1
                    lfp.data= [LFP(probe_no).average_V1(:,time_idx)'];
                    probe_hemisphere = probe_no*ones(1,length(LFP(probe_no).average_V1_shank_id));
                    slow_waves(probe_no).shank_id = [LFP(probe_no).average_V1_shank_id];
                else
                    lfp.data= [LFP(1).average_V1(:,time_idx)' LFP(2).average_V1(:,time_idx)'];
                    probe_hemisphere = [ones(1,length(LFP(1).average_V1_shank_id)) 2*ones(1,length(LFP(2).average_V1_shank_id))];
                    if size(LFP(probe_no).average_V1_shank_id,1)==1
                        slow_waves(probe_no).shank_id = [LFP(1).average_V1_shank_id LFP(2).average_V1_shank_id];
                    else
                        slow_waves(probe_no).shank_id = [LFP(1).average_V1_shank_id' LFP(2).average_V1_shank_id'];
                    end
                end

                slow_waves(probe_no).probe_hemisphere = probe_hemisphere;


                if nprobe == 1 % Only need to grab once
                    % grab delta LFP
                    deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
                    zscored_LFP = zscore(deltaLFP.data);
                    SO_phase_LFP = deltaLFP.phase;
                    SO_amplitude_LFP = zscore(deltaLFP.amp);
                    deltaLFP = [];

                    % grab spindles LFP
                    filter_type  = 'bandpass';
                    passband = [20 50];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_spindle = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_spindle,1, lfp.data(:,nShank));
                    end
                    spindle_amplitude_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
                    spindle_amplitude_LFP = smoothdata(spindle_amplitude_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
                    spindle_phase_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];

                    % grab ripples LFP
                    filter_type  = 'bandpass';
                    passband = [300 600];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_ripple = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_ripple,1, lfp.data(:,nShank));
                    end
                    ripple_amplitude_cortex_LFP = zscore(abs(hilbert(signal))); % z scored amplitude
                    ripple_phase_cortex_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];
                end

                %%%%% Calculate Delta peak amplitude and timestamp per event
                for nevent = 1:size(slow_waves(probe_no).DOWN_ints,1)
                    % for nevent = 600:640
                    midpoint = mean(slow_waves(probe_no).timestamps(nevent));
                    % midpoint = mean(slow_waves(probe_no).DOWN_ints(nevent,:));

                    tidx = FindInInterval(tvec,[midpoint-0.1 midpoint+0.1]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(midpoint-tvec));

                    for nShank=1:length(probe_hemisphere) % across shanks from both probes if using two probes

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        [~,temp]=min(abs(midpoint-tvec(tidx(peak_id))));
                        if ~isempty(temp)
                            DOWN_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            DOWN_peaks_zscore(nShank,nevent) = SO_amplitude_LFP(tidx(peak_id(temp)),nShank);
                        else
                            DOWN_peaks_shank(nShank,nevent) = nan;
                            DOWN_peaks_zscore(nShank,nevent) = nan;
                        end
                    end

                    % (diff(ordered_xcoord)/1000000)'./diff(DOWN_peaks_shank(:,nevent))
                    %
                    % diff(DOWN_peaks_shank(:,425))

                    % Putative calculation of average latency across shanks
                    % for each hemisphere
                    % for h = unique(probe_hemisphere)
                    %     DOWN_peaks_shank_temp = DOWN_peaks_shank(probe_hemisphere==h,nevent);
                    %
                    %     if sum(~isnan(DOWN_peaks_shank_temp))==4 % if delta peaks on four shanks
                    %         peaks_latency(h,nevent) = mean(diff(DOWN_peaks_shank_temp),'omitnan');
                    %     elseif sum(~isnan(DOWN_peaks_shank_temp))==3 % if delta peaks on three shanks
                    %         skipped_shank= diff(LFP(h).average_V1_shank_id(~isnan(DOWN_peaks_shank_temp)))>1;
                    %         if sum(skipped_shank)>0 % if delta peak skipped one shank
                    %
                    %             if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                    %                 peaks_latency(h,nevent) = DOWN_peaks_shank_temp(1)-DOWN_peaks_shank_temp(2);
                    %             elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                    %                 peaks_latency(h,nevent) = DOWN_peaks_shank_temp(2)-DOWN_peaks_shank_temp(3);
                    %             end
                    %         else
                    %             peaks_latency(h,nevent) = mean(diff(DOWN_peaks_shank_temp),'omitnan');
                    %         end
                    %     elseif sum(~isnan(DOWN_peaks_shank_temp))==2 % if delta peaks on two shanks (latency divide by number of shanks skipped to caluclate mean latency per 250 micron)
                    %         peaks_latency(h,nevent) = diff(DOWN_peaks_shank_temp(1:2))/abs(LFP(h).average_V1_shank_id(2)-LFP(h).average_V1_shank_id(1));
                    %     else % only one peak. Can't calculate latency
                    %         peaks_latency(h,nevent) =nan;
                    %     end
                    %
                    %     % From https://www.jneurosci.org/content/34/26/8875#sec-2
                    %     % speed is ~40 milimeter per seconds
                    %     % 6.25 miliseconds to travel 250 micrometer (rough shank spacing)
                    %     % putatively set the minimum delay threshold to be 3
                    %     % miliseconds.
                    %     if h==1
                    %
                    %         % if direction of mean latency is consistent with the latency more than half
                    %         % of the shank latency
                    %         % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                    %         % (if three shanks, then still 2 jumps needed)
                    %         if peaks_latency(nevent)>0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)>0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = 1; % anterior to posterior
                    %         elseif peaks_latency(nevent)<-0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)<0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = -1; % posterior to anterior
                    %         else
                    %             DOWN_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %         % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                    %     elseif h==2
                    %         if peaks_latency(h,nevent)>0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)>0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = -1; % posterior to anterior
                    %         elseif peaks_latency(h,nevent)<-0.003 & sum(DOWN_peaks_shank_temp(2:end)-DOWN_peaks_shank_temp(1)<0)>=length(LFP(h).average_V1_shank_id)/2
                    %             DOWN_traveling(h,nevent) = 1; % anterior to posterior
                    %         else
                    %             DOWN_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %     end
                    % end
                    % nexttile
                    % hold on;xline(tvec(idx));
                    % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % xline(DOWN_peaks_shank(:,nevent)','r')
                end


            end


            slow_waves(probe_no).DOWN_peaktimes = DOWN_peaks_shank;
            slow_waves(probe_no).DOWN_peaks_zscore = DOWN_peaks_zscore;
            % slow_waves(probe_no).DOWN_peaks_latency = peaks_latency;
            % slow_waves(probe_no).DOWN_traveling = DOWN_traveling;
        end
        disp('cortical wave analysis finished')
        toc
        %


        disp('Sharp wave ripple analysis started')
        tic
        %%%%%%%%%%%% sharp wave direction during ripple peaktimes
        % -1 is posterior -> anterior, 0 is no delay or noisy delay and 1 is anterior -> posterior
        filterparms.deltafilter = [9 17];

        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;


            if length(LFP(1).best_HPC) < length(LFP(2).best_HPC)
                time_idx = 1:length(LFP(1).best_HPC);
            elseif length(LFP(1).best_HPC) > length(LFP(2).best_HPC)

                time_idx = 1:length(LFP(2).best_HPC);
            else
                time_idx = 1:length(LFP(1).best_HPC) ;
            end

            lfp.samplingRate = round(1/mean(diff(LFP(probe_no).tvec)));
            lfp.timestamps = LFP(probe_no).tvec(time_idx);
            tvec = LFP(probe_no).tvec(time_idx);
            lfp.data=[];
            % DOWN_peaks_shank = [];
            sharp_wave_peaks_shank=[];
            peaks_latency = [];
            wave_traveling = [];
            sharp_wave_zscore_shank=[];


            ripples(probe_no).sharp_wave_peaktimes = [];
            ripples(probe_no).sharp_wave_zscore = [];

            ripples(probe_no).SWR_peaktimes = [];
            ripples(probe_no).SWR_zscore = [];
            % ripples(probe_no).SWR_traveling = [];

            ripples(probe_no).probe_hemisphere = [];
            ripples(probe_no).shank_id = [];

            if isfield(LFP(probe_no),'best_HPC') % if exist best HPC channels
                probe_id = [];

                if length(LFP)==1
                    lfp.data= [LFP(probe_no).best_HPC(:,time_idx)'];
                    probe_hemisphere = probe_no*ones(1,length(LFP(probe_no).best_HPC_shank_id));
                    ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id LFP(2).best_HPC_shank_id];
                else
                    lfp.data= [LFP(1).best_HPC(:,time_idx)' LFP(2).best_HPC(:,time_idx)'];
                    probe_hemisphere = [ones(1,length(LFP(1).best_HPC_shank_id)) 2*ones(1,length(LFP(2).best_HPC_shank_id))];
                    if size(LFP(probe_no).best_HPC_shank_id,1)==1
                        ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id LFP(2).best_HPC_shank_id];
                    else
                        ripples(probe_no).shank_id = [LFP(1).best_HPC_shank_id' LFP(2).best_HPC_shank_id'];
                    end
                end

                ripples(probe_no).probe_hemisphere = probe_hemisphere;

                if nprobe == 1 % Only need to grab once
                    % grab delta LFP
                    deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
                    zscored_LFP = zscore(deltaLFP.data);
                    SO_phase_HPC_LFP = deltaLFP.phase;
                    SO_amplitude_HPC_LFP = zscore(deltaLFP.amp);
                    deltaLFP = [];

                    % grab spindles LFP
                    filter_type  = 'bandpass';
                    passband = [20 50];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_spindle = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_spindle,1, lfp.data(:,nShank));
                    end
                    spindle_amplitude_HPC_LFP = zscore(envelope(signal,round(lfp.samplingRate/5),'rms'));
                    spindle_amplitude_HPC_LFP = smoothdata(spindle_amplitude_HPC_LFP,'gaussian',round(lfp.samplingRate/5)); % envelop amplitude of spindles
                    spindle_phase_HPC_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];

                    % grab ripples LFP
                    filter_type  = 'bandpass';
                    passband = [300 600];
                    filter_order = round(6*lfp.samplingRate/(max(passband)-min(passband)));  % creates filter for ripple
                    norm_freq_range = passband/(lfp.samplingRate/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
                    b_ripple = fir1(filter_order, norm_freq_range,filter_type);

                    signal = [];
                    for nShank = 1:length(probe_hemisphere)
                        signal(:,nShank) = filtfilt(b_ripple,1, lfp.data(:,nShank));
                    end
                    zscored_LFP = zscore(abs(hilbert(signal)));
                    ripple_amplitude_LFP = zscored_LFP; % z scored amplitude
                    ripple_phase_LFP = angle(hilbert(signal)); % spindle phase
                    signal = [];
                end

                % [ordered_xcoord,~]=sort(LFP(probe_no).best_V1_xcoord);
                zscored_LFP = [];
                zscored_LFP = ripple_amplitude_LFP;
                for nevent = 1:length(ripples(probe_no).onset)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));


                    for nShank=1:length(probe_hemisphere)

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        if ~isempty(peak_id)
                            [~,temp]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec(tidx(peak_id))));
                            sharp_wave_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            sharp_wave_zscore_shank(nShank,nevent) = zscored_LFP(tidx(peak_id(temp)),nShank);
                        else
                            sharp_wave_peaks_shank(nShank,nevent) = nan;
                            sharp_wave_zscore_shank(nShank,nevent) = nan;
                        end
                    end

                    % % (diff(ordered_xcoord)/1000000)'./diff(DOWN_peaks_shank(:,nevent))
                    % %
                    % % diff(DOWN_peaks_shank(:,425))
                    %
                    % % Putatively
                    % for h = unique(probe_hemisphere)
                    %     sharp_wave_peaks_shank_temp = sharp_wave_peaks_shank(probe_hemisphere==h,nevent);
                    %
                    %     if sum(~isnan(sharp_wave_peaks_shank_temp))==4 % if delta peaks on four shanks
                    %         peaks_latency(h,nevent) = mean(diff(sharp_wave_peaks_shank_temp),'omitnan');
                    %     elseif sum(~isnan(sharp_wave_peaks_shank_temp))==3 % if delta peaks on three shanks
                    %         skipped_shank= diff(LFP(probe_no).best_HPC_shank_id(~isnan(sharp_wave_peaks_shank_temp)))>1;
                    %         if sum(skipped_shank)>0 % if delta peak skipped one shank
                    %
                    %             if skipped_shank(1)==1 % if shanks [1 2 4] then delay using 1 -> 2
                    %                 peaks_latency(h,nevent) = sharp_wave_peaks_shank_temp(1)-sharp_wave_peaks_shank_temp(2);
                    %             elseif skipped_shank(2) == 1 % if shanks [1 3 4] then delay using 3 -> 4
                    %                 peaks_latency(h,nevent) = sharp_wave_peaks_shank_temp(2)-sharp_wave_peaks_shank_temp(3);
                    %             end
                    %         else
                    %             peaks_latency(h,nevent) = mean(diff(sharp_wave_peaks_shank_temp),'omitnan');
                    %         end
                    %     elseif sum(~isnan(sharp_wave_peaks_shank_temp))==2 % if delta peaks on two shanks
                    %         peaks_latency(h,nevent) = diff(sharp_wave_peaks_shank_temp(1:2))/abs(LFP(h).best_HPC_shank_id(2)-LFP(h).best_HPC_shank_id(1));
                    %     else % only one peak. Can't calculate latency
                    %         peaks_latency(h,nevent) =nan;
                    %     end
                    %
                    %     % From Patel et al. (2013) https://pmc.ncbi.nlm.nih.gov/articles/PMC3807028/#sec2
                    %     % ripple propogation speed roughly 0.35 m/s or 0.7
                    %     % ms per 250 micron
                    %     % putatively set the minimum delay threshold to be 2
                    %     % miliseconds.
                    %     if h==1
                    %
                    %         % if direction of mean latency is consistent with the latency more than half
                    %         % of the shank latency
                    %         % (i.e. if four shanks, at least 2 jumps between three shanks should be in the same direction as the mean latency)
                    %         % (if three shanks, then still 2 jumps needed)
                    %         if peaks_latency(h,nevent)>0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)>0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = 1; % anterior to posterior
                    %         elseif peaks_latency(h,nevent)<-0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)<0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = -1; % posterior to anterior
                    %         else
                    %             wave_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %         % DOWN_peaks_shank(:,nevent) - DOWN_peaks_shank(1,nevent)
                    %     elseif h==2
                    %         if peaks_latency(h,nevent)>0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)>0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = -1; % posterior to anterior
                    %         elseif peaks_latency(h,nevent)<-0.001 & sum(sharp_wave_peaks_shank_temp(2:end)-sharp_wave_peaks_shank_temp(1)<0)>=length(LFP(h).best_HPC_shank_id)/2
                    %             wave_traveling(h,nevent) = 1; % anterior to posterior
                    %         else
                    %             wave_traveling(h,nevent) = 0; % noisy or standing wave?
                    %         end
                    %     end
                    % end
                    % % nexttile
                    % % hold on;xline(tvec(idx));
                    % % plot(tvec(tidx),zscored_LFP(tidx,:));
                    % % % hold on;xline(median(tvec(tidx(1):tidx(end))))
                    % % xline(sharp_wave_peaks_shank(:,nevent)')
                end

                ripples(probe_no).SWR_peaktimes = sharp_wave_peaks_shank;
                % ripples(probe_no).SWR_latency = peaks_latency;
                ripples(probe_no).SWR_zscore = sharp_wave_zscore_shank;
                % ripples(probe_no).SWR_traveling = wave_traveling;


                sharp_wave_zscore_shank = [];
                sharp_wave_peaks_shank = [];


                zscored_LFP = [];
                zscored_LFP = SO_amplitude_HPC_LFP;

                for nevent = 1:length(ripples(probe_no).onset)
                    %                 ripples(probe_no).peaktimes(nevent)

                    tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent)-0.05 ripples(probe_no).offset(nevent)]);
                    % tidx = FindInInterval(tvec,[slow_waves(probe_no).ints.UP(nevent,1)-0.3 slow_waves(probe_no).ints.UP(nevent,1)+0.3]);
                    tidx=tidx(1):tidx(end);
                    [~,idx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));


                    for nShank=1:length(probe_hemisphere)

                        [~,peak_id] = findpeaks(zscored_LFP(tidx(1):tidx(end),nShank));

                        if ~isempty(peak_id)
                            [~,temp]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec(tidx(peak_id))));
                            sharp_wave_peaks_shank(nShank,nevent) = tvec(tidx(peak_id(temp)));
                            sharp_wave_zscore_shank(nShank,nevent) = zscored_LFP(tidx(peak_id(temp)),nShank);
                        else
                            sharp_wave_peaks_shank(nShank,nevent) = nan;
                            sharp_wave_zscore_shank(nShank,nevent) = nan;
                        end
                    end
                end

                ripples(probe_no).sharp_wave_peaktimes = sharp_wave_peaks_shank;
                ripples(probe_no).sharp_wave_zscore = sharp_wave_zscore_shank;
            end
        end

        disp('Sharp wave ripple analysis started')
        toc

        zscored_LFP = [];



        tic
        disp('UP/DOWN and ripples and spindles phase coupling and amplitude correlation')
        for nprobe = 1:length(session_info(n).probe)
            probe_no = session_info(n).probe(nprobe).probe_id+1;

            %%%%% Ripples
            ref_shank = find(slow_waves(probe_no).shank_id == slow_waves(probe_no).shank(slow_waves(probe_no).channel == slow_waves(probe_no).best_channel)...
                &slow_waves(probe_no).probe_hemisphere == probe_no);
            cortex_ref_shank = ref_shank;
            spindles(nprobe).best_channel = cortex_ref_shank;


            shank_id = find(ripples(probe_no).probe_hemisphere == probe_no);
            HPC_ref_shank = shank_id(ripples(probe_no).best_channel);
            ref_shank= HPC_ref_shank;

            % ripple amplitude
            step_s = 0.02;
            win_s = step_s*2;
            winSamples = round(win_s * lfp.samplingRate);
            stepSamples = round(step_s * lfp.samplingRate);
            event_tidx = [];
            event_index = [];
            nSteps = floor((size(ripple_amplitude_LFP,1) - winSamples) / stepSamples) + 1;
            tvec_interp1 = tvec(1):step_s:tvec(end);

            meanAmp = [];
            for ntidx = 1:nSteps
                idx = (ntidx-1)*stepSamples + (1:winSamples);
                meanAmp(:,ntidx) = mean(ripple_amplitude_LFP(idx, :));
            end

            % Grab 20ms downsampled idx when events happned
            for nevent = 1:size(ripples(probe_no).peaktimes,1)
                tidx = FindInInterval(tvec_interp1,[ripples(probe_no).onset(nevent) ripples(probe_no).offset(nevent)]);

                event_tidx = [event_tidx tidx(1):tidx(end)];
                event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
            end

            amp_corr = [];
            for nchannel = 1:size(meanAmp,1)
                for mchannel = 1:size(meanAmp,1)
                    amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                end
            end

            ripples(probe_no).amp_corr = amp_corr;


            %%%%%%%%% Phase Locking Value (PLV) during ripples
            %%%%%%%%%% transition and amplitude cross correlation PER EVENT
            plv_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
            pd_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
            xcorr_r_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
            xcorr_lag_ripples = nan(length(ripples(probe_no).shank_id),length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));

            cortex_speed = nan(max(ripples(probe_no).probe_hemisphere), length(ripples(probe_no).peaktimes));
            HPC_speed = nan(max(ripples(probe_no).probe_hemisphere), length(ripples(probe_no).peaktimes));

            for nevent = 1:size(ripples(probe_no).onset,1)
                tidx = FindInInterval(tvec,[ripples(probe_no).onset(nevent) ripples(probe_no).offset(nevent)]);

                for nchannel = 1:length(ripples(probe_no).shank_id)
                    for mchannel = 1:length(ripples(probe_no).shank_id)
                        amp1 = ripple_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                        amp2 = ripple_amplitude_LFP(tidx(1):tidx(end),mchannel);
                        [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                        [~,idx] = max(r);
                        xcorr_lag_ripples(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                        xcorr_r_ripples(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag

                        phi1 = ripple_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                        phi2 = ripple_phase_LFP(tidx(1):tidx(end),mchannel);
                        dphi = phi1 - phi2;
                        % Circular mean of phase differences
                        pd_ripples(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                        % phase locking values
                        plv_ripples(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                    end
                end

                for mprobe = 1:max( ripples(probe_no).probe_hemisphere)
                    HPC_phase = angle(mean(exp(1i * ripple_phase_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                    HPC_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(ripple_phase_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                        mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                    cortex_phase = angle(mean(exp(1i * ripple_phase_cortex_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                    cortex_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(ripple_phase_cortex_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                        mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)
                end

            end

            % right hemisphere
            if max(slow_waves(probe_no).probe_hemisphere)==2
                HPC_speed(2,:) = -HPC_speed(2,:);
                cortex_speed(2,:) = -cortex_speed(2,:);
            end

            ripples(probe_no).plv = plv_ripples;
            ripples(probe_no).pd = pd_ripples;
            ripples(probe_no).xcorr_r = xcorr_r_ripples;
            ripples(probe_no).xcorr_lag = xcorr_lag_ripples;
            ripples(probe_no).cortex_speed = cortex_speed;
            ripples(probe_no).HPC_speed = HPC_speed;



            %%%%% Spindles
            ref_shank = cortex_ref_shank;
            % cortex_ref_shank = ref_shank;
            % spindle amplitude
            step_s = 0.02;
            win_s = step_s*2;
            winSamples = round(win_s * lfp.samplingRate);
            stepSamples = round(step_s * lfp.samplingRate);
            event_tidx = [];
            event_index = [];
            nSteps = floor((size(spindle_amplitude_LFP,1) - winSamples) / stepSamples) + 1;
            tvec_interp1 = tvec(1):step_s:tvec(end);

            meanAmp = [];
            for ntidx = 1:nSteps
                idx = (ntidx-1)*stepSamples + (1:winSamples);
                meanAmp(:,ntidx) = mean(spindle_amplitude_LFP(idx, :));
            end

            if ~isempty(spindles(probe_no).peaktimes)
                % Grab 20ms downsampled idx when events happned
                for nevent = 1:size(spindles(probe_no).peaktimes,1)
                    tidx = FindInInterval(tvec_interp1,[spindles(probe_no).onset(nevent) spindles(probe_no).offset(nevent)]);

                    event_tidx = [event_tidx tidx(1):tidx(end)];
                    event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                end

                amp_corr = [];
                for nchannel = 1:size(meanAmp,1)
                    for mchannel = 1:size(meanAmp,1)
                        amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                    end
                end

                spindles(probe_no).amp_corr = amp_corr;



                %%%%%%%%% Phase Locking Value (PLV) during spindles
                %%%%%%%%%% transition and amplitude cross correlation PER EVENT
                plv_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                pd_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                xcorr_r_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                xcorr_lag_spindles = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));

                cortex_speed = nan(max(slow_waves(probe_no).probe_hemisphere), length(spindles(probe_no).peaktimes));
                HPC_speed = nan(max(slow_waves(probe_no).probe_hemisphere), length(spindles(probe_no).peaktimes));

                for nevent = 1:size(spindles(probe_no).onset,1)
                    tidx = FindInInterval(tvec,[spindles(probe_no).onset(nevent) spindles(probe_no).offset(nevent)]);

                    for nchannel = 1:length(slow_waves(probe_no).shank_id)
                        for mchannel = 1:length(slow_waves(probe_no).shank_id)
                            amp1 = spindle_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                            amp2 = spindle_amplitude_LFP(tidx(1):tidx(end),mchannel);
                            [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                            [~,idx] = max(r);
                            xcorr_lag_spindles(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                            xcorr_r_spindles(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag

                            phi1 = spindle_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                            phi2 = spindle_phase_LFP(tidx(1):tidx(end),mchannel);
                            dphi = phi1 - phi2;
                            % Circular mean of phase differences
                            pd_spindles(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                            % phase locking values
                            plv_spindles(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                        end
                    end

                    for mprobe = 1:max( ripples(probe_no).probe_hemisphere)
                        HPC_phase = angle(mean(exp(1i * spindle_phase_HPC_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                        HPC_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(spindle_phase_HPC_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                            mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                        cortex_phase = angle(mean(exp(1i * spindle_phase_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                        cortex_speed(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(spindle_phase_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                            mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                    end

                end

                % right hemisphere
                if max(slow_waves(probe_no).probe_hemisphere)==2
                    HPC_speed(2,:) = -HPC_speed(2,:);
                    cortex_speed(2,:) = -cortex_speed(2,:);
                end

                spindles(probe_no).plv = plv_spindles;
                spindles(probe_no).pd = pd_spindles;
                spindles(probe_no).xcorr_r = xcorr_r_spindles;
                spindles(probe_no).xcorr_lag = xcorr_lag_spindles;

                spindles(probe_no).cortex_speed = cortex_speed;
                spindles(probe_no).HPC_speed = HPC_speed;
            else
                spindles(probe_no).amp_corr = [];
                spindles(probe_no).plv = [];
                spindles(probe_no).pd = [];
                spindles(probe_no).xcorr_r = [];
                spindles(probe_no).xcorr_lag = [];

                spindles(probe_no).cortex_speed = [];
                spindles(probe_no).HPC_speed = [];
            end

            %%%%% UP/DOWN phase locking and amplitude correlation
            %%%%% Ripple and spindle coupling with UP DOWN

            if ~isempty(slow_waves(probe_no).timestamps)
                if size(slow_waves(probe_no).timestamps,1) > 10

                    UP_ints = slow_waves(probe_no).UP_ints;
                    DOWN_ints = slow_waves(probe_no).DOWN_ints;

                    ref_shank = cortex_ref_shank;
                    % Amplitude at D-U transition and U-D transition
                    step_s = 0.02;
                    win_s = step_s*2;
                    winSamples = round(win_s * lfp.samplingRate);
                    stepSamples = round(step_s * lfp.samplingRate);
                    event_tidx = [];
                    event_index = [];
                    nSteps = floor((size(SO_amplitude_LFP,1) - winSamples) / stepSamples) + 1;
                    tvec_interp1 = tvec(1):step_s:tvec(end);

                    meanAmp = [];
                    for ntidx = 1:nSteps
                        idx = (ntidx-1)*stepSamples + (1:winSamples);
                        meanAmp(:,ntidx) = mean(SO_amplitude_LFP(idx, :));
                    end

                    %%%% DOWN-UP transition
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(UP_ints,1)
                        % Default window: 50 ms before UP onset to 100 ms after
                        t_start = UP_ints(nevent,1) - 0.05;
                        t_end   = UP_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous UP offset
                        if nevent > 1 && t_start < UP_ints(nevent-1,2)
                            t_start = UP_ints(nevent-1,2);% Clip to previous UP offset
                        end

                        % Adjust end if overlapping with next UP onset
                        if nevent < size(UP_ints,1) && t_end > UP_ints(nevent+1,1)
                            t_end = UP_ints(nevent,2); % Clip to current UP offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);

                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end

                    slow_waves(probe_no).amp_corr_DU = amp_corr;

                    %%%% UP-DOWN transition
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(DOWN_ints,1)
                        % Default window: 50 ms before DOWN onset to 100 ms after
                        t_start = DOWN_ints(nevent,1) - 0.05;
                        t_end   = DOWN_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous DOWN offset
                        if nevent > 1 && t_start < DOWN_ints(nevent-1,2)
                            t_start = DOWN_ints(nevent-1,2);% Clip to previous DOWN offset
                        end

                        % Adjust end if overlapping with next DOWN onset
                        if nevent < size(DOWN_ints,1) && t_end > DOWN_ints(nevent+1,1)
                            t_end = DOWN_ints(nevent,2); % Clip to current DOWN offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);
                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end

                    slow_waves(probe_no).amp_corr_UD = amp_corr;


                    % UP and DOWN mean Phase
                    phase_UP = nan(length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    phase_DOWN = nan(length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));

                    % DOWN
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(DOWN_ints,1)
                        % Default window: 50 ms before DOWN onset to 100 ms after
                        t_start = DOWN_ints(nevent,1);
                        t_end   = DOWN_ints(nevent,1) + 0.1;

                        % Adjust end if overlapping with next DOWN onset
                        if nevent < size(DOWN_ints,1) && t_end > DOWN_ints(nevent+1,1)
                            t_end = DOWN_ints(nevent,2); % Clip to current DOWN offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);
                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            phase_DOWN(nchannel,nevent)=angle(mean(exp(1i*SO_phase_LFP(tidx(1):tidx(end),nchannel)))); % phase
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);
                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end

                    slow_waves(probe_no).amp_corr_DOWN = amp_corr;
                    slow_waves(probe_no).mean_phase_DOWN = phase_DOWN;

                    % UP
                    event_tidx = [];
                    event_index = [];
                    % Grab 20ms downsampled idx when events happned
                    for nevent = 1:size(UP_ints,1)
                        % Default window: UP onset to 100 ms after
                        t_start = UP_ints(nevent,1);
                        t_end   = UP_ints(nevent,1) + 0.1;

                        % Adjust end if overlapping with next UP onset
                        if nevent < size(UP_ints,1) && t_end > UP_ints(nevent+1,1)
                            t_end = UP_ints(nevent,2); % Clip to current UP offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);
                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            phase_UP(nchannel,nevent)=angle(mean(exp(1i*SO_phase_LFP(tidx(1):tidx(end),nchannel)))); % phase
                        end

                        tidx = FindInInterval(tvec_interp1, [t_start t_end]);
                        event_tidx = [event_tidx tidx(1):tidx(end)];
                        event_index = [event_index nevent*ones(size(tidx(1):tidx(end)))];
                    end

                    amp_corr = [];
                    for nchannel = 1:size(meanAmp,1)
                        for mchannel = 1:size(meanAmp,1)
                            amp_corr(nchannel,mchannel) = corr(meanAmp(nchannel,event_tidx)',meanAmp(mchannel,event_tidx)');
                        end
                    end
                    slow_waves(probe_no).amp_corr_UP = amp_corr;
                    slow_waves(probe_no).mean_phase_UP = phase_UP;
                    % polarhistogram(phase_DOWN(1,:));hold on;polarhistogram(phase_UP(1,:))

                    %%%%%%%%% Phase Locking Value (PLV) at D-U transition and U-D
                    %%%%%%%%%% transition and amplitude cross correlation PER EVENT
                    plv_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    plv_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));
                    pd_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    pd_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));
                    xcorr_r_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    xcorr_r_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));
                    xcorr_lag_DU = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(UP_ints(:,1)));
                    xcorr_lag_UD = nan(length(slow_waves(probe_no).shank_id),length(slow_waves(probe_no).shank_id), length(DOWN_ints(:,1)));

                    cortex_speed_UD = nan(max(slow_waves(probe_no).probe_hemisphere), length(DOWN_ints(:,1)));
                    HPC_speed_UD = nan(max(slow_waves(probe_no).probe_hemisphere), length(DOWN_ints(:,1)));
                    cortex_speed_DU = nan(max(slow_waves(probe_no).probe_hemisphere), length(UP_ints(:,1)));
                    HPC_speed_DU = nan(max(slow_waves(probe_no).probe_hemisphere), length(UP_ints(:,1)));

                    for nevent = 1:size(UP_ints,1)
                        % Default window: 50 ms before UP onset to 100 ms after
                        t_start = UP_ints(nevent,1) - 0.05;
                        t_end   = UP_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous UP offset
                        if nevent > 1 && t_start < UP_ints(nevent-1,2)
                            t_start = UP_ints(nevent-1,2);% Clip to previous UP offset
                        end

                        % Adjust end if overlapping with next UP onset
                        if nevent < size(UP_ints,1) && t_end > UP_ints(nevent+1,1)
                            t_end = UP_ints(nevent,2); % Clip to current UP offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);

                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            for mchannel = 1:length(slow_waves(probe_no).shank_id)
                                amp1 = SO_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                amp2 = SO_amplitude_LFP(tidx(1):tidx(end),mchannel);
                                [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                                [~,idx] = max(r);
                                xcorr_lag_DU(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                                xcorr_r_DU(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag

                                phi1 = SO_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                phi2 = SO_phase_LFP(tidx(1):tidx(end),mchannel);
                                dphi = phi1 - phi2;
                                % Circular mean of phase differences
                                pd_DU(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                                % phase locking values
                                plv_DU(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                            end
                        end


                        for mprobe = 1:max(slow_waves(probe_no).probe_hemisphere)
                            HPC_phase = angle(mean(exp(1i * SO_phase_HPC_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                            HPC_speed_DU(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_HPC_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                                mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)

                            cortex_phase = angle(mean(exp(1i * SO_phase_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                            cortex_speed_DU(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                                mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)
                        end

                    end

                    % right hemisphere
                    if max(slow_waves(probe_no).probe_hemisphere)==2
                        HPC_speed_DU(2,:) = -HPC_speed_DU(2,:);
                        cortex_speed_DU(2,:) = -cortex_speed_DU(2,:);
                    end

                    for nevent = 1:size(DOWN_ints,1)
                        % Default window: 50 ms before DOWN onset to 100 ms after
                        t_start = DOWN_ints(nevent,1) - 0.05;
                        t_end   = DOWN_ints(nevent,1) + 0.1;

                        % Adjust start if overlapping with previous DOWN offset
                        if nevent > 1 && t_start < DOWN_ints(nevent-1,2)
                            t_start = DOWN_ints(nevent-1,2);% Clip to previous DOWN offset
                        end

                        % Adjust end if overlapping with next DOWN onset
                        if nevent < size(DOWN_ints,1) && t_end > DOWN_ints(nevent+1,1)
                            t_end = DOWN_ints(nevent,2); % Clip to current DOWN offset
                        end

                        % Use final t_start and t_end to extract indices
                        tidx = FindInInterval(tvec, [t_start t_end]);

                        for nchannel = 1:length(slow_waves(probe_no).shank_id)
                            for mchannel = 1:length(slow_waves(probe_no).shank_id)
                                amp1 = SO_amplitude_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                amp2 = SO_amplitude_LFP(tidx(1):tidx(end),mchannel);
                                [r,lag] = xcorr(zscore(amp1),zscore(amp2),'coeff');
                                [~,idx] = max(r);
                                xcorr_lag_UD(nchannel,mchannel,nevent) = lag(idx)/lfp.samplingRate; % where is max lag
                                xcorr_r_UD(nchannel,mchannel,nevent) = r(lag==0); % xcorr at zero lag


                                phi1 = SO_phase_LFP(tidx(1):tidx(end),nchannel); % reference channel
                                phi2 = SO_phase_LFP(tidx(1):tidx(end),mchannel);
                                dphi = phi1 - phi2;
                                % Circular mean of phase differences
                                pd_UD(nchannel,mchannel,nevent) = angle(mean(exp(1i * dphi)));
                                % phase locking values
                                plv_UD(nchannel,mchannel,nevent) = abs(mean(exp(1i * dphi)));
                            end
                        end


                        for mprobe = 1:max(slow_waves(probe_no).probe_hemisphere)
                            HPC_phase = angle(mean(exp(1i * SO_phase_HPC_LFP(tidx(1):tidx(end),(ripples(probe_no).probe_hemisphere == mprobe))), 1));
                            HPC_speed_UD(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_HPC_LFP(tidx(1):tidx(end),HPC_ref_shank))))/...
                                mean(diff(unwrap(HPC_phase))./diff(0.25*ripples(probe_no).shank_id(ripples(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)


                            cortex_phase = angle(mean(exp(1i * SO_phase_LFP(tidx(1):tidx(end),(slow_waves(probe_no).probe_hemisphere == mprobe))), 1));
                            cortex_speed_UD(mprobe,nevent) = 2*pi*mean(lfp.samplingRate / (2*pi) * diff(unwrap(SO_phase_LFP(tidx(1):tidx(end),cortex_ref_shank))))/...
                                mean(diff(unwrap(cortex_phase))./diff(0.25*slow_waves(probe_no).shank_id(slow_waves(probe_no).probe_hemisphere == mprobe)));  % speed =  inst frequency / gradient(radians per mm)
                        end
                    end

                    % right hemisphere
                    if max(slow_waves(probe_no).probe_hemisphere)==2
                        HPC_speed_UD(2,:) = -HPC_speed_UD(2,:);
                        cortex_speed_UD(2,:) = -cortex_speed_UD(2,:);
                    end

                    slow_waves(probe_no).xcorr_lag_UD = xcorr_lag_UD;
                    slow_waves(probe_no).xcorr_lag_DU = xcorr_lag_DU;

                    slow_waves(probe_no).xcorr_r_UD = xcorr_r_UD;
                    slow_waves(probe_no).xcorr_r_DU = xcorr_r_DU;

                    slow_waves(probe_no).plv_UD = plv_UD;
                    slow_waves(probe_no).pd_UD = plv_UD;

                    slow_waves(probe_no).plv_DU = plv_DU;
                    slow_waves(probe_no).pd_DU = plv_DU;

                    slow_waves(probe_no).cortex_speed_UD = cortex_speed_UD;
                    slow_waves(probe_no).HPC_speed_UD = HPC_speed_UD;

                    slow_waves(probe_no).cortex_speed_DU = cortex_speed_DU;
                    slow_waves(probe_no).HPC_speed_DU = HPC_speed_DU;

                    % phase and amplitude
                    SO_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    SO_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));

                    spindle_phase_ripple_peaktime = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    spindle_phase_ripple_onset = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).onset));

                    SO_phase_spindle_onset = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    SO_phase_spindle_peaktime = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

                    % SO_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    % spindle_ripples_coupling = nan(length(slow_waves(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    % SO_spindles_coupling = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));
                    ripple_peak_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).peaktimes));
                    ripple_onset_amplitude = nan(length(ripples(probe_no).shank_id), length(ripples(probe_no).onset));
                    spindle_peak_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).peaktimes));
                    spindle_onset_amplitude = nan(length(slow_waves(probe_no).shank_id), length(spindles(probe_no).onset));

                    for nevent = 1:length(ripples(probe_no).onset)
                        [~,tidx]=min(abs(ripples(probe_no).peaktimes(nevent)-tvec));
                        % Phase
                        SO_phase_ripple_peaktime(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_peaktime(:,nevent) = spindle_phase_LFP(tidx,:);
                        ripple_peak_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);

                        [~,tidx]=min(abs(ripples(probe_no).onset(nevent)-tvec));
                        % Phase
                        SO_phase_ripple_onset(:,nevent) = SO_phase_LFP(tidx,:);
                        spindle_phase_ripple_onset(:,nevent) = spindle_phase_LFP(tidx,:);
                        ripple_onset_amplitude(:,nevent) = ripple_amplitude_LFP(tidx,:);
                    end

                    ripples(probe_no).SO_phase_ripple_peaktime = SO_phase_ripple_peaktime;
                    ripples(probe_no).spindle_phase_ripple_peaktime = spindle_phase_ripple_peaktime;
                    ripples(probe_no).ripple_peak_amplitude = ripple_peak_amplitude;

                    ripples(probe_no).SO_phase_ripple_onset = SO_phase_ripple_onset;
                    ripples(probe_no).spindle_phase_ripple_onset = spindle_phase_ripple_onset;
                    ripples(probe_no).ripple_onset_amplitude = ripple_onset_amplitude;

                    if ~isempty(spindles(probe_no).onset)
                        for nevent = 1:length(spindles(probe_no).onset)
                            [~,tidx]=min(abs(spindles(probe_no).peaktimes(nevent)-tvec));
                            % Phase
                            SO_phase_spindle_peaktime(:,nevent) = SO_phase_LFP(tidx,:);
                            spindle_peak_amplitude(:,nevent) = spindle_amplitude_LFP(tidx,:);

                            [~,tidx]=min(abs(spindles(probe_no).onset(nevent)-tvec));
                            % Phase
                            SO_phase_spindle_onset(:,nevent) = SO_phase_LFP(tidx,:);
                            spindle_onset_amplitude(:,nevent) = spindle_amplitude_LFP(tidx,:);
                        end

                        spindles(probe_no).SO_phase_spindle_peaktime = SO_phase_spindle_peaktime;
                        spindles(probe_no).spindle_peak_amplitude = spindle_peak_amplitude;

                        spindles(probe_no).SO_phase_spindle_onset = SO_phase_spindle_onset;
                        spindles(probe_no).spindle_onset_amplitude = spindle_onset_amplitude;
                    end
                else
                    spindles(probe_no).SO_phase_spindle_peaktime = [];
                    spindles(probe_no).spindle_peak_amplitude = [];

                    spindles(probe_no).SO_phase_spindle_onset = [];
                    spindles(probe_no).spindle_onset_amplitude = [];
                end
            end



        end


        toc
        % histogram(slow_waves(1).DOWN_peaks_latency,100,'Normalization','cdf');
        % hold on; xline(prctile(slow_waves(1).DOWN_peaks_latency,50));
        % histogram(slow_waves(2).DOWN_peaks_latency,100,'Normalization','cdf');
        % hold on; xline(prctile(slow_waves(2).DOWN_peaks_latency,50));
        %
        % histogram(slow_waves(1).DOWN_travling,3,'Normalization','probability');hold on;
        % histogram(slow_waves(2).DOWN_travling,3,'Normalization','probability');

        % histogram(slow_waves(1).DOWN_peaks_latency,100,'Normalization','probability');
        % hold on; xline(prctile(slow_waves(1).DOWN_peaks_latency,50));
        % histogram(slow_waves(2).DOWN_peaks_latency,100,'Normalization','probability');
        % hold on; xline(prctile(slow_waves(2).DOWN_peaks_latency,50));
        if contains(stimulus_name{n},'Masa2tracks')
            % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_ripple_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'ripples');
            % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_spindle_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'spindles');
            % % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_markov_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves_markov');
            % save(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_slow_wave_events%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'slow_waves');
        else
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events_band_control.mat'),'slow_waves');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events_band_control.mat'),'spindles');
            % save(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'),'slow_waves_markov');
            save(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events_band_control.mat'),'ripples');
        end
    end
end


%%
addpath(genpath('C:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis'))
addpath(genpath('C:\Users\masah\Documents\GitHub\VR_NPX_analysis'))


clear all
SUBJECTS={'M24016','M24017','M24018','M24062','M24064','M24065'};
option = 'bilateral';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,option);
% Famililar 
% experiment_info=experiment_info([4 5 6 ]);
experiment_info=experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
Stimulus_type = 'SleepChronic';
% experiment_info=experiment_info([6 9 14 19 21 22 27 35 38 40]);
% 1:length(experiment_info)
% [1 2 3 4 6 7 8 9 10 12 14]
session_count = 0;

slow_waves_all = struct();
ripples_all = struct();
spindles_all = struct();
behavioural_state_merged_all = struct();

for nsession =1:length(experiment_info)
    session_info = experiment_info(nsession).session(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    stimulus_name = experiment_info(nsession).StimulusName(contains(experiment_info(nsession).StimulusName,Stimulus_type));
    SUBJECT_experiment_info = subject_session_stimuli_mapping({session_info(1).probe(1).SUBJECT},option);
    % find right date number based on all experiment dates of the subject
    iDate = find([SUBJECT_experiment_info(:).date] == str2double(session_info(1).probe(1).SESSION));
    if isempty(stimulus_name)
        continue
    end
    load(fullfile(session_info(1).probe(1).ANALYSIS_DATAPATH,'..','best_channels.mat'));


    if length(stimulus_name)>1
        if contains(Stimulus_type,'PRE')
            disp('Same stimuli multiple recordings. Will take _2')
            n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
        else

            session_info = session_info(~contains(stimulus_name,'PRE'));
            stimulus_name = stimulus_name(~contains(stimulus_name,'PRE'));

            if length(stimulus_name)>1
                disp('Same stimuli multiple recordings. Will take _2')
                n = find(contains(stimulus_name,'_2')); % Usually because first stimulus is crashed.
            else
                n =1;
            end
        end
    else
        n = 1;
    end

    session_count = session_count + 1;
    options = session_info(n).probe(1);

    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters*.mat'));
    if isempty(DIR)
        continue
    end

    DIR = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
    DIR1 = dir(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));

    if ~isempty(DIR)
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN.mat'));
        session_clusters_RUN=session_clusters;
        clear session_clusters
    end

    if ~isempty(DIR1)
        load(fullfile(options.ANALYSIS_DATAPATH,'..','session_clusters_RUN1.mat'));
        session_clusters_RUN=session_clusters;
        clear session_clusters
    end

    if contains(stimulus_name{n},'Masa2tracks')
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_PSD%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
% %         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_LFP%s.mat',erase(stimulus_name{n},'Masa2tracks'))),'LFP');
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_task_info%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_behaviour%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
% 
%         load(fullfile(options.ANALYSIS_DATAPATH,sprintf('extracted_clusters_ks4%s.mat',erase(stimulus_name{n},'Masa2tracks'))));
%         clusters=clusters_ks4;
    elseif contains(stimulus_name{n},'Sleep')
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        %             load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events_V1.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_candidate_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'behavioural_state_merged.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_ripple_events_band_control.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_events_band_control.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'extracted_slow_wave_markov_events.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_spindle_events_band_control.mat'));
        % load(fullfile(options.ANALYSIS_DATAPATH,'decoded_ripple_events.mat'));
        %             load(fullfile(options.ANALYSIS_DATAPATH,'reactivation_strength.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'..',sprintf('session_clusters_%s.mat',erase(stimulus_name{n},'Chronic'))),'session_clusters');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    else
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_PSD.mat'));
%         load(fullfile(options.ANALYSIS_DATAPATH,'extracted_LFP.mat'),'LFP');
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_task_info.mat'));
        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_behaviour.mat'));

        load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks4.mat'));
        clusters=clusters_ks4;
    end
    if isfield(slow_waves,'detectorinfo')
        ripples = rmfield(ripples, 'detectorinfo');
        spindles = rmfield(spindles, 'detectorinfo');
        slow_waves = rmfield(slow_waves, 'detectorinfo');
    end
    if isfield(slow_waves,'DOWN_peaks_latency')
        slow_waves = rmfield(slow_waves, 'DOWN_peaks_latency');
        slow_waves = rmfield(slow_waves, 'DOWN_traveling');
        slow_waves = rmfield(slow_waves, 'gammaspikecorr');
        slow_waves = rmfield(slow_waves, 'deltaspikecorr');
        slow_waves = rmfield(slow_waves, 'deltagammacorr');        
    end


    if isfield(ripples,'SWR_traveling')
        ripples = rmfield(ripples, 'SWR_traveling');
        ripples = rmfield(ripples, 'SWR_latency');
    end
    % slow_waves_markov = rmfield(slow_waves_markov, 'NREM_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'spike_count');
    % slow_waves_markov = rmfield(slow_waves_markov, 'alpha_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'gamma_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'beta_t');
    % slow_waves_markov = rmfield(slow_waves_markov, 'xi_t');
    % % 
    % for nprobe = 1:2
    %     slow_waves(nprobe).power = [];
    %     slow_waves(nprobe).frequency = [];
    %     slow_waves(nprobe).PSD_slope = [];
    %     slow_waves(nprobe).timebin_edges = [];
    %     slow_waves(nprobe).DOWN_PSD_slope = [];
    %     slow_waves(nprobe).UP_PSD_slope = [];
    %     slow_waves(nprobe).DOWN_delta_power = [];
    %     slow_waves(nprobe).DOWN_peaks_zscore = [];
    %     slow_waves(nprobe).DOWN_peaktimes = [];
    %     slow_waves(nprobe).DOWN_peaks_latency = [];
    %     slow_waves(nprobe).DOWN_traveling = [];
    %     slow_waves(nprobe).probe_hemisphere = [];
    %     slow_waves(nprobe).shank_id = [];
    %     % slow_waves.power = [];
    % 
    % 
    %     slow_waves(nprobe).power = [];
    %     slow_waves(nprobe).frequency = [];
    %     slow_waves(nprobe).PSD_slope = [];
    %     slow_waves(nprobe).timebin_edges = [];
    %     slow_waves(nprobe).DOWN_PSD_slope = [];
    %     slow_waves(nprobe).UP_PSD_slope = [];
    %     slow_waves(nprobe).DOWN_delta_power = [];
    %     slow_waves(nprobe).DOWN_peaks_zscore = [];
    %     slow_waves(nprobe).DOWN_peaktimes = [];
    %     slow_waves(nprobe).DOWN_peaks_latency = [];
    %     slow_waves(nprobe).DOWN_traveling = [];
    %     slow_waves(nprobe).probe_hemisphere = [];
    %     slow_waves(nprobe).shank_id = [];
    % end
    % 
    % From cell structure back to spike times and spike id
    session_clusters_RUN.spike_id=vertcat(session_clusters_RUN.spike_id{:});
    session_clusters_RUN.spike_times=vertcat(session_clusters_RUN.spike_times{:});
    [session_clusters_RUN.spike_times,index] =sort(session_clusters_RUN.spike_times);
    session_clusters_RUN.spike_id=session_clusters_RUN.spike_id(index);

    session_clusters.spike_id=vertcat(session_clusters.spike_id{:});
    session_clusters.spike_times=vertcat(session_clusters.spike_times{:});
    [session_clusters.spike_times,index] =sort(session_clusters.spike_times);
    session_clusters.spike_id=session_clusters.spike_id(index);

    clusters_combined= session_clusters; % SUA from both probes

    if length(clusters) > 1
        MUA_combined = combine_clusters_from_multiple_probes(clusters(1),clusters(2)); % all clusters for MUA
    else
        MUA_combined = clusters;
    end

    clear selected_clusters
    %         spatial_cell_index = find((session_clusters_RUN.peak_percentile(:,1)>0.95&session_clusters_RUN.odd_even_stability(:,1)>0.95) ...
    %             | (session_clusters_RUN.peak_percentile(:,2)>0.95&session_clusters_RUN.odd_even_stability(:,2)>0.95));
    spatial_cell_index = find(session_clusters_RUN.odd_even_stability(:,1)>0.95 ...
        | session_clusters_RUN.odd_even_stability(:,2)>0.95);

    metric_param =[];
    metric_param.cluster_id = @(x) ismember(x,session_clusters_RUN.cluster_id(spatial_cell_index));
    [selected_clusters,cluster_id] = select_clusters(session_clusters_RUN,metric_param);
    x_window = [0 140];
    x_bin_width = 2;
    place_fields = calculate_spatial_cells(selected_clusters,selected_clusters.tvec{1},...
        selected_clusters.position{1},selected_clusters.speed{1},selected_clusters.track_ID_all{1},selected_clusters.start_time_all{1},selected_clusters.end_time_all{1},x_window,x_bin_width);


    % put all behavioural states in one struct
    field_names = fieldnames(behavioural_state_merged);

    %     for nprobe = 1:length(ripples)
    %         probe_no = session_info(n).probe(nprobe).probe_hemisphere;
    for iField =1:length(field_names)
        behavioural_state_merged_all.(field_names{iField}){session_count} = behavioural_state_merged.(field_names{iField});
    end
    %     end

    % put all ripples in one struct
    field_names = fieldnames(ripples);

    for nprobe = 1:length(ripples)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        ripples_all(probe_no).subject(session_count,:) = options.SUBJECT;
        ripples_all(probe_no).session(session_count,:) = options.SESSION;
        ripples_all(probe_no).session_day(session_count) = iDate;

        for iField =1:length(field_names)
            if session_count == 1
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_ripple_peaktime','spindle_phase_ripple_peaktime','ripple_peak_amplitude',...
                        'SO_phase_ripple_onset','spindle_phase_ripple_onset','ripple_onset_amplitude'})
                    ripples_all(probe_no).(field_names{iField}){session_count} = ripples(nprobe).(field_names{iField});
                else
                    ripples_all(probe_no).(field_names{iField}) = ripples(nprobe).(field_names{iField});
                end


            else
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_ripple_peaktime','spindle_phase_ripple_peaktime','ripple_peak_amplitude',...
                        'SO_phase_ripple_onset','spindle_phase_ripple_onset','ripple_onset_amplitude'})
                    ripples_all(probe_no).(field_names{iField}){session_count} = ripples(nprobe).(field_names{iField});
                else
                    A = ripples_all(probe_no).(field_names{iField});
                    B = ripples(nprobe).(field_names{iField});
                    try
                        ripples_all(probe_no).(field_names{iField}) = [A;B];
                    catch
                        ripples_all(probe_no).(field_names{iField}) = [A B];
                    end
                end
            end
        end

        session_count_events = repmat(session_count,size(ripples(nprobe).onset));
        if session_count == 1
            ripples_all(probe_no).session_count = session_count_events;
        else
            ripples_all(probe_no).session_count = [ripples_all(probe_no).session_count;session_count_events];
        end
    end

    % put all spindles in one struct
    field_names = fieldnames(spindles);

    for nprobe = 1:length(spindles)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        spindles_all(probe_no).subject(session_count,:) = options.SUBJECT;
        spindles_all(probe_no).session(session_count,:) = options.SESSION;
        spindles_all(probe_no).session_day(session_count) = iDate;
        
        for iField =1:length(field_names)
            if session_count == 1
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_spindle_peaktime','SO_phase_spindle_onset','spindle_peak_amplitude','spindle_onset_amplitude'})
                    spindles_all(probe_no).(field_names{iField}){session_count} = spindles(nprobe).(field_names{iField});
                else
                    spindles_all(probe_no).(field_names{iField}) = spindles(nprobe).(field_names{iField});

                end
            else
                if  ismember(field_names{iField},{'sharp_wave_zscore','sharp_wave_peaktimes','SWR_zscore','SWR_peaktimes','shank_id','probe_hemisphere',...
                        'amp_corr','xcorr_lag','xcorr_r','pd','plv','SO_phase_spindle_peaktime','SO_phase_spindle_onset','spindle_peak_amplitude','spindle_onset_amplitude'})
                    spindles_all(probe_no).(field_names{iField}){session_count} = spindles(nprobe).(field_names{iField});
                else
                    A = spindles_all(probe_no).(field_names{iField});
                    B = spindles(nprobe).(field_names{iField});
                    try
                        spindles_all(probe_no).(field_names{iField}) = [A;B];
                    catch
                        spindles_all(probe_no).(field_names{iField}) = [A B];
                    end
                end
            end
        end

        session_count_events = repmat(session_count,size(spindles(nprobe).onset));
        if session_count == 1
            spindles_all(probe_no).session_count = session_count_events;
        else
            spindles_all(probe_no).session_count = [spindles_all(probe_no).session_count;session_count_events];
        end
    end

    % put all slow waves in one struct
    field_names = fieldnames(slow_waves);
%     field_names = fieldnames(slow_waves_markov);
%     slow_waves = slow_waves_markov;
    
    for nprobe = 1:length(slow_waves)
        probe_no = session_info(n).probe(nprobe).probe_hemisphere;
        slow_waves_all(probe_no).subject(session_count,:) = options.SUBJECT;
        slow_waves_all(probe_no).session(session_count,:) = options.SESSION;
        slow_waves_all(probe_no).session_day(session_count) = iDate;

        for iField =1:length(field_names)
            if session_count == 1
                if strcmp(field_names{iField},'ints')
                    slow_waves_all(probe_no).UP_ints = slow_waves(nprobe).ints.UP;
                    slow_waves_all(probe_no).DOWN_ints = slow_waves(nprobe).ints.DOWN;

                elseif  ismember(field_names{iField},{'power','timebin_edges','PSD_slope','frequency',...
                        'deltaspikecorr','gammaspikecorr','deltagammacorr','channel','shank','depth',...
                        'xcoord','gamma_t','viterbi_states','p','DOWN_peaktimes','DOWN_peaks_zscore','shank_id','probe_hemisphere',...
                        'amp_corr_DU','amp_corr_UD','amp_corr_UP','amp_corr_DOWN','mean_phase_UP','mean_phase_DOWN',...
                        'xcorr_lag_UD','xcorr_lag_DU','xcorr_r_UD','xcorr_r_DU','plv_UD','plv_DU','pd_UD','pd_DU'})
                    slow_waves_all(probe_no).(field_names{iField}){session_count} = slow_waves(nprobe).(field_names{iField});
                else
                    slow_waves_all(probe_no).(field_names{iField}) = slow_waves(nprobe).(field_names{iField});
                end
            else
                if strcmp(field_names{iField},'ints')
                    A = slow_waves_all(probe_no).UP_ints;
                    B = slow_waves(nprobe).ints.UP;
                    try
                        slow_waves_all(probe_no).UP_ints = [A;B];
                    catch
                        slow_waves_all(probe_no).UP_ints = [A B];
                    end


                    A = slow_waves_all(probe_no).DOWN_ints;
                    B = slow_waves(nprobe).ints.DOWN;
                    try
                        slow_waves_all(probe_no).DOWN_ints = [A;B];
                    catch
                        slow_waves_all(probe_no).DOWN_ints = [A B];
                    end
                elseif  ismember(field_names{iField},{'power','timebin_edges','PSD_slope','frequency',...
                        'deltaspikecorr','gammaspikecorr','deltagammacorr','channel','shank','depth',...
                        'xcoord','gamma_t','viterbi_states','p','DOWN_peaktimes','DOWN_peaks_zscore','shank_id','probe_hemisphere',...
                        'amp_corr_DU','amp_corr_UD','amp_corr_UP','amp_corr_DOWN','mean_phase_UP','mean_phase_DOWN',...
                        'xcorr_lag_UD','xcorr_lag_DU','xcorr_r_UD','xcorr_r_DU','plv_UD','plv_DU','pd_UD','pd_DU'})

                    slow_waves_all(probe_no).(field_names{iField}){session_count} = slow_waves(nprobe).(field_names{iField});
                else
                    A = slow_waves_all(probe_no).(field_names{iField});
                    B = slow_waves(nprobe).(field_names{iField});


                    try
                        slow_waves_all(probe_no).(field_names{iField}) = [A;B];
                    catch
                        slow_waves_all(probe_no).(field_names{iField}) = [A B];
                    end
                end
            end
        end
        % slow_waves_all(probe_no).UP_PSD_slope
        % session_count_events = repmat(session_count,[size(slow_waves(nprobe).ints.DOWN,1),1]);
        session_count_events = repmat(session_count,[size(slow_waves(nprobe).DOWN_ints,1),1]);
        if session_count == 1
            slow_waves_all(probe_no).DOWN_session_count = session_count_events;
        else
            slow_waves_all(probe_no).DOWN_session_count = [slow_waves_all(probe_no).DOWN_session_count;session_count_events];
        end

        % session_count_events = repmat(session_count,[size(slow_waves(nprobe).ints.UP,1),1]);
        session_count_events = repmat(session_count,[size(slow_waves(nprobe).UP_ints,1),1]);
        if session_count == 1
            slow_waves_all(probe_no).UP_session_count = session_count_events;
        else
            slow_waves_all(probe_no).UP_session_count = [slow_waves_all(probe_no).UP_session_count;session_count_events];
        end

%         session_count_events = repmat(session_count,[size(slow_waves(nprobe).UP_DOWN_index,1),1]);
%         if session_count == 1
%             slow_waves_all(probe_no).UP_DOWN_session_count = session_count_events;
%         else
%             slow_waves_all(probe_no).UP_DOWN_session_count = [slow_waves_all(probe_no).UP_DOWN_session_count;session_count_events];
%         end
% 
%         session_count_events = repmat(session_count,[size(slow_waves(nprobe).DOWN_UP_index,1),1]);
%         if session_count == 1
%             slow_waves_all(probe_no).DOWN_UP_session_count = session_count_events;
%         else
%             slow_waves_all(probe_no).DOWN_UP_session_count = [slow_waves_all(probe_no).DOWN_UP_session_count;session_count_events];
%         end


        if probe_no==1
            metric_param =[];
            % metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            metric_param.region = @(x) contains(x,'V1_L');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).V1_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_L');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).HPC_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'V1_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).V1_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).V1_spike_id{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).HPC_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).HPC_spike_id{session_count} = selected_clusters.spike_times;
        elseif probe_no==2
            metric_param =[];
            % metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            metric_param.region = @(x) contains(x,'V1_R');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).V1_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_R');
            [selected_clusters,cluster_id] = select_clusters(MUA_combined,metric_param);
            slow_waves_all(probe_no).HPC_MUA_spiketimes{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'V1_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).V1_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).V1_spike_id{session_count} = selected_clusters.spike_times;

            metric_param.region = @(x) contains(x,'HPC_L');
            metric_param.cluster_id = @(x) ismember(x,clusters_combined.cluster_id(spatial_cell_index));
            [selected_clusters,cluster_id] = select_clusters(clusters_combined,metric_param);
            slow_waves_all(probe_no).HPC_spike_times{session_count} = selected_clusters.spike_times;
            slow_waves_all(probe_no).HPC_spike_id{session_count} = selected_clusters.spike_times;
        end
    end
end

% 
% for nsession = 1:max(slow_waves_all(1).DOWN_session_count)
%     UP_index=[];
%     DOWN_index=[];
%     % get events index this session
%     for nprobe = 1:length(slow_waves_all)
%         UP_index{nprobe} = find(slow_waves_all(nprobe).UP_session_count == nsession);
%         DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
% 
%         UP_DOWN_transition = slow_waves_all(nprobe).UP_DOWN_index;
%         DOWN_UP_transition = slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe},1);
% 
%         Exclude_index = find(~ismember(DOWN_UP_transition,UP_DOWN_transition));% Only DOWN followed by UP is included for analysis.
% 
%         % remove DOWN events without UP
%         slow_waves_all(nprobe).DOWN_session_count(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_ints(DOWN_index{nprobe}(Exclude_index),:)=[];
% %         slow_waves_all(nprobe).SWpeakmag(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).timestamps(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_PSD_slope(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_peaks_shank{nsession}(:,DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_peaks_latency(DOWN_index{nprobe}(Exclude_index))=[];
%         slow_waves_all(nprobe).DOWN_traveling(DOWN_index{nprobe}(Exclude_index))=[];
% 
%         % DOWN_index{nprobe} = find(slow_waves_all(nprobe).DOWN_session_count == nsession);
%     end
% end

% ripples_all = rmfield(ripples_all, 'detectorinfo');
% spindles_all = rmfield(spindles_all, 'detectorinfo');
if exist('D:\corticohippocampal_replay')>0
    analysis_folder = 'D:\corticohippocampal_replay';
elseif exist('P:\corticohippocampal_replay')>0
    analysis_folder = 'P:\corticohippocampal_replay';
end

if contains(Stimulus_type,'Sleep') & ~contains(Stimulus_type,'PRE')
    save(fullfile(analysis_folder,'slow_waves_all_POST_band_control.mat'),'slow_waves_all','-v7.3')
    save(fullfile(analysis_folder,'ripples_all_POST_band_control.mat'),'ripples_all','-v7.3')
    save(fullfile(analysis_folder,'spindles_all_POST_band_control.mat'),'spindles_all','-v7.3')
    % save(fullfile(analysis_folder,'behavioural_state_merged_all_POST.mat'),'behavioural_state_merged_all','-v7.3')
elseif contains(Stimulus_type,'PRE')
    % save(fullfile(analysis_folder,'slow_waves_all_PRE.mat'),'slow_waves_all','-v7.3')
    % save(fullfile(analysis_folder,'ripples_all_PRE.mat'),'ripples_all','-v7.3')
    % save(fullfile(analysis_folder,'spindles_all_PRE.mat'),'spindles_all','-v7.3')
    % save(fullfile(analysis_folder,'behavioural_state_merged_all_PRE.mat'),'behavioural_state_merged_all','-v7.3')
else
    % save(fullfile(analysis_folder,'slow_waves_all.mat'),'slow_waves_all','-v7.3')
    % save(fullfile(analysis_folder,'ripples_all.mat'),'ripples_all','-v7.3')
    % save(fullfile(analysis_folder,'spindles_all.mat'),'spindles_all','-v7.3')
    % save(fullfile(analysis_folder,'behavioural_state_merged_all.mat'),'behavioural_state_merged_all','-v7.3')
end