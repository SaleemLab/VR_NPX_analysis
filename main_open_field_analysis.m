SUBJECTS = {'M23087'};
options = 'bilateral';
Stimulus_type = 'OpenFieldChronic';
experiment_info = subject_session_stimuli_mapping(SUBJECTS,options);
extract_and_preprocess_NPX_batch(experiment_info,Stimulus_type)
%    extract_and_preprocess_NPX(session_info,Stimulus_type)

SUBJECT = 'M23087';
SESSION = '20231214';
load(fullfile(ROOTPATH,'DATA','SUBJECTS',SUBJECT,'analysis',SESSION,Stimulus_type,'session_info.mat'));


for nprobe = 1:length(session_info.probe)
    load(fullfile(session_info.probe(nprobe).ANALYSIS_DATAPATH,"extracted_behaviour.mat"));
    load(fullfile(session_info.probe(nprobe).ANALYSIS_DATAPATH,"extracted_clusters.mat"));
    HPC_boundary = [max(clusters.peak_depth)-800 max(clusters.peak_depth)-800-1500];% later replaced by best channel determined based on LFP profile
    cluster_id = clusters.unit_id(find(clusters.peak_depth<= HPC_boundary(1)& clusters.peak_depth>= HPC_boundary(2) & clusters.cell_type == 1));% Find pyramidal neurons

    time_bin_width = mean(diff(Behaviour.sglxTime));

    fig = figure
    fig.Position = [300 50 1400 920]

    count =1
    for unit = 1:length(cluster_id)
        if count == 13

            fig = figure
            fig.Position = [300 50 1400 920]

            count = 1;
        end

        this_cluster_index = find(clusters.unit_id == cluster_id(unit));
        clusters.spike_times(clusters.spike_id == cluster_id(unit));

%        spped = sqrt([0 diff(Behaviour.X)].^2+[0 diff(Behaviour.Y)].^2)
        immobility = kmeans(Behaviour.mobility',2)-2;% find moving and not moving
        immobility(immobility==-1) = 1; 

        X_during_spike = interp1(Behaviour.sglxTime(immobility == 0),Behaviour.X(immobility == 0),clusters.spike_times(clusters.spike_id == cluster_id(unit)),'nearest'); %interpolates position into spike time
        Y_during_spike = interp1(Behaviour.sglxTime(immobility == 0),Behaviour.Y(immobility == 0),clusters.spike_times(clusters.spike_id == cluster_id(unit)),'nearest'); %interpolates position into spike time

        [spike_hist,Xedges,Yedges] = histcounts2(X_during_spike,Y_during_spike,linspace(0,1,21),linspace(0,1,21));%linspace(min(Behaviour.X),max(Behaviour.X),10),linspace(min(Behaviour.Y),max(Behaviour.Y),10)
        [XY_hist] = time_bin_width.*histcounts2(Behaviour.X(immobility == 0),Behaviour.Y(immobility == 0),linspace(0,1,21),linspace(0,1,21));


        spike_hist = spike_hist./XY_hist;
        spike_hist = spike_hist';

        subplot(3,4,count)
        scatter(Behaviour.X,Behaviour.Y,'filled','MarkerFaceAlpha','0.01'); hold on
        scatter(X_during_spike,Y_during_spike,'filled','r')
        count = count + 1;
        title(cluster_id(unit))

        subplot(3,4,count)
        %         spike_hist(isnan(spike_hist))=0;
        %         scal_f = 10; % scale image by this before...
        %         sigma = 3; % ...filtering by this
        %         spike_hist_smoothed = imresize(spike_hist,scal_f);
%         spike_hist_smoothed = imgaussfilt(spike_hist_smoothed,sigma);
%         imagesc(flip(spike_hist_smoothed));

                imagesc(flip(spike_hist))
        %         set(gca,'YDir','normal')
        %         set(gca,'XDir','normal')
        colormap(flip(gray))
        title(cluster_id(unit))

        count = count + 1;
    end

end