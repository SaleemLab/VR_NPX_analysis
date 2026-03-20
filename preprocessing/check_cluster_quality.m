function  check_cluster_quality(session_info,Stimulus_type)

options = session_info.probe(1);

if contains(Stimulus_type,'Masa2tracks')
    % If Masa2tracks, PRE, RUN and/or POST saved in one folder
    load(fullfile(options.ANALYSIS_DATAPATH,...
        sprintf('extracted_clusters_ks2%s.mat',erase(options.StimulusName,Stimulus_type))))
    load(fullfile(options.ANALYSIS_DATAPATH,...
        sprintf('extracted_clusters_ks3%s.mat',erase(options.StimulusName,Stimulus_type))))
    load(fullfile(options.ANALYSIS_DATAPATH,...
        sprintf('extracted_clusters%s.mat',erase(options.StimulusName,Stimulus_type))))
    %     save(fullfile(options.ANALYSIS_DATAPATH,...
    %         sprintf('extracted_spikes%s.mat',erase(stimulus_name,Stimulus_type))),'spikes')
else
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks2.mat'))
    load(fullfile(options.ANALYSIS_DATAPATH,'extracted_clusters_ks3.mat'))
%     save(fullfile(options.ANALYSIS_DATAPATH,'extracted_spikes.mat'),'spikes')
end


% isi violation
figure
count = 1
for nprobe = 1:length(clusters_ks2)
    subplot(3,2,count)
    h(1) = histogram(clusters_ks2(nprobe).isi_violations_ratio,0:0.05:1.2,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none'); hold on
    h(2) = histogram(clusters_ks3(nprobe).isi_violations_ratio,0:0.05:1.2,'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none')
    h(3) = histogram(clusters(nprobe).isi_viol,0:0.05:1.2,'FaceColor','m','FaceAlpha',0.5,'EdgeColor','none')
    xlabel('ISI violation ratio')
    count = count + 1;
    legend([h(1) h(2) h(3)],{'KS2','KS3','old KS3'})
    title(sprintf('Probe %i ISI violation ratio',nprobe))
end

% amplitude cutoff
for nprobe = 1:length(clusters_ks2)
    subplot(3,2,count)
    h(1) = histogram(clusters_ks2(nprobe).amplitude_cutoff,0:0.01:0.5,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none'); hold on
    h(2) = histogram(clusters_ks3(nprobe).amplitude_cutoff,0:0.01:0.5,'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none')
    h(3) = histogram(clusters(nprobe).amplitude_cutoff,0:0.01:0.5,'FaceColor','m','FaceAlpha',0.5,'EdgeColor','none')
    count = count + 1;
    xlabel('amplitude cutoff')
    legend([h(1) h(2) h(3)],{'KS2','KS3','old KS3'})
    title(sprintf('Probe %i amplitude cutoff',nprobe))
end

% amplitude
for nprobe = 1:length(clusters_ks2)
    subplot(3,2,count)
    h(1) = histogram(clusters_ks2(nprobe).amplitude_median,0:5:400,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none'); hold on
    h(2) = histogram(clusters_ks3(nprobe).amplitude_median,0:5:400,'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none')
    h(3) = histogram(clusters(nprobe).amplitude,0:5:400,'FaceColor','m','FaceAlpha',0.5,'EdgeColor','none')
    count = count + 1;
    xlabel('amplitude')
    legend([h(1) h(2) h(3)],{'KS2','KS3','old KS3'})
    title(sprintf('Probe %i amplitude cutoff',nprobe))
end




for nprobe = 1:length(clusters_ks2)
    fig = figure;
    fig.Position = [680 185 1000 800];
    subplot(2,2,1)
    good_unit_index = (clusters_ks2(nprobe).amplitude_cutoff <= 0.1...
        &clusters_ks2(nprobe).isi_violations_ratio <= 0.1...
        &clusters_ks2(nprobe).amplitude_median >=50);

    h(1) = scatter3(clusters_ks2(nprobe).isi_violations_ratio(good_unit_index),clusters_ks2(nprobe).amplitude_cutoff(good_unit_index),clusters_ks2(nprobe).amplitude_median(good_unit_index),'filled','MarkerFaceAlpha',0.1)
    xlabel('ISI violation ratio')
    xlim([0 0.2])
    ylabel('amplitude cutoff')
    ylim([0 0.2])
    zlabel('amplitude')
    title('KS2')


    subplot(2,2,2)
    good_unit_index = (clusters_ks3(nprobe).amplitude_cutoff <= 0.1...
        &clusters_ks3(nprobe).isi_violations_ratio <= 0.1...
        &clusters_ks3(nprobe).amplitude_median >=50);

    h(1) = scatter3(clusters_ks3(nprobe).isi_violations_ratio(good_unit_index),clusters_ks3(nprobe).amplitude_cutoff(good_unit_index),clusters_ks3(nprobe).amplitude_median(good_unit_index),'filled','MarkerFaceAlpha',0.1)
    xlabel('ISI violation ratio')
    xlim([0 0.2])
    ylabel('amplitude cutoff')
    ylim([0 0.2])
    zlabel('amplitude')
    title('KS3')

    subplot(2,2,3)
    good_unit_index = (clusters(nprobe).amplitude_cutoff <= 0.1...
        &clusters(nprobe).isi_viol <= 0.1...
        &clusters(nprobe).amplitude >=50);

    h(1) = scatter3(clusters(nprobe).isi_viol(good_unit_index),clusters(nprobe).amplitude_cutoff(good_unit_index),clusters(nprobe).amplitude(good_unit_index),'filled','MarkerFaceAlpha',0.1)
    xlabel('ISI violation ratio')
    xlim([0 0.2])
    ylabel('amplitude cutoff')
    ylim([0 0.2])
    zlabel('amplitude')
    title('Old KS3')

    sgtitle(sprintf('probe %i',nprobe))
end