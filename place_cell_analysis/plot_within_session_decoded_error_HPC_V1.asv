function decoding_performance = plot_within_session_decoded_error_HPC_V1(estimated_position_lap_CV_V1,decoded_position_V1,decoded_position_HPC,decoded_position_HPC_combined,decoded_position_V1_combined,place_fields,options)

%         decoding_error_confusion_matrix_HPC{nsession} = [];
position_bin = estimated_position_lap_CV_V1(1).track(1).lap(1).track(1).position_bin_centres;

confusion_matrix = [];
for track_id = 1:length(place_fields)
    for nbin =1:length(position_bin)
        true_pos = actual_position{nsession}{track_id};
        confusion_matrix{track_id}(:,nbin) = histcounts(decoded_position_HPC_combined{nsession}{track_id}(true_pos == position_bin(nbin)...
            & VR_speed{nsession}{track_id}>5),length(position_bin_across_tracks));

        confusion_matrix{track_id}(:,nbin) = histcounts(decoded_position_HPC_combined{nsession}{track_id}(true_pos == position_bin(nbin)...
            & VR_speed{nsession}{track_id}>5),length(position_bin_across_tracks));

        confusion_matrix{track_id}(:,nbin) = confusion_matrix{track_id}(:,nbin)/max(confusion_matrix{track_id}(:,nbin));
    end
end
decoding_performance.confusion_matrix.HPC = confusion_matrix;

nfigure = 1;
fig = figure(nfigure)
fig.Position = [300 150 1250 830];
fig.Name = sgtitle(sprintf('%s %s CV decoding confusion matrix and decoding error',options.SUBJECT,options.SESSION));
subplot(4,4,1)
imagesc(flip([confusion_matrix{1}...
    confusion_matrix{2}]))
hold on
xline(14.5,'LineWidth',1)
yline(14.5,'LineWidth',1)
%         imagesc(decoding_confusion_matrix{2})
colorbar
colormap(flip(gray))
clim([0 1])
%         clim([prctile(reshape([confusion_matrix{1} confusion_matrix{2}],1,[]),10) prctile(reshape([confusion_matrix{1} confusion_matrix{2}],1,[]),90)])
xticks(1:2:length(confusion_matrix{1}))
xticklabels([position_bin(1:2:end) position_bin(1:2:end)])
yticks(1:2:length(confusion_matrix{1}))
yticklabels([position_bin(1:2:end) position_bin(1:2:end)])
xlabel('True Position (cm)')
ylabel('Decoded Position (cm)')
set(gca,"TickDir","out",'box', 'off','Color','none')
title('HPC combined decoding confusion matrix')


if length(session_info(n).probe)>1
    confusion_matrix = [];
    for track_id = 1:length(place_fields)
        for nbin =1:length(position_bin)
            true_pos = actual_position{nsession}{track_id};
            confusion_matrix{track_id}(:,nbin) = histcounts(decoded_position_V1_combined{nsession}{track_id}(true_pos == position_bin(nbin)...
                & VR_speed{nsession}{track_id}>5),length(position_bin_across_tracks));
            confusion_matrix{track_id}(:,nbin) = confusion_matrix{track_id}(:,nbin)/max(confusion_matrix{track_id}(:,nbin));
        end
    end
    decoding_performance.confusion_matrix.V1_combined = confusion_matrix;

    subplot(4,4,2)
    imagesc(flip([confusion_matrix{1}...
        confusion_matrix{2}]))
    hold on
    xline(14.5,'LineWidth',1)
    yline(14.5,'LineWidth',1)
    %         imagesc(decoding_confusion_matrix{2})
    colorbar
    colormap(flip(gray))
    xticks(1:2:length(confusion_matrix{1}))
    xticklabels([position_bin(1:2:end) position_bin(1:2:end)])
    yticks(1:2:length(confusion_matrix{1}))
    yticklabels([position_bin(1:2:end) position_bin(1:2:end)])
    xlabel('True Position (cm)')
    ylabel('Decoded Position (cm)')
    set(gca,"TickDir","out",'box', 'off','Color','none')
    title('V1 combined decoding confusion matrix')
end

confusion_matrix = [];

for nprobe = 1:length(session_info(n).probe)
    probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

    for track_id = 1:length(place_fields)
        for nbin =1:length(position_bin)
            true_pos = actual_position{nsession}{track_id};
            confusion_matrix{probe_hemisphere}{track_id}(:,nbin) = histcounts(decoded_position_V1{probe_hemisphere}{nsession}{track_id}(true_pos == position_bin(nbin)...
                & VR_speed{nsession}{track_id}>5),length(position_bin_across_tracks));
            confusion_matrix{probe_hemisphere}{track_id}(:,nbin) = confusion_matrix{probe_hemisphere}{track_id}(:,nbin)/max(confusion_matrix{probe_hemisphere}{track_id}(:,nbin));
        end
    end
end

decoding_performance.confusion_matrix.V1 = confusion_matrix;


if ~isempty(confusion_matrix{1})
    subplot(4,4,3)
    imagesc(flip([confusion_matrix{1}{1}...
        confusion_matrix{1}{2}]))
    hold on
    xline(14.5,'LineWidth',1)
    yline(14.5,'LineWidth',1)
    %         imagesc(decoding_confusion_matrix{2})
    colorbar
    colormap(flip(gray))
    xticks(1:2:length(confusion_matrix{1}{1}))
    xticklabels([position_bin(1:2:end) position_bin(1:2:end)])
    yticks(1:2:length(confusion_matrix{1}{1}))
    yticklabels([position_bin(1:2:end) position_bin(1:2:end)])
    xlabel('True Position (cm)')
    ylabel('Decoded Position (cm)')
    set(gca,"TickDir","out",'box', 'off','Color','none')
    title('V1 Left')
end

if ~isempty(confusion_matrix{2})
    subplot(4,4,4)
    imagesc(flip([confusion_matrix{2}{1}...
        confusion_matrix{2}{2}]))
    hold on
    xline(14.5,'LineWidth',1)
    yline(14.5,'LineWidth',1)
    %         imagesc(decoding_confusion_matrix{2})
    colorbar
    colormap(flip(gray))
    xticks(1:2:length(confusion_matrix{2}{2}))
    xticklabels([position_bin(1:2:end) position_bin(1:2:end)])
    yticks(1:2:length(confusion_matrix{2}{2}))
    yticklabels([position_bin(1:2:end) position_bin(1:2:end)])
    title('V1 Right')
end

median_lap_decoding_error = [];

for track_id = 1:length(place_fields)
    for temp_track = 1:length(place_fields)
        for lap_id = 1:length(estimated_position_lap_CV_V1(nprobe).track(track_id).lap)

            median_lap_decoding_error{track_id}{temp_track}(lap_id) = nanmedian(decoded_error_HPC_combined{nsession}{track_id}{temp_track}...
                (decoded_position_lap_id{nsession}{track_id} == lap_id & VR_speed{nsession}{track_id}>5));
        end
        %                 scatter(median_lap_decoding_error{track_id}{temp_track})
    end
end


subplot(4,4,5)
data = [median_lap_decoding_error{1}{1} median_lap_decoding_error{1}{2} median_lap_decoding_error{2}{1} median_lap_decoding_error{2}{2}];
label = [10*ones(1,length(median_lap_decoding_error{1}{1})) 20*ones(1,length(median_lap_decoding_error{1}{2}))...
    30*ones(1,length(median_lap_decoding_error{2}{1})) 40*ones(1,length(median_lap_decoding_error{2}{2}))];
beeswarm(label',data','sort_style','rand','overlay','sd'); hold on
decoding_performance.decoding_error.HPC = median_lap_decoding_error;
xlim([0 50])
title('HPC decoding error')
xticks([10 20 30 40])
xticklabels(["Track 1 by T1 template","Track 1 by T2 template","Track 2 by T1 template","Track 2 by T2 template"])

if length(session_info(n).probe)>1
    median_lap_decoding_error = [];

    for track_id = 1:length(place_fields)
        for temp_track = 1:length(place_fields)
            for lap_id = 1:length(estimated_position_lap_CV_V1(nprobe).track(track_id).lap)

                median_lap_decoding_error{track_id}{temp_track}(lap_id) = nanmedian(decoded_error_V1_combined{nsession}{track_id}{temp_track}...
                    (decoded_position_lap_id{nsession}{track_id} == lap_id & VR_speed{nsession}{track_id}>5));
            end

            %                 scatter(median_lap_decoding_error{track_id}{temp_track})
        end
    end

    subplot(4,4,6)
    data = [median_lap_decoding_error{1}{1} median_lap_decoding_error{1}{2} median_lap_decoding_error{2}{1} median_lap_decoding_error{2}{2}];
    label = [10*ones(1,length(median_lap_decoding_error{1}{1})) 20*ones(1,length(median_lap_decoding_error{1}{2}))...
        30*ones(1,length(median_lap_decoding_error{2}{1})) 40*ones(1,length(median_lap_decoding_error{2}{2}))];
    beeswarm(label',data','sort_style','rand','overlay','sd'); hold on
    decoding_performance.decoding_error.V1_combined = median_lap_decoding_error;
    xlim([0 50])
    title('V1 combined decoding error')
    xticks([10 20 30 40])
    xticklabels(["Track 1 by T1 template","Track 1 by T2 template","Track 2 by T1 template","Track 2 by T2 template"])
end

median_lap_decoding_error = [];
for nprobe = 1:length(session_info(n).probe)
    probe_hemisphere = session_info(n).probe(nprobe).probe_hemisphere;

    for track_id = 1:length(place_fields)
        for temp_track = 1:length(place_fields)
            for lap_id = 1:length(estimated_position_lap_CV_V1(nprobe).track(track_id).lap)

                median_lap_decoding_error{probe_hemisphere}{track_id}{temp_track}(lap_id) = nanmedian(decoded_error_V1{probe_hemisphere}{nsession}{track_id}{temp_track}...
                    (decoded_position_lap_id{nsession}{track_id} == lap_id & VR_speed{nsession}{track_id}>5));
            end

            %                 scatter(median_lap_decoding_error{track_id}{temp_track})
        end
    end
end
decoding_performance.decoding_error.V1 = median_lap_decoding_error;

if ~isempty(median_lap_decoding_error{1})

    subplot(4,4,7)
    data = [median_lap_decoding_error{1}{1}{1} median_lap_decoding_error{1}{1}{2} median_lap_decoding_error{1}{2}{1} median_lap_decoding_error{1}{2}{2}];
    label = [10*ones(1,length(median_lap_decoding_error{1}{1}{1})) 20*ones(1,length(median_lap_decoding_error{1}{1}{2}))...
        30*ones(1,length(median_lap_decoding_error{1}{2}{1})) 40*ones(1,length(median_lap_decoding_error{1}{2}{2}))];
    beeswarm(label',data','sort_style','rand','overlay','sd'); hold on
    xlim([0 50])
    title('V1 left decoding error')
    xticks([10 20 30 40])
    xticklabels(["Track 1 by T1 template","Track 1 by T2 template","Track 2 by T1 template","Track 2 by T2 template"])

end

if ~isempty(median_lap_decoding_error{2})

    subplot(4,4,8)
    data = [median_lap_decoding_error{2}{1}{1} median_lap_decoding_error{2}{1}{2} median_lap_decoding_error{2}{2}{1} median_lap_decoding_error{2}{2}{2}];
    label = [10*ones(1,length(median_lap_decoding_error{1}{1}{1})) 20*ones(1,length(median_lap_decoding_error{2}{1}{2}))...
        30*ones(1,length(median_lap_decoding_error{1}{2}{1})) 40*ones(1,length(median_lap_decoding_error{2}{2}{2}))];
    beeswarm(label',data','sort_style','rand','overlay','sd'); hold on
    xlim([0 50])
    title('V1 Right decoding error')
    xticks([10 20 30 40])
    xticklabels(["Track 1 by T1 template","Track 1 by T2 template","Track 2 by T1 template","Track 2 by T2 template"])

end


subplot(4,4,9)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{1}{1}(VR_speed{nsession}{1}>5),decoded_error_V1{1}{nsession}{1}{1}(VR_speed{nsession}{1}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
clim([0 0.5])
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Left V1 left-sided T1')

subplot(4,4,10)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{2}{2}(VR_speed{nsession}{2}>5),decoded_error_V1{2}{nsession}{2}{2}(VR_speed{nsession}{2}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Right V1 right-sided T2')


subplot(4,4,11)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{1}{1}(VR_speed{nsession}{1}>5),decoded_error_V1{2}{nsession}{1}{1}(VR_speed{nsession}{1}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Right V1 left-sided T1')

subplot(4,4,12)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{2}{2}(VR_speed{nsession}{2}>5),decoded_error_V1{1}{nsession}{2}{2}(VR_speed{nsession}{2}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Left V1 right-sided T2')



subplot(4,4,9)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{1}{1}(VR_speed{nsession}{1}>5),decoded_error_V1{1}{nsession}{1}{1}(VR_speed{nsession}{1}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
clim([0 0.5])
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Left V1 left-sided T1')

subplot(4,4,10)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{2}{2}(VR_speed{nsession}{2}>5),decoded_error_V1{2}{nsession}{2}{2}(VR_speed{nsession}{2}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Right V1 right-sided T2')


subplot(4,4,11)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{1}{1}(VR_speed{nsession}{1}>5),decoded_error_V1{2}{nsession}{1}{1}(VR_speed{nsession}{1}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Right V1 left-sided T1')

subplot(4,4,12)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{2}{2}(VR_speed{nsession}{2}>5),decoded_error_V1{1}{nsession}{2}{2}(VR_speed{nsession}{2}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Left V1 right-sided T2')


% Left HPC1 and Left V1
subplot(4,4,9)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{1}{1}(VR_speed{nsession}{1}>5),decoded_error_V1{1}{nsession}{1}{1}(VR_speed{nsession}{1}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
clim([0 0.5])
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Left V1 left-sided T1')

subplot(4,4,10)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{2}{2}(VR_speed{nsession}{2}>5),decoded_error_V1{2}{nsession}{2}{2}(VR_speed{nsession}{2}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Right V1 right-sided T2')


subplot(4,4,11)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{1}{1}(VR_speed{nsession}{1}>5),decoded_error_V1{2}{nsession}{1}{1}(VR_speed{nsession}{1}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Right V1 left-sided T1')

subplot(4,4,12)
[N,Xedges,Yedges,binX,binY] = histcounts2(decoded_error_HPC_combined{nsession}{2}{2}(VR_speed{nsession}{2}>5),decoded_error_V1{1}{nsession}{2}{2}(VR_speed{nsession}{2}>5),-140:10:140,-140:10:140);
imagesc((flip(N'))/max(max(N)))
clim([0 0.5])
xticks(1:2:length(Xedges))
xticklabels(Xedges(1:2:end))
yticks(1:2:length(Yedges))
yticklabels(Yedges(1:2:end))
colorbar
xlabel('HPC decoded error (cm)')
ylabel('V1 decoded error (cm)')
title('Left V1 right-sided T2')

sgtitle(sprintf('%s %s CV decoding confusion matrix and decoding error',options.SUBJECT,options.SESSION))