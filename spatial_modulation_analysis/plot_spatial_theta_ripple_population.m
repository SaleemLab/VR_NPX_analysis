function plot_spatial_theta_ripple_population(cell_selection,selected_cells_L,selected_cells_R,place_fields_all_L,place_fields_all_R)

if isempty(cell_selection)
    cell_selection = 'Nan';
end

fig = figure
fig.Position = [931 100 850 860]
fig.Name = sprintf('All sessions spatial map (%s)',cell_selection)
sgtitle(sprintf('All sessions spatial map (%s)',cell_selection))
average_map = normalize([place_fields_all_L(1).average_map(:,selected_cells_L);place_fields_all_L(2).average_map(:,selected_cells_L)],'range')';
average_map = reshape(average_map,size(average_map,1),size(place_fields_all_L(1).average_map,1),[]);
subplot(2,2,1)
imagesc(0:2:140,[],squeeze(average_map(:,:,1)));
colormap(flip(gray))
set(gca,'YDir','normal')
xlabel('position')
title('Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,2)
imagesc(0:2:140,[],squeeze(average_map(:,:,2)));
set(gca,'YDir','normal')
colormap(flip(gray))
xlabel('position')
title('Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

average_map = normalize([place_fields_all_R(1).average_map(:,selected_cells_R);place_fields_all_R(2).average_map(:,selected_cells_R)],'range')';
average_map = reshape(average_map,size(average_map,1),size(place_fields_all_R(1).average_map,1),[]);
subplot(2,2,3)
imagesc(0:2:140,[],squeeze(average_map(:,:,1)));
colormap(flip(gray))
set(gca,'YDir','normal')
xlabel('position')
title('Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
imagesc(0:2:140,[],squeeze(average_map(:,:,2)));
set(gca,'YDir','normal')
colormap(flip(gray))
xlabel('position')
title('Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)





fig = figure
fig.Position = [500 100 400 860]
fig.Name = sprintf('All sessions within block SMI (%s)',cell_selection)
sgtitle(sprintf('All sessions  within block SMI (%s)',cell_selection))
within_block_SMI = [place_fields_all_L(1).within_block_SMI(:,selected_cells_L);place_fields_all_L(2).within_block_SMI(:,selected_cells_L)]';
within_block_SMI = reshape(within_block_SMI,size(within_block_SMI,1),size(place_fields_all_L(1).within_block_SMI,1),[]);
subplot(2,2,1)
imagesc(squeeze(within_block_SMI(:,:,1)));
colorbar
set(gca,'YDir','normal')
clim([0 1])
xlabel('lap')
title('Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,2)
imagesc(squeeze(within_block_SMI(:,:,2)));
colorbar
set(gca,'YDir','normal')
clim([0 1])
xlabel('lap')
title('Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

within_block_SMI = [place_fields_all_R(1).within_block_SMI(:,selected_cells_R);place_fields_all_R(2).within_block_SMI(:,selected_cells_R)]';
within_block_SMI = reshape(within_block_SMI,size(within_block_SMI,1),size(place_fields_all_R(1).within_block_SMI,1),[]);
subplot(2,2,3)
imagesc(squeeze(within_block_SMI(:,:,1)));
colorbar
set(gca,'YDir','normal')
clim([0 1])
xlabel('lap')
title('Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
imagesc(squeeze(within_block_SMI(:,:,2)));
colorbar
set(gca,'YDir','normal')
clim([0 1])
xlabel('lap')
title('Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)



fig = figure
fig.Position = [931 100 850 860]
fig.Name = sprintf('All sessions ripple PSTH (%s)',cell_selection)
sgtitle(sprintf('All sessions ripple PSTH (%s)',cell_selection))
% rippple_PSTH = [place_fields_all_L(1).ripple_PSTH(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH(:,selected_cells_L)]';
if length(selected_cells_L)==length(selected_cells_R)
    rippple_PSTH = normalize([place_fields_all_L(1).ripple_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH_zscored(:,selected_cells_L);...
        place_fields_all_R(1).ripple_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH_zscored(:,selected_cells_R)],'range')';
%     rippple_PSTH = [place_fields_all_L(1).ripple_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH_zscored(:,selected_cells_L);...
%         place_fields_all_R(1).ripple_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH_zscored(:,selected_cells_R)]';
else
    rippple_PSTH = normalize([place_fields_all_L(1).ripple_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH_zscored(:,selected_cells_L)],'range')';
%     rippple_PSTH = [place_fields_all_L(1).ripple_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH_zscored(:,selected_cells_L)]';
    rippple_PSTH = [rippple_PSTH rippple_PSTH];
end
rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_L(1).ripple_PSTH,1),[]);
subplot(2,2,1)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,1)));
hold on;xline(0,'r')
set(gca,'YDir','normal')
xlabel('Time')
title('L ripple on Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% colorbar
% clim([0 5])

subplot(2,2,2)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,2)));
hold on;xline(0,'r')'
set(gca,'YDir','normal')
xlabel('Time')
title('L ripple on Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% colorbar
% clim([0 5])

% rippple_PSTH = [place_fields_all_R(1).ripple_PSTH(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH(:,selected_cells_R)]';
% rippple_PSTH = normalize([place_fields_all_R(1).ripple_PSTH(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH(:,selected_cells_R)],'range')';
% rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_R(1).ripple_PSTH,1),[]);
if length(selected_cells_L)==length(selected_cells_R)
    % rippple_PSTH = [place_fields_all_L(1).ripple_PSTH(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH(:,selected_cells_L)]';
    rippple_PSTH = normalize([place_fields_all_L(1).ripple_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH_zscored(:,selected_cells_L);...
        place_fields_all_R(1).ripple_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH_zscored(:,selected_cells_R)],'range')';
%     rippple_PSTH = [place_fields_all_L(1).ripple_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH_zscored(:,selected_cells_L);...
%         place_fields_all_R(1).ripple_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH_zscored(:,selected_cells_R)]';
else
%     rippple_PSTH = normalize([place_fields_all_R(1).ripple_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH_zscored(:,selected_cells_R)],'range')';
    rippple_PSTH = [place_fields_all_R(1).ripple_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH_zscored(:,selected_cells_R)]';
    rippple_PSTH = [rippple_PSTH rippple_PSTH];
end
rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_L(1).ripple_PSTH,1),[]);

subplot(2,2,3)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,3)));
hold on;xline(0,'r')
set(gca,'YDir','normal')
xlabel('Time')
title('R ripple on Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% colorbar
% clim([0 5])

subplot(2,2,4)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,4)));
hold on;xline(0,'r')
set(gca,'YDir','normal')
xlabel('Time')
title('R ripple on Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
% colorbar
% clim([0 5])

fig = figure
fig.Position = [931 100 850 860]
fig.Name = sprintf('All sessions  V1 bursting events PSTH (%s)',cell_selection)
sgtitle(sprintf('All sessions V1 bursting events PSTH (%s)',cell_selection))
if length(selected_cells_L)==length(selected_cells_R)
    % rippple_PSTH = [place_fields_all_L(1).ripple_PSTH(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH(:,selected_cells_L)]';
    rippple_PSTH = normalize([place_fields_all_L(1).V1_event_PSTH(:,selected_cells_L);place_fields_all_L(2).V1_event_PSTH(:,selected_cells_L);...
        place_fields_all_R(1).V1_event_PSTH(:,selected_cells_R);place_fields_all_R(2).V1_event_PSTH(:,selected_cells_R)],'range')';
%     rippple_PSTH = [place_fields_all_L(1).V1_event_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).V1_event_PSTH_zscored(:,selected_cells_L);...
%         place_fields_all_R(1).V1_event_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).V1_event_PSTH_zscored(:,selected_cells_R)]';
else
    rippple_PSTH = normalize([place_fields_all_L(1).V1_event_PSTH(:,selected_cells_L);place_fields_all_L(2).V1_event_PSTH(:,selected_cells_L)],'range')';
%     rippple_PSTH = [place_fields_all_L(1).V1_event_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).V1_event_PSTH_zscored(:,selected_cells_L)]';
    rippple_PSTH = [rippple_PSTH rippple_PSTH];
end
rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_L(1).ripple_PSTH,1),[]);
subplot(2,2,1)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,1)));
hold on;xline(0,'r')
set(gca,'YDir','normal')
xlabel('Time')
title('L V1 bursting events on Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,2)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,2)));
hold on;xline(0,'r')
xlabel('Time')
title('L V1 bursting events on Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)
if length(selected_cells_L)==length(selected_cells_R)
    % rippple_PSTH = [place_fields_all_L(1).ripple_PSTH(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH(:,selected_cells_L)]';
        rippple_PSTH = normalize([place_fields_all_L(1).V1_event_PSTH(:,selected_cells_L);place_fields_all_L(2).V1_event_PSTH(:,selected_cells_L);...
            place_fields_all_R(1).V1_event_PSTH(:,selected_cells_R);place_fields_all_R(2).V1_event_PSTH(:,selected_cells_R)],'range')';
%     rippple_PSTH = [place_fields_all_L(1).V1_event_PSTH_zscored(:,selected_cells_L);place_fields_all_L(2).V1_event_PSTH_zscored(:,selected_cells_L);...
%     place_fields_all_R(1).V1_event_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).V1_event_PSTH_zscored(:,selected_cells_R)]';
else
    rippple_PSTH = normalize([place_fields_all_R(1).V1_event_PSTH(:,selected_cells_R);place_fields_all_R(2).V1_event_PSTH(:,selected_cells_R)],'range')';
%     rippple_PSTH = [place_fields_all_R(1).V1_event_PSTH_zscored(:,selected_cells_R);place_fields_all_R(2).V1_event_PSTH_zscored(:,selected_cells_R)]';
    rippple_PSTH = [rippple_PSTH rippple_PSTH];
end
% rippple_PSTH = [place_fields_all_R(1).ripple_PSTH(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH(:,selected_cells_R)]';
rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_R(1).ripple_PSTH,1),[]);
subplot(2,2,3)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,3)));
hold on;xline(0,'r')
set(gca,'YDir','normal')
xlabel('Time')
title('R V1 bursting events on Track Left')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

subplot(2,2,4)
imagesc(-1.5:0.02:1.5,[],squeeze(rippple_PSTH(:,25:175,4)));
hold on;xline(0,'r')
set(gca,'YDir','normal')
xlabel('Time')
title('R V1 bursting events on Track Right')
set(gca,"TickDir","out",'box', 'off','Color','none','FontSize',12)

if isfield(place_fields_all_L,'theta_phase_map') & ~isempty(place_fields_all_L(1).theta_phase_map)
    fig = figure
    fig.Position = [931 100 850 860]
    fig.Name = sprintf('All sessions theta modulation (%s)',cell_selection)
    sgtitle(sprintf('All sessions theta modulation (%s)',cell_selection))
    subplot(2,2,1)
    imagesc(0:20:360,[],normalize(reshape([place_fields_all_L(1).theta_phase_map{selected_cells_L}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
    set(gca,'YDir','normal')
    title('L theta on Track Left')
    subplot(2,2,2)
    imagesc(0:20:360,[],normalize(reshape([place_fields_all_L(2).theta_phase_map{selected_cells_L}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
    set(gca,'YDir','normal')
    title('L theta on Track Right')
    subplot(2,2,3)
    imagesc(0:20:360,[],normalize(reshape([place_fields_all_R(1).theta_phase_map{selected_cells_R}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
    set(gca,'YDir','normal')
    title('R theta on Track Left')
    subplot(2,2,4)
    imagesc(0:20:360,[],normalize(reshape([place_fields_all_R(2).theta_phase_map{selected_cells_R}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
    set(gca,'YDir','normal')
    title('R theta on Track Right')
    sgtitle('Theta')
end