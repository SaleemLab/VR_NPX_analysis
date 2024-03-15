function plot_spatial_theta_ripple_population(selected_cells_L,selected_cells_R,place_fields_all_L,place_fields_all_R)

fig = figure
average_map = normalize([place_fields_all_L(1).average_map(:,selected_cells_L);place_fields_all_L(2).average_map(:,selected_cells_L)],'range')';
average_map = reshape(average_map,size(average_map,1),size(place_fields_all_L(1).average_map,1),[]);
subplot(2,2,1)
imagesc(squeeze(average_map(:,:,1)));
xlabel('position')
title('Track Left')

subplot(2,2,2)
imagesc(squeeze(average_map(:,:,2)));
xlabel('position')
title('Track Right')

average_map = normalize([place_fields_all_R(1).average_map(:,selected_cells_R);place_fields_all_R(2).average_map(:,selected_cells_R)],'range')';
average_map = reshape(average_map,size(average_map,1),size(place_fields_all_R(1).average_map,1),[]);
subplot(2,2,3)
imagesc(squeeze(average_map(:,:,1)));
xlabel('position')
title('Track Left')

subplot(2,2,4)
imagesc(squeeze(average_map(:,:,2)));
xlabel('position')
title('Track Right')


figure
rippple_PSTH = normalize([place_fields_all_L(1).ripple_PSTH(:,selected_cells_L);place_fields_all_L(2).ripple_PSTH(:,selected_cells_L)],'range')';
rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_L(1).ripple_PSTH,1),[]);
subplot(2,2,1)
imagesc(squeeze(rippple_PSTH(:,:,1)));
xline(50)
xlabel('Time')
title('L ripple on Track Left')

subplot(2,2,2)
imagesc(squeeze(rippple_PSTH(:,:,2)));
xline(50)
xlabel('position')
title('L ripple on Track Right')

rippple_PSTH = normalize([place_fields_all_R(1).ripple_PSTH(:,selected_cells_R);place_fields_all_R(2).ripple_PSTH(:,selected_cells_R)],'range')';
rippple_PSTH = reshape(rippple_PSTH,size(rippple_PSTH,1),size(place_fields_all_R(1).ripple_PSTH,1),[]);
subplot(2,2,3)
imagesc(squeeze(rippple_PSTH(:,:,1)));
xline(50)
xlabel('Time')
title('R ripple on Track Left')

subplot(2,2,4)
imagesc(squeeze(rippple_PSTH(:,:,2)));
xline(50)
xlabel('position')
title('R ripple on Track Right')


figure
subplot(2,2,1)
imagesc(normalize(reshape([place_fields_all_L(1).theta_phase_map{selected_cells_L}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
title('L theta on Track Left')
subplot(2,2,2)
imagesc(normalize(reshape([place_fields_all_L(2).theta_phase_map{selected_cells_L}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
title('L theta on Track Right')
subplot(2,2,3)
imagesc(normalize(reshape([place_fields_all_R(1).theta_phase_map{selected_cells_R}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
title('R theta on Track Left')
subplot(2,2,4)
imagesc(normalize(reshape([place_fields_all_R(2).theta_phase_map{selected_cells_R}],size(place_fields_all_R(1).theta_phase_map{1},2),[]),'range')');
title('R theta on Track Right')
sgtitle('Theta')