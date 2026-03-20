addpath(genpath('/Users/atom/Documents/GitHub/V1_MEC_acute_MAT'))

mouse = {'M23032','M23034','M23037','M23038'};
date = {['20230718';'20230719';'20230720';'20230721';'20230722'];
    ['20230804';'20230805';'20230806';'20230807'];
    ['20230810';'20230811';'20230812';'20230813'];
    ['20230816';'20230817']};
session_count = 0;
V1_spatial_tuning = cell(15,4);
V1_SMI = cell(15,2);
HPC_spatial_tuning = cell(15,4);
MEC_spatial_tuning = cell(15,4);
HVA_spatial_tuning = cell(15,4);

for iMouse = 1:4
    for iDate = 1:size(date{iMouse,1},1)
        session_count = session_count+1;
        base_folder = '/Users/atom/Library/CloudStorage/OneDrive-UniversityCollegeLondon/data';
        session_folder = fullfile(base_folder,mouse{iMouse},'analysis',date{iMouse}(iDate,:));

        DIR = dir(session_folder);
        stimulus_list = {DIR(3:6).name};

        stimulus = 'Track';

        %
        load(fullfile(session_folder,stimulus,'extracted_behaviour.mat'))
        load(fullfile(session_folder,stimulus,'extracted_peripherals.mat'))
        load(fullfile(session_folder,stimulus,'extracted_task_info.mat'))
        load(fullfile(session_folder,stimulus,'session_info.mat'))
        load(fullfile(session_folder,'best_channels.mat'))
        load(fullfile(session_folder,stimulus,'extracted_place_fields.mat'))

        place_fields = spatial_modulation_calculation(place_fields,Task_info);

        %w = gausswin(21);
        w = gausswin(5);
        w = w / sum(w);
        no_unit = size(place_fields(1).raw,1);
        for i = 1:no_unit
            place_fields(1).raw{i} = filtfilt(w,1,place_fields(1).raw{i}')';
            place_fields(2).raw{i} = filtfilt(w,1,place_fields(2).raw{i}')';

        end
        t1_laps = 1:size(place_fields(1).raw{1},1);
        t2_laps = 1:size(place_fields(2).raw{1},1);
        
        t1_odd_laps = t1_laps(mod(t1_laps, 2) ~= 0);
        t1_even_laps = t1_laps(mod(t1_laps, 2) == 0);

        t2_odd_laps = t2_laps(mod(t2_laps, 2) ~= 0);
        t2_even_laps = t2_laps(mod(t2_laps, 2) == 0);

        putative_place_cells = cell(2,1);
        for i =1:2

            putative_place_cells{i} = find(place_fields(i).odd_even_stability >= 0.99 ...
                & place_fields(i).peak_percentile >= 0.99);

        end
        ANALYSIS_DATAPATH = session_info(1).probe(1).ANALYSIS_DATAPATH;
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH, 'Z:\ibn-vision\DATA\SUBJECTS', base_folder);
        ANALYSIS_DATAPATH = strrep(ANALYSIS_DATAPATH,'\', '/');
        session_info.probe(1).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
        session_info.probe(2).ANALYSIS_DATAPATH = ANALYSIS_DATAPATH;
        options = session_info.probe(2);
        V1_channels = determine_region_channels(best_channels{2},options,'region','V1','group','by probe');
        HPC_channels = determine_region_channels(best_channels{2},options,'region','HPC','group','by probe');
        options = session_info.probe(1);
        MEC_channels = determine_region_channels(best_channels{1},options,'region','MEC_entry','group','by probe');
        HVA_channels = determine_region_channels(best_channels{1},options,'region','HVA','group','by probe');
        %only cells with good spatial tuning
        V1_units = intersect(find(ismember(place_fields(1).peak_channel,V1_channels)),find(place_fields(1).cluster_id<10000));
        HPC_units = intersect(find(ismember(place_fields(1).peak_channel,HPC_channels)),find(place_fields(1).cluster_id<10000));
        MEC_units = intersect(find(ismember(place_fields(1).peak_channel,MEC_channels)),find(place_fields(1).cluster_id>10000));
        HVA_units = intersect(find(ismember(place_fields(1).peak_channel,HVA_channels)),find(place_fields(1).cluster_id>10000));

        V1_spatial_tuning{session_count,1} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_odd_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',V1_units)),'UniformOutput',false),'UniformOutput',false));
        V1_spatial_tuning{session_count,2} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_even_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',V1_units)),'UniformOutput',false),'UniformOutput',false));
        V1_spatial_tuning{session_count,3} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_odd_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',V1_units)),'UniformOutput',false),'UniformOutput',false));
        V1_spatial_tuning{session_count,4} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_even_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',V1_units)),'UniformOutput',false),'UniformOutput',false));



        HPC_spatial_tuning{session_count,1} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_odd_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',HPC_units)),'UniformOutput',false),'UniformOutput',false));
        HPC_spatial_tuning{session_count,2} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_even_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',HPC_units)),'UniformOutput',false),'UniformOutput',false));
        HPC_spatial_tuning{session_count,3} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_odd_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',HPC_units)),'UniformOutput',false),'UniformOutput',false));
        HPC_spatial_tuning{session_count,4} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_even_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',HPC_units)),'UniformOutput',false),'UniformOutput',false));



        MEC_spatial_tuning{session_count,1} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_odd_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',MEC_units)),'UniformOutput',false),'UniformOutput',false));
        MEC_spatial_tuning{session_count,2} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_even_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',MEC_units)),'UniformOutput',false),'UniformOutput',false));
        MEC_spatial_tuning{session_count,3} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_odd_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',MEC_units)),'UniformOutput',false),'UniformOutput',false));
        MEC_spatial_tuning{session_count,4} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_even_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',MEC_units)),'UniformOutput',false),'UniformOutput',false));
    

        HVA_spatial_tuning{session_count,1} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_odd_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',HVA_units)),'UniformOutput',false),'UniformOutput',false));
        HVA_spatial_tuning{session_count,2} = cell2mat(cellfun(@mean,cellfun(@(x) x(t1_even_laps,:),place_fields(1).raw(intersect(putative_place_cells{1}',HVA_units)),'UniformOutput',false),'UniformOutput',false));
        HVA_spatial_tuning{session_count,3} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_odd_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',HVA_units)),'UniformOutput',false),'UniformOutput',false));
        HVA_spatial_tuning{session_count,4} = cell2mat(cellfun(@mean,cellfun(@(x) x(t2_even_laps,:),place_fields(2).raw(intersect(putative_place_cells{2}',HVA_units)),'UniformOutput',false),'UniformOutput',false));


        V1_SMI{session_count,2} = sort(place_fields(2).SMI(V1_units));
        V1_SMI{session_count,1} = sort(place_fields(1).SMI(V1_units));
    end
end

MEC_across_t1_odd = normalize(cell2mat(MEC_spatial_tuning(:,1)),2,'range');
MEC_across_t1_even = normalize(cell2mat(MEC_spatial_tuning(:,2)),2,'range');
[~,MEC_t1_peak_odd] = max(MEC_across_t1_odd,[],2);
[~,MEC_t1_peak_odd_sort] = sort(MEC_t1_peak_odd);
[~,MEC_t1_peak_even] = max(MEC_across_t1_even,[],2);
[~,MEC_t1_peak_even_sort] = sort(MEC_t1_peak_even);


MEC_across_t2_odd = normalize(cell2mat(MEC_spatial_tuning(:,3)),2,'range');
MEC_across_t2_even = normalize(cell2mat(MEC_spatial_tuning(:,4)),2,'range');
[~,MEC_t2_peak_odd] = max(MEC_across_t2_odd,[],2);
[~,MEC_t2_peak_odd_sort] = sort(MEC_t2_peak_odd);
[~,MEC_t2_peak_even] = max(MEC_across_t2_even,[],2);
[~,MEC_t2_peak_even_sort] = sort(MEC_t2_peak_even);
figure;subplot(2,2,1)
imagesc(place_fields(1).x_bin_centres,1:length(MEC_t1_peak_odd_sort),flip(MEC_across_t1_even(MEC_t1_peak_odd_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('MEC t1 even sorted by odd')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,2)
imagesc(place_fields(1).x_bin_centres,1:length(MEC_t1_peak_even_sort),flip(MEC_across_t1_odd(MEC_t1_peak_even_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('MEC t1 odd sorted by even')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,3)
imagesc(place_fields(1).x_bin_centres,1:length(MEC_t2_peak_odd_sort),flip(MEC_across_t2_even(MEC_t2_peak_odd_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('MEC t2 even sorted by odd')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
subplot(2,2,4)
imagesc(place_fields(1).x_bin_centres,1:length(MEC_t2_peak_even_sort),flip(MEC_across_t2_odd(MEC_t2_peak_even_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('MEC t2 odd sorted by even')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

V1_across_t1_odd = normalize(cell2mat(V1_spatial_tuning(:,1)),2,'range');
V1_across_t1_even = normalize(cell2mat(V1_spatial_tuning(:,2)),2,'range');
[~,V1_t1_peak_odd] = max(V1_across_t1_odd,[],2);
[~,V1_t1_peak_odd_sort] = sort(V1_t1_peak_odd);
[~,V1_t1_peak_even] = max(V1_across_t1_even,[],2);
[~,V1_t1_peak_even_sort] = sort(V1_t1_peak_even);


V1_across_t2_odd = normalize(cell2mat(V1_spatial_tuning(:,3)),2,'range');
V1_across_t2_even = normalize(cell2mat(V1_spatial_tuning(:,4)),2,'range');
[~,V1_t2_peak_odd] = max(V1_across_t2_odd,[],2);
[~,V1_t2_peak_odd_sort] = sort(V1_t2_peak_odd);
[~,V1_t2_peak_even] = max(V1_across_t2_even,[],2);
[~,V1_t2_peak_even_sort] = sort(V1_t2_peak_even);
figure;subplot(2,2,1)
imagesc(place_fields(1).x_bin_centres(1,20:60),1:length(V1_t1_peak_odd_sort),flip(V1_across_t1_even(V1_t1_peak_odd_sort,20:60)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t1 even sorted by odd')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,2)
imagesc(place_fields(1).x_bin_centres(1,20:60),1:length(V1_t1_peak_even_sort),flip(V1_across_t1_odd(V1_t1_peak_even_sort,20:60)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t1 odd sorted by even')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,3)
imagesc(place_fields(1).x_bin_centres(1,20:60),1:length(V1_t2_peak_odd_sort),flip(V1_across_t2_even(V1_t2_peak_odd_sort,20:60)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t2 even sorted by odd')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,4)
imagesc(place_fields(1).x_bin_centres(1,20:60),1:length(V1_t2_peak_even_sort),flip(V1_across_t2_odd(V1_t2_peak_even_sort,20:60)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t2 odd sorted by even')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)


figure;subplot(2,2,1)
imagesc(place_fields(1).x_bin_centres(1,:),1:length(V1_t1_peak_odd_sort),flip(V1_across_t1_even(V1_t1_peak_odd_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t1 even sorted by odd')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,2)
imagesc(place_fields(1).x_bin_centres(1,:),1:length(V1_t1_peak_even_sort),flip(V1_across_t1_odd(V1_t1_peak_even_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t1 odd sorted by even')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,3)
imagesc(place_fields(1).x_bin_centres(1,:),1:length(V1_t2_peak_odd_sort),flip(V1_across_t2_even(V1_t2_peak_odd_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t2 even sorted by odd')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)

subplot(2,2,4)
imagesc(place_fields(1).x_bin_centres(1,:),1:length(V1_t2_peak_even_sort),flip(V1_across_t2_odd(V1_t2_peak_even_sort,:)))
clim([0.25 0.75]);colormap(flip(gray(50)));colorbar
title('V1 t2 odd sorted by even')
xlabel('position(cm)')
ylabel('units')
set(gca,'TickDir','out','box','off','Color','none','FontSize',12)
V1_SMI_t1 = cell2mat(V1_SMI(:,1));
V1_SMI_t1 = V1_SMI_t1(V1_t1_peak_odd<=60 & V1_t1_peak_odd >=20);
V1_SMI_t1 = V1_SMI_t1(~isnan(V1_SMI_t1));
V1_SMI_t2 = cell2mat(V1_SMI(:,2));
V1_SMI_t2 = V1_SMI_t2(V1_t2_peak_odd<=60 & V1_t2_peak_odd >=20);
V1_SMI_t2 = V1_SMI_t2(~isnan(V1_SMI_t2));


[V1_SMI_t1_edge, V1_SMI_t1_hist] = hist(V1_SMI_t1, 100);
[V1_SMI_t2_edge, V1_SMI_t2_hist] = hist(V1_SMI_t2, 100);
figure;subplot(1,2,1)
plot(V1_SMI_t1_hist,cumsum(V1_SMI_t1_edge)/length(V1_SMI_t1))
hold on; xline(0);yline(0.5)
subplot(1,2,2)
plot(V1_SMI_t2_hist,cumsum(V1_SMI_t2_edge)/length(V1_SMI_t2))
hold on; xline(0);yline(0.5)
