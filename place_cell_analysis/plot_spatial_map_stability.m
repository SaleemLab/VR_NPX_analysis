function plot_spatial_map_stability(place_fields_all,options) 

all_regions={'V1_L','V1_R'};

for n = 1:2
    ia = find((place_fields_all(1).odd_even_stability>0.95 ...
        | place_fields_all(2).odd_even_stability>0.95)&contains(place_fields_all(1).region,all_regions{n}));

    % ia = find(place_fields_all(1).odd_even_stability>0.95 ...
    %     | place_fields_all(2).odd_even_stability>0.95);
    ratemap_matrix=[];
    average_maps=[];
    lap_correlation_pv=[];
    within_track_corr=[];
    lap_correlation=[];
    pv_corr=[];
    for track_id = 1:2
        %             ratemap_matrix{track_id} = [conv(place_fields_all(track_id).raw{good_cell_index}gaussianWindow,'same')];
        ratemap_matrix{track_id} = [place_fields_all(track_id).raw{ia}];
        ratemap_matrix{track_id} = reshape( ratemap_matrix{track_id},size(place_fields_all(track_id).within_track_corr{1},1),[],length(ia));%laps X position bins X cells


        within_track_corr{track_id} = [place_fields_all(track_id).within_track_corr{ia}];
        %             ratemaps_track1(nevent,:) = conv(ratemaps_track1(nevent,:),gaussianWindow,'same');
        within_track_corr{track_id} = reshape( within_track_corr{track_id},size(place_fields_all(track_id).within_track_corr{1},1),size(place_fields_all(track_id).within_track_corr{1},1),[]);
        lap_correlation{track_id}= squeeze(nanmean(within_track_corr{track_id},3));

        for nlap = 1:size(place_fields_all(track_id).raw{1},1)
            for mlap = 1:size(place_fields_all(track_id).raw{1},1)
                map1 = zscore(squeeze(ratemap_matrix{track_id}(nlap,:,:)));
                map2 = zscore(squeeze(ratemap_matrix{track_id}(mlap,:,:)));
                temp = corr(map1',map2');
                x_bins = 1:size(map1,1);

                lap_correlation_pv{track_id}(nlap,mlap) = nanmean(diag(temp));
                %                     lap_correlation_pv{track_id}(nlap,mlap) = nanmean(nanmean(corr(map1,map2)));
            end
        end

        for nlap = 1:size(place_fields_all(track_id).raw{1},1)
            map1 = zscore(squeeze(ratemap_matrix{track_id}(nlap,:,:)));
            map2 = zscore(squeeze(nanmean(ratemap_matrix{track_id}(:,:,:))));
            temp = corr(map1',map2');
            x_bins = 1:size(map1,1);

            pv_corr(track_id,nlap) = nanmean(diag(temp));
        end
    end

    fig = figure
    fig.Position = [500 70 830 600];
    sgtitle(sprintf('%s %s %s cell spatial map stability probe %s',options.SUBJECT,options.SESSION,all_regions{n},place_fields_all(1).probe_hemisphere(1)))
    fig.Name = sprintf('%s %s %s cell spatial map stability probe %s',options.SUBJECT,options.SESSION,all_regions{n},place_fields_all(1).probe_hemisphere(1));

    subplot(3,2,1)
    percentile_95 = prctile(reshape(lap_correlation{1},1,[]),95);
    imagesc(lap_correlation{1}')
    clim([0 percentile_95])
    title('Track 1 mean lap by lap single cell corr')
    colorbar

    subplot(3,2,2)
    percentile_95 = prctile(reshape(lap_correlation{2},1,[]),95);
    imagesc(lap_correlation{2}')
    clim([0 percentile_95])
    colorbar
    title('Track 2 mean lap by lap single cell corr')


    subplot(3,2,3)
    percentile_95 = prctile(reshape(lap_correlation_pv{1},1,[]),95);
    imagesc(lap_correlation_pv{1}')
    clim([0 percentile_95])
    title('Track 1 mean lap to lap PV corr')
    colorbar

    subplot(3,2,4)
    percentile_95 = prctile(reshape(lap_correlation_pv{2},1,[]),95);
    imagesc(lap_correlation_pv{2}')
    clim([0 percentile_95])
    colorbar
    title('Track 2 mean lap to lap PV corr')

    subplot(3,2,5)
    %         percentile_95 = prctile(reshape(pv_corr{1},1,[]),95);
    plot(movmedian(pv_corr(1,:)',3))
    %       ylim([-percentile_95 percentile_95])
    title('Track 1 mean lap to template PV corr')
    box off
    %         colorbar

    subplot(3,2,6)
    %         percentile_95 = prctile(reshape(pv_corr(2,:),1,[]),95);
    plot(movmedian(pv_corr(2,:)',3))
    %         ylim([-percentile_95 percentile_95])
    title('Track 2 mean lap to template PV corr')
    box off
end

