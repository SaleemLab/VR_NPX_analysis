    
function [laps_corr,singleCell_corr] = place_field_stability_masa(clusters,place_fields,position,lap_times,tracks_compared,number_of_laps)
% Marta Huelin_ 2019
% Assess place field stability between same track exposures or between
% different tracks. It does so both at population level and at single cell level. 
% It looks at the correlation and difference (respectively) of the following parameters at population level: centre of mass, peak and mean firing rate;
% INPUTS: 
    % tracks_compared input options are: 
         % within_track: compare laps within the same track exposure
         % between_exposures:compare laps between first and second exposure to the tracks
         % between_exposures_FULLT1: compares full ratemap of first exposure (T1/T2) to each lap of second exposure
         % between_exposures_REexp: compare the last laps from second exposure to laps in first exposure
    % number_of_laps: number of laps that are compared 

%     load('extracted_clusters.mat');
%     load('extracted_position.mat');
%     load('extracted_laps.mat');
%     load('extracted_place_fields.mat')
       
   % Set tracks that need to be compared and type of comparisons based on the number of laps 
    if strcmp(tracks_compared,'within_track')
        compare_trackIDs = {[1,1], [2,2]};  
%         compare_trackIDs = {[1,1], [2,2], [3,3], [4,4]};  
        count = 1;
        consecutive = 1;
        % Set lap indices to be compared
        for track = 1: length(compare_trackIDs)
            % For comparing consecutive laps (e.g. laps 1 to 6 vs laps 7 to 12):
            c = 1;
            for i = 1:length(lap_times(track).total_number_of_laps)
                if track == 1 || track == 3 || track == 4      %&& lap_times(track).number_completeLaps<16
                    if i <= 12 %do comparisons only up to lap 12
                        consecutive_laps(track).lap_ids(c,1) = {[i i+number_of_laps-1]};  % start and end lap ID for the first set of laps
                        consecutive_laps(track).lap_ids(c,2) = {[i+number_of_laps i+2*number_of_laps-1]}; % start and end lap ID for the second set of laps
                    end
                else
                    consecutive_laps(track).lap_ids(c,1) = {[i i+number_of_laps-1]};  % start and end lap ID for the first set of laps
                    consecutive_laps(track).lap_ids(c,2) = {[i+number_of_laps i+2*number_of_laps-1]}; % start and end lap ID for the second set of laps
                end
                c = c +1;
            end
            % For comparing sets of laps vs final laps of the exposure (e.g. 1-3 vs 9-12; 2-4 vs 9-12; etc):
            % Added a change, where T1,T3 & T4 are now comparing to last 2-4 laps(depending if run 14 or 16 laps)
            c = 1;
            for i = lap_times(track).completeLap_id(1):1:(lap_times(track).completeLap_id(end)-number_of_laps)
                if track == 1 %&& lap_times(track).number_completeLaps<16
                    if i <= 12 %do comparisons only up to lap 12
                        ends_laps(track).lap_ids(c,1) = {[i i+number_of_laps-1]};  % start and end lap ID for the first set of laps
                        if lap_times(track).number_completeLaps >= 16
                            end_lap_ID = 16;
                        elseif lap_times(track).number_completeLaps < 16
                            end_lap_ID = lap_times(track).completeLap_id(end);
                        end
                        ends_laps(track).lap_ids(c,2) = {[13 end_lap_ID]}; % start and end lap ID for the second set of laps
                    end
                elseif track == 3 || track == 4 %&& lap_times(track).number_completeLaps>15
                    if i <= 12 %do comparisons only up to lap 12
                        ends_laps(track).lap_ids(c,1) = {[i i+number_of_laps-1]};  % start and end lap ID for the first set of laps
                        ends_laps(track).lap_ids(c,2) = {[13 16]}; % start and end lap ID for the second set of laps
                    end
                else %regular (and original) way
                    ends_laps(track).lap_ids(c,1) = {[i i+number_of_laps-1]};  % start and end lap ID for the first set of laps
                    ends_laps(track).lap_ids(c,2) = {[lap_times(track).completeLap_id(end)-number_of_laps+1 lap_times(track).completeLap_id(end)]}; % start and end lap ID for the second set of laps
                end
                c = c +1;
            end
        end
    elseif strcmp(tracks_compared,'between_exposures')
        compare_trackIDs = {[1,3], [2,4]};
        consecutive = 0;
        count = 1;
        % Compare last 4 laps of first exposure to each lap of second exposure (starting from lap #1)
        for track = 1 : length(compare_trackIDs)
            c = 1;
            for i = lap_times(track+2).completeLap_id(1):number_of_laps:lap_times(track+2).completeLap_id(end)
                if track == 1
                    if lap_times(track).number_completeLaps >= 16
                        end_lap_ID = 16;
                    elseif lap_times(track).number_completeLaps < 16
                        end_lap_ID = lap_times(track).completeLap_id(end);
                    end
                    ends_laps(track).lap_ids(c,1) = {[13 end_lap_ID]};  % start and end lap ID for track 1 (last 2/4 laps depends if total is 14 or 16 laps)
                    ends_laps(track).lap_ids(c,2) = {[i i+number_of_laps-1]}; % start and end lap ID for track 3 
                else
                    if lap_times(track).number_completeLaps >= 8
                        start_lap = lap_times(track).completeLap_id(end-1);
                        end_lap_ID = 8;
                    elseif lap_times(track).number_completeLaps == 1
                        end_lap_ID = 1;
                        start_lap = 1;
                    else 
                        end_lap_ID = lap_times(track).completeLap_id(end);
                        start_lap = lap_times(track).completeLap_id(end-1);
                    end
                    ends_laps(track).lap_ids(c,1) = {[start_lap end_lap_ID]};  % start and end lap ID for track 2 (last lap)
                    ends_laps(track).lap_ids(c,2) = {[i i+number_of_laps-1]}; % start and end lap ID for track 4
                end
                c = c +1;
            end
        end
    elseif strcmp(tracks_compared,'between_exposures_FULLT1')
        compare_trackIDs = {[1,3], [2,4]};
        consecutive = 0;
        count = 1;
        % Compare spatial map of first exposure to each lap of second exposure (starting from lap #1)
        for track = 1 : length(compare_trackIDs)
            c = 1;
            for i = lap_times(track+2).completeLap_id(1):number_of_laps:lap_times(track+2).completeLap_id(end)
                if track == 1
                    ends_laps(track).lap_ids(c,2) = {[i i+number_of_laps-1]}; % start and end lap ID for track 3
                else
                    ends_laps(track).lap_ids(c,2) = {[i i+number_of_laps-1]}; % start and end lap ID for track 4
                end
                c = c +1;
            end
        end
    elseif strcmp(tracks_compared,'between_exposures_REexp')
        compare_trackIDs = {[1,3], [2,4]};
        consecutive = 0;
        count = 1;
        % Compare laps 12-16 from second exposure to first exposure
        for track = 1 : length(compare_trackIDs)
            c = 1;
            for i = lap_times(track).completeLap_id(1):number_of_laps:lap_times(track).completeLap_id(end)
                if track == 1
                    if i+number_of_laps-1 > 16
                        continue
                    else
                        end_lap_ID = i+number_of_laps-1;
                    end
                    ends_laps(track).lap_ids(c,1) = {[i i+number_of_laps-1]};  % start and end lap ID for track 1 (last 2/4 laps depends if total is 14 or 16 laps)
                    ends_laps(track).lap_ids(c,2) = {[13 16]}; % start and end lap ID for track 3
                else
                    if i+number_of_laps-1 > 8
                        continue
                    else
                        end_lap_ID = i+number_of_laps-1;
                    end
                    ends_laps(track).lap_ids(c,1) = {[i end_lap_ID]};  % start and end lap ID for track 2 (last lap)
                    ends_laps(track).lap_ids(c,2) = {[13 16]}; % start and end lap ID for track 4
                end
                c = c +1;
            end
        end
    end
    
   % Start correlation analysis for each type of comparison     
   
   while count == 1   %If both consecutive laps and end laps are being analysed
       
        if consecutive == 1 % If consecutive laps are also being analysed
            laps_compared = consecutive_laps;
            comparison_type = 'consecutive_laps';
        else
            laps_compared = ends_laps;
            comparison_type = 'ends_laps';          
        end
        
        for t = 1 : length(compare_trackIDs)
            track_id = cell2mat(compare_trackIDs(t));          
            
            if strcmp(tracks_compared,'within_track')&& consecutive == 0 % to avoid overwriting when running two comparisons
                j = t + length(compare_trackIDs);
            else
                j = t;
            end 
            
            singleCell_corr(j).tracks_compared = mat2str(track_id);  % Saves loop info in single cell correlation struct
            singleCell_corr(j).comparison_type = strcat(tracks_compared,'-',comparison_type);
            singleCell_corr(j).laps_jump = number_of_laps;

            for i = 1 : size(laps_compared(t).lap_ids,1)
                
                laps1 = cell2mat(laps_compared(t).lap_ids(i,1));
                laps2 = cell2mat(laps_compared(t).lap_ids(i,2));
                
                % Calculates place field inside lap(s) selected for laps 1 and laps 2
                if  strcmp(tracks_compared,'between_exposures_FULLT1')
                    place_fields1 = place_fields.track(track_id(1)); %full rate map of T1/T2
                    place_fields1.mean_rate_lap = place_fields1.mean_rate_track;
                    place_fields2 = get_lap_place_fields(track_id(2),laps2(1),laps2(2),0,'complete');%Laps 2
                else
                    place_fields1 = get_lap_place_fields(track_id(1),laps1(1),laps1(2),0,'complete'); %Laps 1
                    place_fields2 = get_lap_place_fields(track_id(2),laps2(1),laps2(2),0,'complete');%Laps 2
                end
                
                % Find all good cells active in both set of laps and save those ones that are only active in one of the set of laps
                common_goodCells = intersect(place_fields1.good_cells,place_fields2.good_cells);
                common_activeCells = intersect(find(~isnan(place_fields1.centre_of_mass)),find(~isnan(place_fields2.centre_of_mass))); % if a centre of mass is NaN it means that the cell wasn't active in that set of laps
                common_lapCells_indices = intersect(common_goodCells,common_activeCells);
                
                % Save information in structure
                laps_corr(j).tracks_compared = mat2str(track_id);
                laps_corr(j).comparison_type = strcat(tracks_compared,'-',comparison_type);
                laps_corr(j).laps_jump = number_of_laps;
                laps_corr(j).LapsID_1(i) = {laps1};
                laps_corr(j).LapsID_2(i) = {laps2}; 
                laps_corr(j).pyr_cells_Laps1(i) = {place_fields1.good_cells};
                laps_corr(j).pyr_cells_Laps2(i) = {place_fields2.good_cells};
                laps_corr(j).cells_analysed(i) = {common_lapCells_indices};
                laps_corr(j).cells_active_only_1(i) = {setdiff(place_fields1.good_cells,common_lapCells_indices)};
                laps_corr(j).cells_active_only_2(i) = {setdiff(place_fields2.good_cells,common_lapCells_indices)};
                
                % Calculate correlation between center of mass
                [laps_corr(j).centremass_r(i),laps_corr(j).centremass_p(i)] = corr(place_fields1.centre_of_mass(common_lapCells_indices)',place_fields2.centre_of_mass(common_lapCells_indices)','type','Spearman','Rows','pairwise'); %assumes non normal distribution
                
                % Calculate correlation between normalized center of mass
                [laps_corr(j).normcentremass_r(i),laps_corr(j).normcentremass_p(i)] = corr((place_fields1.centre_of_mass(common_lapCells_indices)./place_fields1.half_max_width(common_lapCells_indices))',...
                    (place_fields2.centre_of_mass(common_lapCells_indices)./place_fields2.half_max_width(common_lapCells_indices))','type','Spearman','Rows','pairwise'); %assumes non normal distribution
                
                % Calculate correlation between peak of place field
                [laps_corr(j).peak_r(i),laps_corr(j).peak_p(i)] = corr(place_fields1.peak(common_lapCells_indices)', place_fields2.peak(common_lapCells_indices)','type','Spearman','Rows','pairwise');
                
                % Find the timestamps from chosen lap_start and lap_End
                if  strcmp(tracks_compared,'between_exposures_FULLT1')
                    timeWindow_Laps1 = [place_fields1.time_window(1) place_fields1.time_window(2)]; %full rate map of T1/T2
                    timeWindow_Laps2 = [lap_times(track_id(2)).completeLaps_start(laps2(2)==lap_times(track_id(2)).completeLap_id), lap_times(track_id(2)).completeLaps_stop(laps2(2)==lap_times(track_id(2)).completeLap_id)];
                else
                    timeWindow_Laps1 = [lap_times(track_id(1)).completeLaps_start(laps1(1)==lap_times(track_id(1)).completeLap_id), lap_times(track_id(1)).completeLaps_stop(laps1(1)==lap_times(track_id(1)).completeLap_id)];
                    timeWindow_Laps2 = [lap_times(track_id(2)).completeLaps_start(laps2(2)==lap_times(track_id(2)).completeLap_id), lap_times(track_id(2)).completeLaps_stop(laps2(2)==lap_times(track_id(2)).completeLap_id)];
                end
                
                % For each common good cell, calculate difference between centre of mass, peak and mean firing rate (includes ones active only in 1 of the laps set)
                for kk = 1:length(common_goodCells)
                                                            
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).cell= common_goodCells(kk); % cluster_id - To know original cluster_id, need to check clusters.id_conversion
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).LapsID_1 = laps1;
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).LapsID_2= laps2; 
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).laps1_duration = timeWindow_Laps1(2)- timeWindow_Laps1(1); %time window analyzed for first set of laps
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).laps2_duration = timeWindow_Laps2(2)- timeWindow_Laps2(1); %time window analyzed for second set of laps
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).peakFR = [place_fields1.peak(common_goodCells(kk)),place_fields2.peak(common_goodCells(kk))]; % peak firing rate for laps 1 and laps 2
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).meanFR = [place_fields1.mean_rate_lap(common_goodCells(kk)),place_fields2.mean_rate_lap(common_goodCells(kk))]; % mean firing rate for laps 1 and laps 2
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).centreMass = [place_fields1.centre_of_mass(common_goodCells(kk)),place_fields2.centre_of_mass(common_goodCells(kk))]; %centre of mass for laps 1 and laps 2
                    
                    if kk < length(common_goodCells) && consecutive == 0
                        all_peakFR(common_goodCells(kk),i) = place_fields1.peak(common_goodCells(kk)); %Rows - cells/ columns - laps
                        all_meanFR(common_goodCells(kk),i) = place_fields1.mean_rate_lap(common_goodCells(kk));
                        all_centreMass(common_goodCells(kk),i)= place_fields1.centre_of_mass(common_goodCells(kk));
                    elseif kk == length(common_goodCells) && consecutive == 0
                        all_peakFR(common_goodCells(kk),i) = place_fields2.peak(common_goodCells(kk));
                        all_meanFR(common_goodCells(kk),i) = place_fields2.mean_rate_lap(common_goodCells(kk));
                        all_centreMass(common_goodCells(kk),i) = place_fields2.centre_of_mass(common_goodCells(kk));
                    else
                        all_peakFR(common_goodCells(kk),i) = place_fields1.peak(common_goodCells(kk));
                        all_meanFR(common_goodCells(kk),i) = place_fields1.mean_rate_lap(common_goodCells(kk));
                        all_centreMass(common_goodCells(kk),i) = place_fields1.centre_of_mass(common_goodCells(kk));                        
                    end
                    
                    cell_spikeTimes = clusters.spike_times(clusters.spike_id == common_goodCells(kk));
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).numSpikes_Laps1 = length(find(cell_spikeTimes>timeWindow_Laps1(1) & cell_spikeTimes<timeWindow_Laps1(2))); %number of spikes for laps set 1
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).numSpikes_Laps2 = length(find(cell_spikeTimes>timeWindow_Laps2(1) & cell_spikeTimes<timeWindow_Laps2(2))); %number of spikes for laps set 2
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).numSpikes_diff = singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).numSpikes_Laps2 - singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).numSpikes_Laps1;
                    
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).peakFR_diff = place_fields2.peak(common_goodCells(kk)) - place_fields1.peak(common_goodCells(kk)); % peak firing rate difference
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).meanFR_diff = place_fields2.mean_rate_lap(common_goodCells(kk)) - place_fields1.mean_rate_lap(common_goodCells(kk)); % mean firing rate difference
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).centreMass_diff = place_fields2.centre_of_mass(common_goodCells(kk)) - place_fields1.centre_of_mass(common_goodCells(kk)); %centre of mass difference
                    singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).norm_centreMass_diff = (place_fields2.centre_of_mass(common_goodCells(kk)) - place_fields1.centre_of_mass(common_goodCells(kk)))/(place_fields1.half_max_width(common_goodCells(kk))+place_fields2.half_max_width(common_goodCells(kk))); %centre of mass difference
                    
                    [singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).xcorr_raw_Plfields , singleCell_corr(j).(sprintf('%s','Laps_',num2str(i)))(kk).xcorr_lags] = xcorr(place_fields1.raw{common_goodCells(kk)},place_fields2.raw{common_goodCells(kk)}); 

                end
            end
            
            % For each cell, linear regression between number of laps and peak/mean firing rate/centre of mass
            
            for cell = 1 : length(all_centreMass)
                if ~isempty(all_centreMass(cell,:))
                    laps_corr(j).singleCell_laps_corr(cell).cell = cell;
                    [laps_corr(j).singleCell_laps_corr(cell).centreMass_LR , laps_corr(j).singleCell_laps_corr(cell).centreMass_slope,~] = regression([1:1:length(all_centreMass(cell,:))],all_centreMass(cell,:));
                    [laps_corr(j).singleCell_laps_corr(cell).peakFR_LR , laps_corr(j).singleCell_laps_corr(cell).peakFR_slope,~] = regression([1:1:length(all_peakFR(cell,:))],all_peakFR(cell,:));
                    [laps_corr(j).singleCell_laps_corr(cell).meanFR_LR , laps_corr(j).singleCell_laps_corr(cell).meanFR_slope,~] = regression([1:1:length(all_meanFR(cell,:))],all_meanFR(cell,:));
                end
            end
        end
        
        % Update count and comparison_type for the next loop
        if consecutive == 1
            consecutive = 0;
            count = 1;
        else
            count = 0;
        end      
   end
   
   if strcmp(tracks_compared,'within_track')
       save plfield_LapCompare_WithinTrack laps_corr singleCell_corr
   elseif strcmp(tracks_compared,'between_exposures')
       save plfield_LapCompare_BetweenExposures laps_corr singleCell_corr
   elseif strcmp(tracks_compared,'between_exposures_REexp')
       save plfield_LapCompare_BetweenExposures_REexp laps_corr singleCell_corr
   elseif strcmp(tracks_compared,'between_exposures_FULLT1')
       save plfield_LapCompare_BetweenExposures_FULLT1 laps_corr singleCell_corr
   end

end