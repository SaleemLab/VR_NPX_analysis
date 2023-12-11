function estimated_position_test = bayesian_decoding_error_masa(varargin)
% INPUTS [ALL OPTIONAL, will run defaults if no inputs]: 
%               input as argument ('name', variable) pairs
%               'method': 'leave one out', 'cross_tracks'
%               'control': 'global_remap', 'rate_remap'

% this code extracts RUN periods (v > 5cm/s) during laps (back and forth) 
% then creates a training and test set for the bayesian decoder
% creates cumulative distribution curves and confusion matrices
parameters= list_of_parameters;

p = inputParser;
addParameter(p,'method','leave one out',@ischar);
addParameter(p,'control',[],@ischar);
addParameter(p,'bin_size',parameters.x_bins_width_bayesian,@ismatrix);
addParameter(p,'save_figures',0,@ismatrix);
addParameter(p,'place_fields',[],@isstruct);

parse(p,varargin{:});
% makes variables a bit shorter to call...
method= p.Results.method;
control= p.Results.control;
bin_size= p.Results.bin_size;
place_fields = p.Results.bin_size;

switch method
    case 'leave one out'
        disp('running leave one out')
        load('lap_times.mat');
%         if p.Results.bin_size == parameters.x_bins_width_bayesian
%             load('extracted_place_fields_BAYESIAN.mat');
%             place_fields= place_fields_BAYESIAN;
%             clear place_fields_BAYESIAN;
%         elseif p.Results.bin_size == parameters.x_bins_width
%             load('extracted_place_fields.mat');
%         end
        disp('finding laps...')
        timestep= mean(diff(lap_times(1).lap(1).t));
        for this_track=1:length(lap_times)
            typ= mod(length(lap_times(this_track).lap),2);
            if typ % odd number of laps
                start_idx= 1:2:length(lap_times(this_track).lap)-1;
                stop_idx= 2:2:length(lap_times(this_track).lap)-1;
            else % even number of laps
                start_idx= 1:2:length(lap_times(this_track).lap);
                stop_idx= 2:2:length(lap_times(this_track).lap);
            end
            lap_idx{this_track}= [start_idx' stop_idx'];
            
            % initialise variables
            t_segments_test{this_track}= [];
            test_track_id{this_track}= [];
            test_lap_id{this_track}= [];
            pos1= []; t1= [];
            t_test{this_track}= []; pos_test{this_track}= [];
            
            for this_lap=1:length(start_idx)
                % find RUN periods
                RUN_idx1= lap_times(this_track).lap(start_idx(this_lap)).v_cm > parameters.speed_threshold_laps; % one way
                RUN_idx2= lap_times(this_track).lap(stop_idx(this_lap)).v_cm > parameters.speed_threshold_laps; % return
                % of those find start and stops
                end_gaps1= find(diff(RUN_idx1)== -1);
                start_gaps1= find(diff(RUN_idx1)== 1)+1;
                if ~isempty(end_gaps1) && ~isempty(start_gaps1)
                    if (~isempty(end_gaps1) && isempty(start_gaps1)) || (start_gaps1(1) > end_gaps1(1))
                        start_gaps1= [find(RUN_idx1 == 1,1,'first') start_gaps1];
                    end
                    if start_gaps1(end) > end_gaps1(end)
                        end_gaps1= [end_gaps1 find(RUN_idx1 == 1,1,'last')];
                    end
                    seg1=[]; pos1=[]; t1=[]; k1=1;
                    for this_seg=1:length(start_gaps1)
                        if (start_gaps1(this_seg)+ ceil(parameters.run_bin_width/timestep)) <= end_gaps1(this_seg) % needs to be at least one bin for decoding
                            % keep start:stop
                            seg1= [seg1; start_gaps1(this_seg) end_gaps1(this_seg)]; 
                            % store position and timestamps
                            pos1{k1}= lap_times(this_track).lap(start_idx(this_lap)).x(start_gaps1(this_seg):end_gaps1(this_seg));
                            t1{k1}= lap_times(this_track).lap(start_idx(this_lap)).t(start_gaps1(this_seg):end_gaps1(this_seg));
                            k1= k1+1;
                        end
                    end
                else
                    seg1=[]; pos1=[]; t1=[];
                end
                
                % same for return bit of the lap
                end_gaps2= find(diff(RUN_idx2)== -1);
                start_gaps2= find(diff(RUN_idx2)== 1)+1;
                if ~isempty(end_gaps2) && ~isempty(start_gaps2)
                    if (~isempty(end_gaps2) && isempty(start_gaps2)) || (start_gaps2(1) > end_gaps2(1)) 
                        start_gaps2= [find(RUN_idx2 == 1,1,'first') start_gaps2];
                    end
                    if start_gaps2(end) > end_gaps2(end)
                        end_gaps2= [end_gaps2 find(RUN_idx2 == 1,1,'last')];
                    end
                    seg2=[]; pos2=[]; t2=[]; k2=1;
                    for this_seg=1:length(start_gaps2)
                        if (start_gaps2(this_seg)+ ceil(parameters.run_bin_width/timestep)) <= end_gaps2(this_seg) % needs to be at least one bin for decoding
                            % keep start:stop
                            seg2= [seg2; start_gaps2(this_seg) end_gaps2(this_seg)];
                            % store position and timestamps
                            pos2{k2}= lap_times(this_track).lap(stop_idx(this_lap)).x(start_gaps2(this_seg):end_gaps2(this_seg));
                            t2{k2}= lap_times(this_track).lap(stop_idx(this_lap)).t(start_gaps2(this_seg):end_gaps2(this_seg));
                            k2= k2+1;
                        end
                    end
                else 
                    seg2=[]; pos2=[]; t2=[];
                end
               
                ts1= lap_times(this_track).lap(start_idx(this_lap)).t(seg1);
                ts2= lap_times(this_track).lap(stop_idx(this_lap)).t(seg2);
                t_segments_test{this_track}= [t_segments_test{this_track}; ts1 ; ts2];
                test_track_id{this_track}= [test_track_id{this_track}; this_track*ones(size(ts1,1),1) ;  this_track*ones(size(ts2,1),1)];
                test_lap_id{this_track}= [test_lap_id{this_track}; start_idx(this_lap)*ones(size(ts1,1),1)  ; stop_idx(this_lap)*ones(size(ts2,1),1)];
                t_test{this_track}= [t_test{this_track} t1 t2];
                pos_test{this_track}= [pos_test{this_track} pos1 pos2];

            end
        end
        
        % calculate place_fields for training set (and test set)
        disp('calculating lap place fields..')
        for this_track=1:length(lap_times)
            this_track
            for this_lap=1:length(lap_idx{this_track})
                this_lap
                test_indices= test_lap_id{this_track}== lap_idx{this_track}(this_lap,1) | test_lap_id{this_track}== lap_idx{this_track}(this_lap,2);
                lap_test= cell(1,length(lap_times));
                laps_training= lap_test;
                lap_test{this_track}= t_segments_test{this_track}(test_indices,:);
                laps_training{this_track}= t_segments_test{this_track}(~test_indices,:);
                
                if ~isempty(control)
                    switch control
                        case 'global_remap'
                            place_fields_test= calculate_place_fields_epochs(bin_size,lap_test,'global_remap');
                            place_fields_test.track= place_fields_test.track(this_track);
                            place_fields_training= calculate_place_fields_epochs(bin_size,laps_training,'global_remap');
                            place_fields_training.track= place_fields_training.track(this_track);
                        case 'rate_remap'
                            place_fields_test= calculate_place_fields_epochs(bin_size,lap_test,'rate_remap');
                            place_fields_test.track= place_fields_test.track(this_track);
                            place_fields_training= calculate_place_fields_epochs(bin_size,laps_training,'rate_remap');
                            place_fields_training.track= place_fields_training.track(this_track);
                    end
                else
                    place_fields_test= calculate_place_fields_epochs(bin_size,lap_test);
                    place_fields_test.track= place_fields_test.track(this_track);
                    place_fields_training= calculate_place_fields_epochs(bin_size,laps_training);
                    place_fields_training.track(this_track)= place_fields_training.track(this_track);
                    %
                    remaining_idx= [1:length(place_fields.track)]; 
                    remaining_idx(this_track)= [];
                    fnames= fieldnames(place_fields_training.track);
                    fnames(strcmp(fnames,'field_boundaries'))=[];
                    for ii=1:length(remaining_idx)
                        place_fields.track(remaining_idx(ii)).total_time= diff(place_fields.track(remaining_idx(ii)).time_window);
                        for jj=1:length(fnames)
                        place_fields_training.track(remaining_idx(ii)).(fnames{jj})= place_fields.track(remaining_idx(ii)).(fnames{jj});
                        end
                    end
                end
                
                
                if ~isfield(place_fields_training.track,'x_bin_edges')
                    keyboard;
                end
                 % spike count
                bayesian_spike_count_test= spike_count([],t_segments_test{this_track}(test_indices,1)',t_segments_test{this_track}(test_indices,2)','N','laps');
              
                % decoding
                est_temp = bayesian_decoding(place_fields_training,bayesian_spike_count_test,'N');
                
                % find position on track
                pos= [pos_test{this_track}{test_indices}];
                t= [t_test{this_track}{test_indices}];

                % discretise to the position bins used in the decoding
                discrete_position = NaN(size(pos));
                discrete_position = discretize(pos,place_fields_test.track.x_bin_edges); %group position points in bins delimited by edges
                index = find(~isnan(discrete_position));
                discrete_position(index) = est_temp(this_track).position_bin_centres(discrete_position(index)); %creates new positions based on centre of bins
                interp_discrete_pos= interp1(t,discrete_position, [bayesian_spike_count_test.run_epochs.run_time_centered], 'nearest');
                interp_pos= interp1(t,pos, [bayesian_spike_count_test.run_epochs.run_time_centered], 'nearest');

                % calculate error
                for hh=1:length(place_fields.track)
                estimated_position_test(this_track).run_epochs(this_lap).run_time_edges(hh,:)= [est_temp(hh).run_epochs.run_time_edges];
                estimated_position_test(this_track).run_epochs(this_lap).run_time_centered(hh,:)= [est_temp(hh).run_epochs.run_time_centered];
                estimated_position_test(this_track).run_epochs(this_lap).run{hh,1}= [est_temp(hh).run_epochs.run];
                estimated_position_test(this_track).run_epochs(this_lap).peak_position(hh,:)= [est_temp(hh).run_epochs.peak_position];
                estimated_position_test(this_track).run_epochs(this_lap).run_bias(hh,:)= [est_temp(hh).run_epochs.run_bias];
                estimated_position_test(this_track).run_epochs(this_lap).max_prob(hh,:)= [est_temp(hh).run_epochs.max_prob];
                estimated_position_test(this_track).run_epochs(this_lap).discrete_run_error(hh,:)= abs([est_temp(hh).run_epochs.peak_position]-interp_discrete_pos);
                estimated_position_test(this_track).run_epochs(this_lap).run_error(hh,:)= abs([est_temp(hh).run_epochs.peak_position]-interp_pos);
                estimated_position_test(this_track).run_epochs(this_lap).pos_interp(hh,:)= interp_pos;
                estimated_position_test(this_track).run_epochs(this_lap).discrete_pos_interp(hh,:)= interp_discrete_pos;
                end

            end
            
            estimated_position_test(this_track).position_bin_centres= est_temp.position_bin_centres;
            estimated_position_test(this_track).error_thresh =20; % 20cm - 1 bin
            errors_tmp= [estimated_position_test(this_track).run_epochs.run_error];
            errors_tmp= errors_tmp(this_track,:);
            estimated_position_test(this_track).percent_small_error= numel(find(errors_tmp<20))/numel(errors_tmp);
            estimated_position_test(this_track).median_error= nanmedian(errors_tmp);
        end
        
%        switch bin_size
%            case parameters.x_bins_width_bayesian
%                 save('estimated_position_leave_one_out.mat','estimated_position_test');  
%            case parameters.x_bins_width
%                 save('estimated_position_leave_one_out_SMALL.mat','estimated_position_test');  
%            otherwise
%                 save(['estimated_position_leave_one_out_' num2str(bin_size) 'cm.mat'],'estimated_position_test');  
%        end 
       
     if p.Results.save_figures  
     % create cumulative distr plots
     f1= figure;
     for track_id=1:length(estimated_position_test)
        x_bins= place_fields_test.track.x_bin_centres;
        pd= cellfun(@nanmean, {estimated_position_test(track_id).run_epochs.run_error});
        y_distr= cdf_fraction(pd,x_bins);
        plot(x_bins, y_distr); hold on;
        median_error(track_id)= nanmedian([estimated_position_test(track_id).run_epochs.discrete_run_error]);
        estimated_position_test(track_id).median_error= median_error(track_id);
        total_decoding_time(track_id)= sum(cellfun(@numel, {estimated_position_test(track_id).run_epochs.run_time_centered}))*parameters.run_bin_width;
     end
%      legend({['track1, median: ' num2str(median_error(1))],['track2, median: ' num2str(median_error(2))],['track3, median: ' num2str(median_error(3))]},'Location','southeast')
     xlabel('cm');
     ylabel('cumulative fraction')
     
     % confusion matrices
     f3=figure('Position',[154,284,889,209]);
     for track_id=1:length(estimated_position_test)
        subplot(1,length(estimated_position_test),track_id);
        true_pos= [estimated_position_test(track_id).run_epochs.discrete_pos_interp];
        decoded_pos= [estimated_position_test(track_id).run_epochs.peak_position];
        n= hist3([true_pos' decoded_pos'],'Ctrs',{estimated_position_test(track_id).position_bin_centres estimated_position_test(track_id).position_bin_centres}); 
        sum_prob= sum(n,2);
        imagesc(n./sum_prob);
        set(gca,'YDir','Normal');
        colormap gray
        map= colormap;
        colormap(flipud(map));
        title(['Track ' num2str(track_id)]);
        colorbar
        caxis([0 1]);
        xt= xticks;
        xticklabels(bin_size*xt);
        yt= yticks;
        yticklabels(bin_size*yt);
        xlabel('laps true position (cm)');
        ylabel('estimated position (cm)');
     end

      
        if exist('figures\decoding') ~= 7
            mkdir('figures\decoding')
        end
        switch bin_size
           case parameters.x_bins_width_bayesian
%                 saveas(f1,'figures\decoding\decoding_error_LeaveOneOut_laps.png');
                savefig(f1,'figures\decoding\decoding_error_LeaveOneOut_laps');
%                 saveas(f3,'figures\decoding\confusion_matrices_LeaveOneOut_laps.png');
                savefig(f3,'figures\decoding\confusion_matrices_LeaveOneOut_laps');
            case parameters.x_bins_width
%                 saveas(f1,'figures\decoding\decoding_error_LeaveOneOut_laps_SMALL.png');
                savefig(f1,'figures\decoding\decoding_error_LeaveOneOut_laps_SMALL');
%                 saveas(f3,'figures\decoding\confusion_matrices_LeaveOneOut_laps.png');
                savefig(f3,'figures\decoding\confusion_matrices_LeaveOneOut_laps_SMALL');
            otherwise
                 saveas(f1,['figures\decoding\decoding_error_LeaveOneOut_laps_' num2str(bin_size) 'cm.png']);
%                 savefig(f1,['figures\decoding\decoding_error_LeaveOneOut_laps_' num2str(bin_size) 'cm']);
                saveas(f3,['figures\decoding\decoding_error_LeaveOneOut_laps_' num2str(bin_size) 'cm.png']);
%                 savefig(f3,['figures\decoding\decoding_error_LeaveOneOut_laps_' num2str(bin_size) 'cm']);
        end
    end

    case 'cross_tracks'
        disp('running cross tracks')
        switch bin_size
            case parameters.x_bins_width_bayesian
               load('extracted_place_fields_BAYESIAN.mat');
               place_fields_training= place_fields_BAYESIAN;
            case parameters.x_bins_width
               load('extracted_place_fields.mat');
               place_fields_training= place_fields;
            otherwise
                disp('cannot run with this size bin at the moment')
                keyboard;
       end
       load('extracted_position.mat');
       timestep= mean(diff(position.t));
       
       for this_track=1:length(place_fields_training.track)
            % keep only run epochs
                RUN_idx= (position.t >= min(position.t(~isnan(position.linear(this_track).linear))) &...
                         position.t <= max(position.t(~isnan(position.linear(this_track).linear))) &...
                         position.v_cm > parameters.speed_threshold_laps);
                                
                end_gaps= find(diff(RUN_idx)== -1);
                start_gaps= find(diff(RUN_idx)== 1)+1;
                if ~isempty(end_gaps) && ~isempty(start_gaps)
                    if (~isempty(end_gaps) && isempty(start_gaps)) || (start_gaps(1) > end_gaps(1))
                        start_gaps= [find(RUN_idx == 1,1,'first') start_gaps];
                    end
                    if start_gaps(end) > end_gaps(end)
                        end_gaps= [end_gaps find(RUN_idx == 1,1,'last')];
                    end
                    seg=[]; pos=[]; t=[]; k=1;
                    for this_seg=1:length(start_gaps)
                        if (start_gaps(this_seg)+ ceil(parameters.run_bin_width/timestep)) <= end_gaps(this_seg) % needs to be at least one bin for decoding
                            % keep start:stop
                            seg= [seg; start_gaps(this_seg) end_gaps(this_seg)]; 
                            % store position and timestamps
                            pos= [pos position.linear(this_track).linear(start_gaps(this_seg):end_gaps(this_seg))];
                            t= [t position.t(start_gaps(this_seg):end_gaps(this_seg))];
                            k= k+1;
                        end
                    end
                else
                    seg=[]; pos=[]; t=[];
                end
                
                ts_seg_training= position.t(seg);
                % spike count
                switch bin_size
                    case parameters.x_bins_width_bayesian
                        spike_count_track= spike_count([],ts_seg_training(:,1)',ts_seg_training(:,2)','N','laps');
                    case parameters.x_bins_width
                        load('extracted_place_fields.mat');
                        spike_count_track= spike_count(place_fields,ts_seg_training(:,1)',ts_seg_training(:,2)','N','laps');
                    otherwise
                        disp('running with bayesian place fields')
                        spike_count_track= spike_count([],ts_seg_training(:,1)',ts_seg_training(:,2)','N','laps');
                end
                % use place fields from other tracks to decode spikes of
                % current track
                est_temp = bayesian_decoding(place_fields_training,spike_count_track,'N');
                
                % discretise to the position bins used in the decoding
                discrete_position = NaN(size(pos));
                discrete_position = discretize(pos,place_fields_training.track(this_track).x_bin_edges); %group position points in bins delimited by edges
                index = find(~isnan(discrete_position));
                discrete_position(index) = est_temp(this_track).position_bin_centres(discrete_position(index)); %creates new positions based on centre of bins
                interp_pos= interp1(t,discrete_position, [spike_count_track.run_epochs.run_time_centered], 'nearest');
                
                % calculate error
                for this_other_track=1:length(est_temp)
                    estimated_position_test(this_track).run_epochs(this_other_track).run_time_edges= [est_temp(this_other_track).run_epochs.run_time_edges];
                    estimated_position_test(this_track).run_epochs(this_other_track).run_time_centered= [est_temp(this_other_track).run_epochs.run_time_centered];
                    estimated_position_test(this_track).run_epochs(this_other_track).run= [est_temp(this_other_track).run_epochs.run];
                    estimated_position_test(this_track).run_epochs(this_other_track).peak_position= [est_temp(this_other_track).run_epochs.peak_position];
                    estimated_position_test(this_track).run_epochs(this_other_track).run_bias= [est_temp(this_other_track).run_epochs.run_bias];
                    estimated_position_test(this_track).run_epochs(this_other_track).max_prob= [est_temp(this_other_track).run_epochs.max_prob];
                    estimated_position_test(this_track).run_epochs(this_other_track).run_error= abs([est_temp(this_other_track).run_epochs.peak_position]-interp_pos);
                    estimated_position_test(this_track).run_epochs(this_other_track).true_pos_interp= interp_pos;
                    estimated_position_test(this_track).run_epochs(this_other_track).error_thresh =20; % 20cm - 1 bin
                    errors_tmp= [estimated_position_test(this_track).run_epochs(this_other_track).run_error];
                    estimated_position_test(this_track).run_epochs(this_other_track).percent_small_error= numel(find(errors_tmp<20))/numel(errors_tmp);
                end
                % get classification accuracy
        all_prob= vertcat(estimated_position_test(this_track).run_epochs.run); % look at all probabilities
        [max_prob,max_prob_pos]= max(all_prob,[],1); % get index, first track will be first rows etc.. 
        no_spk_bins= find(sum(all_prob,1)==0); % if there are no spikes, get a column of zeros, but find gives you an index of 1
        max_prob_pos(no_spk_bins)=[];
        max_prob(no_spk_bins)=[];

        num_bins_tracks= cellfun(@(x) size(x,1), {estimated_position_test(this_track).run_epochs.run});
        delim_bins_tracks(:,1)= [1 cumsum(num_bins_tracks(1:end-1))+1];
        delim_bins_tracks(:,2)= cumsum(num_bins_tracks);

        track_classification = NaN(size(max_prob_pos));
        for that_track=1:length(estimated_position_test(this_track).run_epochs)
            % classify which track max prob falls in
            indices= max_prob_pos >= delim_bins_tracks(that_track,1) & max_prob_pos <= delim_bins_tracks(that_track,2);
            track_classification(indices)= that_track;
        end
   
        for that_track=1:length(estimated_position_test(this_track).run_epochs)
            track_classification_tmp= track_classification;       
            t_bins= estimated_position_test(this_track).run_epochs(that_track).run_time_centered;
            t_bins(no_spk_bins)= [];
            % to be fair, only take times where v_cm > 5cm/s
            speed_discrete= interp1(position.t, position.v_cm, t_bins, 'nearest');
            pos_discrete=estimated_position_test(this_track).run_epochs(that_track).true_pos_interp;
            pos_discrete(no_spk_bins)=[];
            low_speed_idx= speed_discrete < parameters.speed_threshold_laps;
            track_classification_tmp(low_speed_idx)= NaN;
            pos_discrete(low_speed_idx)= NaN;

            good_classification= numel(find(track_classification_tmp(~isnan(pos_discrete)) == that_track));
            all_classification= numel(find(~isnan(track_classification_tmp(~isnan(pos_discrete)))));
            estimated_position_test(this_track).run_epochs(that_track).classification_accuracy = good_classification/all_classification;
        end          
        
       end
      
       % save output
%         switch bin_size
%             case parameters.x_bins_width_bayesian
%                 save('estimated_position_cross_tracks.mat','estimated_position_test');
%             case parameters.x_bins_width
%                save('estimated_position_cross_tracks_SMALL.mat','estimated_position_test');
%             otherwise
%                save(['estimated_position_cross_tracks_' num2str(bin_size) 'cm.mat'],'estimated_position_test');
%         end
        
        
        if p.Results.save_figures
            f1= figure;
            for this_track=1:length(estimated_position_test)
                subplot(1,length(estimated_position_test),this_track);
                run_error= {estimated_position_test(this_track).run_epochs.run_error};
                track_group= cellfun(@(x,y) y*ones(size(x)) ,run_error,num2cell(1:length(run_error)),'UniformOutput',false);
                run_error= [run_error{:}];
                track_group= [track_group{:}];
                boxplot(run_error,track_group,'PlotStyle','compact','ColorGroup',track_group,'LabelOrientation','horizontal','Symbol','');
                ylabel('decoding error while running');
                xlabel('track used to decode');
                title(['decoded track: ' num2str(this_track)]);
            end

            if exist('figures\decoding') ~= 7
                mkdir('figures\decoding')
            end

            switch bin_size
                case parameters.x_bins_width_bayesian
                    saveas(f1,'figures\decoding\decoding_error_cross_tracks.png');
                    savefig(f1,'figures\decoding\decoding_error_cross_tracks');
                case parameters.x_bins_width_bayesian
                    saveas(f1,'figures\decoding\decoding_error_cross_tracks_SMALL.png');
                    savefig(f1,'figures\decoding\decoding_error_cross_tracks_SMALL');
                otherwise
                    saveas(f1,['figures\decoding\decoding_error_cross_tracks_' num2str(bin_size) 'cm.png']);
                    savefig(f1,['figures\decoding\decoding_error_cross_tracks_' num2str(bin_size) 'cm']);
            end
        end

end
end
