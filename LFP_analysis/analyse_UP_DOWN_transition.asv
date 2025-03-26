function analyse_UP_DOWN_transition(event_info,ripples_all,spindles_all,slow_waves_all)

% if exist('D:\corticohippocampal_replay')>0
%     analysis_folder = 'D:\corticohippocampal_replay';
% elseif exist('P:\corticohippocampal_replay')>0
%     analysis_folder = 'P:\corticohippocampal_replay';
% end
% load(fullfile(analysis_folder,'slow_waves_all_POST.mat'))

hemisphere = {'L','R'};
normalised_duration = 1/40:1/20:(1-1/40);

unique_duration = unique(event_info(1).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{1}))(:,3));
if length(unique_duration) == length(normalised_duration)
    for nprobe = 1:2
        for mprobe = 1:2
            unique_duration = unique(event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(:,3));
            for i = 1:length(unique_duration)
                event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))...
                    (event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(:,3) == unique_duration(i),3) = normalised_duration(i);

                event_info(nprobe).(sprintf('%s_ripple_normalised_DOWN_duration',hemisphere{mprobe}))...
                    (event_info(nprobe).(sprintf('%s_ripple_normalised_DOWN_duration',hemisphere{mprobe}))(:,3) == unique_duration(i),3) = normalised_duration(i);
            end
        end
    end
end

duration_threshold = median(unique(event_info(1).L_ripple_normalised_UP_duration(:,3)));


for nprobe = 1:2
    for mprobe = 1:2

        % Define survival response variables
        timeToTransition = event_info(nprobe).UP_duration; % Time until UP to DOWN transition

        L_V1_MUA_activity{nprobe}{nprobe} = [];
        L_HPC_MUA_activity{nprobe}{nprobe} = [];
        R_V1_MUA_activity{nprobe}{nprobe} = [];
        R_HPC_MUA_activity{nprobe}{nprobe} = [];

        % find predictor variables for each UP event
        normalised_ripple_duration{nprobe}{mprobe} = [];
        ripple_duration{nprobe}{mprobe} = [];
        ripple_count{nprobe}{mprobe} = [];
        ripple_power{nprobe}{mprobe} = [];
        ripple_rate{nprobe}{mprobe} = [];

        early_ripple_count{nprobe}{mprobe} = [];
        late_ripple_count{nprobe}{mprobe} = [];

        ripple_L_V1_MUA{nprobe}{mprobe} = [];
        ripple_L_HPC_MUA{nprobe}{mprobe}= [];
        ripple_R_V1_MUA{nprobe}{mprobe} = [];
        ripple_R_HPC_MUA{nprobe}{mprobe} = [];

        time_to_first_ripples{nprobe}{mprobe} = [];
        time_from_last_ripples{nprobe}{mprobe} = [];
        time_from_mean_ripples{nprobe}{mprobe}= [];

        if mprobe ==1
            unique_UP_index = unique(event_info(nprobe).L_ripple_normalised_UP_duration(:,2));
        else
            unique_UP_index = unique(event_info(nprobe).R_ripple_normalised_UP_duration(:,2));
        end

        for nevent = 1:length(event_info(nprobe).UP_duration)
            index = find(event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(:,2)==event_info(nprobe).UP_index(nevent));

            UP_ints = slow_waves_all(nprobe).UP_ints(event_info(nprobe).UP_index(nevent),:);

            if isempty(index)
                early_ripple_count{nprobe}{mprobe}(nevent) = 0;
                late_ripple_count{nprobe}{mprobe}(nevent) = 0;
                normalised_ripple_duration{nprobe}{mprobe}(nevent) = 0;
                ripple_duration{nprobe}{mprobe}(nevent) = 0;
                ripple_count{nprobe}{mprobe}(nevent) = 0;
                ripple_power{nprobe}{mprobe}(nevent) = 0;
                ripple_rate{nprobe}{mprobe}(nevent) = 0;
                ripple_L_V1_MUA{nprobe}{mprobe}(nevent) = 0;
                ripple_L_HPC_MUA{nprobe}{mprobe}(nevent) = 0;
                ripple_R_V1_MUA{nprobe}{mprobe}(nevent) = 0;
                ripple_R_HPC_MUA{nprobe}{mprobe}(nevent) = 0;

                time_to_first_ripples{nprobe}{mprobe}(nevent) = 0;
                time_from_last_ripples{nprobe}{mprobe}(nevent) = 0;
                time_from_mean_ripples{nprobe}{mprobe}(nevent) = 0;


            else

                ripples_index = event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(index,1);
                
                time_to_first_ripples{nprobe}{mprobe}(nevent) = min(ripples_all(mprobe).peaktimes(ripples_index))-UP_ints(1);
                time_from_last_ripples{nprobe}{mprobe}(nevent) = UP_ints(2)-max(ripples_all(mprobe).peaktimes(ripples_index));
                time_from_mean_ripples{nprobe}{mprobe}(nevent) = UP_ints(2)-mean(ripples_all(mprobe).peaktimes(ripples_index));

                
                late_ripple_count{nprobe}{mprobe}(nevent) = sum(event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(index,3) > 0.5);
                early_ripple_count{nprobe}{mprobe}(nevent) = sum(event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(index,3) <= 0.5);

                ripple_duration{nprobe}{mprobe}(nevent) = event_info(nprobe).UP_duration(nevent)*...
                    event_info(nprobe).(sprintf('%s_ripple_cumulative_duration_UP',hemisphere{mprobe}))(unique_UP_index == unique(event_info(nprobe).(sprintf('%s_ripple_normalised_UP_duration',hemisphere{mprobe}))(index,2)));

                if length(index)>1
                    normalised_ripple_duration{nprobe}{mprobe}(nevent) = ripple_duration{nprobe}{mprobe}(nevent)...
                        /(max(ripples_all(mprobe).peaktimes(ripples_index))-min(ripples_all(mprobe).peaktimes(ripples_index))); % normalised by time windows where first and last ripple happens

                    ripple_rate{nprobe}{mprobe}(nevent) = length(index)...
                        /(max(ripples_all(mprobe).peaktimes(ripples_index))-min(ripples_all(mprobe).peaktimes(ripples_index))); % normalised by time windows where first and last ripple happens
                else
                    normalised_ripple_duration{nprobe}{mprobe}(nevent) = ripple_duration{nprobe}{mprobe}(nevent);

                    ripple_rate{nprobe}{mprobe}(nevent) = length(index)...
                        /(ripples_all(mprobe).offset(ripples_index)-ripples_all(mprobe).onset(ripples_index)); % normalised by time windows where first and last ripple happens
                end
                ripple_count{nprobe}{mprobe}(nevent) = length(index);

                [ripple_power{nprobe}{mprobe}(nevent),temp] = max(event_info(nprobe).(sprintf('%s_ripple_zscore_UP',hemisphere{mprobe}))(index));

                % peri-ripple MUA peak zscore
                ripple_L_HPC_MUA{nprobe}{mprobe}(nevent) = event_info(nprobe).(sprintf('%s_ripple_HPC_MUA_peak_UP',hemisphere{mprobe}))(index(temp),1);
                ripple_R_HPC_MUA{nprobe}{mprobe}(nevent) = event_info(nprobe).(sprintf('%s_ripple_HPC_MUA_peak_UP',hemisphere{mprobe}))(index(temp),2);
                ripple_L_V1_MUA{nprobe}{mprobe}(nevent) = event_info(nprobe).(sprintf('%s_ripple_V1_MUA_peak_UP',hemisphere{mprobe}))(index(temp),1);
                ripple_R_V1_MUA{nprobe}{mprobe}(nevent) = event_info(nprobe).(sprintf('%s_ripple_V1_MUA_peak_UP',hemisphere{mprobe}))(index(temp),2);
            end

        end
    end


    L_V1_MUA_DU_slope{nprobe} = event_info(nprobe).DU_slope_V1(1,:);
    % L_HPC_MUA_activity{nprobe} = event_info(nprobe).L_HPC_cumulative_activity_UP;
    R_V1_DU_slope{nprobe} = event_info(nprobe).DU_slope_V1(2,:);
    % R_HPC_MUA_activity{nprobe} = event_info(nprobe).R_HPC_cumulative_activity_UP;

    % Cumulative MUA activity (0-99% normalised so that cumulative activity is addictive)
    L_V1_MUA_activity{nprobe} = event_info(nprobe).L_V1_cumulative_activity_UP;
    L_HPC_MUA_activity{nprobe} = event_info(nprobe).L_HPC_cumulative_activity_UP;
    R_V1_MUA_activity{nprobe} = event_info(nprobe).R_V1_cumulative_activity_UP;
    R_HPC_MUA_activity{nprobe} = event_info(nprobe).R_HPC_cumulative_activity_UP;

    % Define predictor variables
    timeToTransition_ripples = time_from_last_ripples{1}{1}(normalised_ripple_duration{1}{1}>0);
    X = [ripple_L_HPC_MUA{1}(normalised_ripple_duration{1}>0); ripple_power{1}(normalised_ripple_duration{1}>0)]';
    % X = [normalised_ripple_duration{1}; ripple_power{1}; ripple_L_V1_MUA{1}; ripple_R_V1_MUA{1};...
    %     L_V1_MUA_activity{1}; R_V1_MUA_activity{1}]';
    % 
    % scatter(normalised_ripple_duration{1}(normalised_ripple_duration{1}>0),timeToTransition(normalised_ripple_duration{1}>0),'filled','MarkerFaceAlpha','0.02')
    % scatter(ripple_power{1}(normalised_ripple_duration{1}>0),timeToTransition(normalised_ripple_duration{1}>0),'filled','MarkerFaceAlpha','0.02')

    % Fit Cox proportional hazards model
    [b, logL, H, stats] = coxphfit(X, timeToTransition_ripples);

    % Display results
    disp('Cox Model Coefficients:');
    disp(coxMdl);



    % Plot survival function (Kaplan-Meier estimate)

    % ripple rate -> time from last ripple to UP end
    figure;
    timeToTransition_ripples = time_from_last_ripples{nprobe}{mprobe};
    ripple_threshold1 = prctile(ripple_rate{nprobe}{mprobe}(ripple_count{nprobe}{mprobe}>1),25);
    ripple_threshold2 = prctile(ripple_rate{nprobe}{mprobe}(ripple_count{nprobe}{mprobe}>1),75);

    ecdf(timeToTransition_ripples(ripple_rate{nprobe}{mprobe}>=ripple_threshold2), 'Function', 'survivor');
    hold on;
    ecdf(timeToTransition_ripples(ripple_rate{nprobe}{mprobe}<=ripple_threshold1&ripple_count{nprobe}{mprobe}>1), 'Function', 'survivor');
    legend('top 25% ripple rate','bottom 25% ripple rate','box off')

    xlabel('Time (s)');
    ylabel('Survival Probability of UP');
    title('Kaplan-Meier Survival Curve');

    
    % ripple count -> time from last ripple to UP end (confound....)
    figure;
    timeToTransition_ripples = time_from_last_ripples{nprobe}{mprobe};
    ripple_threshold1 = prctile(ripple_count{nprobe}{mprobe}(ripple_count{nprobe}{mprobe}>0),25);
    ripple_threshold2 = prctile(ripple_count{nprobe}{mprobe}(ripple_count{nprobe}{mprobe}>0),75);

    ecdf(timeToTransition_ripples(ripple_count{nprobe}{mprobe}>=ripple_threshold2), 'Function', 'survivor');
    hold on;
    ecdf(timeToTransition_ripples(ripple_count{nprobe}{mprobe}<=ripple_threshold1&ripple_count{nprobe}{mprobe}>0), 'Function', 'survivor');
    legend('2 ripples or more','1 ripple')

    xlabel('Time (s)');
    ylabel('Survival Probability of UP');
    title('Kaplan-Meier Survival Curve');


    % ripple power
    figure;
    timeToTransition_ripples = time_from_last_ripples{nprobe}{mprobe};
    ripple_threshold1 = prctile(ripple_power{nprobe}{mprobe}(ripple_count{nprobe}{mprobe}>0),25);
    ripple_threshold2 = prctile(ripple_power{nprobe}{mprobe}(ripple_count{nprobe}{mprobe}>0),75);

    ecdf(timeToTransition_ripples(ripple_power{nprobe}{mprobe}>=ripple_threshold2), 'Function', 'survivor');
    hold on;
    ecdf(timeToTransition_ripples(ripple_power{nprobe}{mprobe}<=ripple_threshold1&ripple_count{nprobe}{mprobe}>0), 'Function', 'survivor');
    legend('top 25% ripple power','bottom 25% ripple power','box off')

    xlabel('Time (s)');
    ylabel('Survival Probability of UP');
    title('Kaplan-Meier Survival Curve');

    







    % Plot survival function (Kaplan-Meier estimate)
    figure;
    ripple_threshold1 = prctile(ripple_R_V1_MUA{1}(normalised_ripple_duration{1}>0),25);
    ripple_threshold2 = prctile(ripple_R_V1_MUA{1}(normalised_ripple_duration{1}>0),75);

    ecdf(timeToTransition(ripple_L_V1_MUA{1}>ripple_threshold2), 'Function', 'survivor');
    hold on;
    ecdf(timeToTransition(ripple_L_V1_MUA{1}<ripple_threshold1&normalised_ripple_duration{1}>0), 'Function', 'survivor');
    legend('top 25% normalised V1 MUA during ripples','bottom 25% normalised V1 MUA during ripples')

    xlabel('Time (s)');
    ylabel('Survival Probability of UP');
    title('Kaplan-Meier Survival Curve');

    % Plot survival function (Kaplan-Meier estimate)
    figure;
    ripple_threshold1 = prctile(R_V1_MUA_activity{1},25);
    ripple_threshold2 = prctile(R_V1_MUA_activity{1},75);

    ecdf(timeToTransition(R_V1_MUA_activity{1}>ripple_threshold2), 'Function', 'survivor');
    hold on;
    ecdf(timeToTransition(R_V1_MUA_activity{1}<ripple_threshold1), 'Function', 'survivor');
    legend('top 25% normalised HPC MUA during ripples','bottom 25% normalised HPC MUA during ripples')

    xlabel('Time (s)');
    ylabel('Survival Probability of UP');
    title('Kaplan-Meier Survival Curve');

    % Hazard ratio interpretation
    HR = exp(b);
    disp('Hazard Ratios:');
    disp(HR);

end
