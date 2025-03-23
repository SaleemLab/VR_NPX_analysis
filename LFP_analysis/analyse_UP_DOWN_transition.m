function analyse_UP_DOWN_transition(slow_waves_all,event_info,ripples_all,UP_DOWN_ripple_PSTH_MUA,sessions_to_process)


for nprobe = 1:2
    timeToTransition = event_info(nprobe).UP_duration;
    % event_info(1).DOWN_duration

    % Define survival response variables
    timeToTransition = UP_duration; % Time until UP to DOWN transition
    event = Transition_UP_DOWN; % 1 if transitioned, 0 if censored

    % Define predictor variables
    X = [Ripple_count, Ripple_rate, Hipp_spike_rate, V1_spike_rate];

    % Fit Cox proportional hazards model
    coxMdl = coxphfit(X, timeToTransition, 'Censoring', 1-event);

    % Display results
    disp('Cox Model Coefficients:');
    disp(coxMdl);

    % Plot survival function (Kaplan-Meier estimate)
    figure;
    ecdf(timeToTransition, 'Function', 'survivor');
    xlabel('Time (s)');
    ylabel('Survival Probability');
    title('Kaplan-Meier Survival Curve');

    % Hazard ratio interpretation
    HR = exp(coxMdl);
    disp('Hazard Ratios:');
    disp(HR);
end
