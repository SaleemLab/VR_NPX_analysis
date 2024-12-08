function slow_waves_markov = detect_UP_DOWN_markov(tvec,slow_waves,spiketimes,behavioural_state)


% % Filter the LFP: delta, high gamma and get MUA spike counts
% filterparms.deltafilter = [0.5 8];%heuristically defined.  room for improvement here.
% filterparms.gammafilter = [100 400];
% filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
% filterparms.gammanormwin = 20; %window for gamma normalization (s)
% 
% SR = round(1/mean(diff(tvec)));
% display('Filtering LFP')
% best_channel = find(LFP.best_V1_channel==slow_waves.best_channel);
% lfp = [];
% lfp.data = LFP.best_V1(best_channel,:)';
% 
% if size(tvec,1)<size(tvec,2)
%     lfp.timestamps = tvec';
% else
%     lfp.timestamps = tvec;
% end
% 
% lfp.samplingRate = SR;
% 
% %%% Delta power
% deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
% deltaLFP.normamp = NormToInt(deltaLFP.data,'modZ',NREMInts,SR);
% 
% %%% Gamma power
% gammaLFP= bz_Filter(lfp,'passband',filterparms.gammafilter,'filter','fir1','order',4);
% gammaLFP.smoothamp = smooth(gammaLFP.amp,round(filterparms.gammasmoothwin.*SR),'moving');
% gammaLFP.normamp = NormToInt(gammaLFP.smoothamp,'modZ',NREMInts,SR,'moving',filterparms.gammanormwin);
% 
% %%% Spikecount smoothed with σ = 30 ms
% if ~isempty(NREMInts)
%     sleep_tvec = Restrict(tvec,NREMInts);
%     sleep_index = ismember(tvec,sleep_tvec);
% end

%%%%%%%%%%%%%%%%% Hidden Markov UP DOWN ;
% Parameters for the detection algorithm
% (based on https://elifesciences.org/articles/22425#s4)
NREMInts = behavioural_state.SWS;

% Define observation sequences (spike rate, delta power, gamma power)
timebin_size = 0.01;
tvec_interp1 = tvec(1):timebin_size:tvec(end);
tvec_edges = [tvec_interp1(1)-1/(1/mean(diff(tvec_interp1))*2) tvec_interp1+1/(1/mean(diff(tvec_interp1))*2)];
w = gausswin(0.03*1/mean(diff(tvec_interp1))); % Smoothed with σ = 30 ms
w = w / sum(w);
% spikeCounts =histcounts(spiketimes,tvec_edges);
spikeCounts = round(filtfilt(w,1,histcounts(spiketimes(:,1),tvec_edges)')');
spikeCounts(spikeCounts>50*(timebin_size/0.01))=50*(timebin_size/0.01);% limit spike count due to noise
% observations = [spikeCounts',interp1(tvec,deltaLFP.normamp,tvec_interp1)',interp1(tvec,gammaLFP.normamp,tvec_interp1)'];
% spikeCounts = observations(:,1); %



%%%%%%% Calculate transition probability based on the output of Buzaki algorithm
sleep_tvec = Restrict(tvec_interp1,NREMInts);

% Initialize state series (0 = DOWN, 1 = UP)
state_series = zeros(size(sleep_tvec));
% Based on Ji and Wilson, the V1 frame should be peak around 0.5s-1s and majority of them will be less than 2-3s
% For quantification of UP DOWN stats for Markov model, only UP
% with less than 2 seconds will be taken as events longer than 2s may be due to missing DOWN states.
UP_ints = slow_waves.ints.UP(diff(slow_waves.ints.UP,[],2)<2,:);
DOWN_ints = slow_waves.ints.DOWN;

% Assign UP states based on up_onsets and up_offsets
for i = 1:size(UP_ints,1)
    % Find indices in tvec corresponding to the UP state
    up_indices = sleep_tvec >= UP_ints(i,1) & sleep_tvec <= UP_ints(i,2);
    state_series(up_indices) = 1; % Set UP state to 1
end
% Compute transitions
transitions = diff(state_series);  % Differences between consecutive tvec elements

% Count transitions
N_DD = sum(transitions == 0 & state_series(1:end-1) == 0); % DOWN -> DOWN
N_DU = sum(transitions == 1);                              % DOWN -> UP
N_UU = sum(transitions == 0 & state_series(1:end-1) == 1); % UP -> UP
N_UD = sum(transitions == -1);                             % UP -> DOWN

% Total states
N_D = N_DD + N_DU;  % Total DOWN states
N_U = N_UU + N_UD;  % Total UP states

% Transition probabilities
P_DD = N_DD / N_D;        % Transition probability from DOWN to UP
P_DU = N_DU / N_D;        % Transition probability from UP to DOWN
P_UU = N_UU / N_U;        % Self-transition in DOWN state
P_UD = N_UD / N_U;        % Self-transition in UP state


% Initialize parameters for EM
J = round(0.02/timebin_size);                   % Number of history bins (corresponds to 20 ms memory)


% mu(1) = log(mean(observations(status,1))); % Average FR

% Initialize emission parameters alpha (UP firing rates), mu (UP-DOWN
% Firing rate difference) and beta (History-dependence weight)
[status,~,index] = InIntervals(tvec_interp1,[UP_ints(:,1) UP_ints(:,1)+0.1]); % take average firing rate
alpha = log(mean(spikeCounts(status))); % UP Firing rate

if alpha>log(50*(timebin_size/0.01)) | alpha<log(3*(timebin_size/0.01))
    alpha = log(20*(timebin_size/0.01));
end

[status,~,index] = InIntervals(tvec_interp1,[DOWN_ints(:,1) DOWN_ints(:,1)+0.1]);
mu = log(mean(spikeCounts(status)))-alpha; % UP-DOWN firing rate difference
% mu = log(mean(spikeCounts(status)))-alpha; % UP-DOWN firing rate difference
if mu >-1 & mu <-5
    mu = -2;
end
% mu = -4.7
beta = 0.01*(0.01/timebin_size);               % History-dependence weight

% alpha = [0, log(3)];      % State-specific rates (DOWN=0, UP=3)
[status,interval,index] = InIntervals(tvec_interp1,NREMInts);
numBins = sum(status);           % Number of time bins during SLEEP
maxIter = 100;            % Maximum EM iterations

% P = [P_DD,P_DU;P_UD,P_UU]; % Transition probabilities
P = [0.9,0.1;0.1,0.9]; % Transition probabilities
Pi = [0.7 0.3];% rough UP and DOWN probabilities
% Pi= [mean(diff(DOWN_ints,[],2)) mean(diff(UP_ints,[],2))]/(mean(diff(DOWN_ints,[],2))+mean(diff(UP_ints,[],2)));


stateProbs = zeros(numBins, 2); % Forward-backward probabilities
epoch_id = interval(interval~=0);

% Initialize parameters
gamma_t = zeros(numBins, 2); % Posterior probabilities of each state (UP/DOWN)
xi_t = zeros(numBins-1, 2, 2); % Joint state probabilities for transitions

% EM Algorithm Setup
maxIter = 100; tol = 1e-5;  % Convergence criteria
logLikelihoods = zeros(maxIter, 1);

% EM Algorithm
for iter = 1:maxIter
    tic
    alpha_t = zeros(numBins, 2); % [P(y1:t, St=DOWN), P(y1:t, St=UP)]
    beta_t = zeros(numBins, 2);
    xi_t = zeros(numBins-1, 2, 2);
    emission_t = [];
    for nEpoch = 1:size(NREMInts,1)
        tidx = find(interval==nEpoch); % tidx out of the entire session
        NREM_tidx = find(epoch_id == nEpoch); % tidx during sleep windows

        %%%%% E-Step: Forward-Backward Algorithm
        % Forward probabilities
        for t = 1:length(tidx)
            if t == 1
                alpha_t(NREM_tidx(t), :) = Pi .* emission_prob(spikeCounts(tidx(1)), mu, alpha, beta, [], spikeCounts); % Initial probabilities
            else
                alpha_t(NREM_tidx(t), :) = (alpha_t(NREM_tidx(t)-1, :) * P) .* emission_prob(spikeCounts(tidx(t)), mu, alpha, beta, tidx(t)-J:tidx(t)-1, spikeCounts);
            end
            alpha_t(NREM_tidx(t), :) = alpha_t(NREM_tidx(t), :) / sum(alpha_t(NREM_tidx(t), :)); % Normalize

            emission_t(NREM_tidx(t),:) = emission_prob(spikeCounts(tidx(t)), mu, alpha, beta, tidx(t)-J:tidx(t)-1, spikeCounts);
        end

        % Backward probabilities

        beta_t(NREM_tidx(end), :) = 1; % Initialization
        for t = length(tidx)-1:-1:1
            beta_t(NREM_tidx(t), :) = (beta_t(NREM_tidx(t)+1, :) .* emission_prob(spikeCounts(tidx(t)+1), mu, alpha, beta, tidx(t):tidx(t)+J, spikeCounts)) * P';
            beta_t(NREM_tidx(t), :) = beta_t(NREM_tidx(t), :) / sum(beta_t(NREM_tidx(t), :)); % Normalize
        end


        for t = 1:length(tidx)-1
            xi_t(NREM_tidx(t), :, :) = (alpha_t(NREM_tidx(t), :)' * (beta_t(NREM_tidx(t)+1, :) .* emission_prob(spikeCounts(tidx(t)+1), mu, alpha, beta, tidx(t):tidx(t)+J, spikeCounts))) .* P;
            xi_t(NREM_tidx(t), :, :) = xi_t(NREM_tidx(t), :, :) ./ sum(sum(xi_t(NREM_tidx(t), :, :)));
        end
    end

    % State probabilities
    gamma_t = alpha_t .* beta_t; % Joint probability
    gamma_t = gamma_t ./ sum(gamma_t, 2);

    % === M-Step: Parameter Re-estimation ===
    % Update transition probabilities
    for i = 1:2
        for j = 1:2
            P(i, j) = sum(xi_t(:, i, j)) / sum(gamma_t(1:end-1, i));
        end
    end

    % Update emission parameters using Newton-Raphson
    for param = {'alpha', 'mu'}
        param = param{1};
        param_value = eval(param);
        f = @(x) compute_log_likelihood(x, param, gamma_t, spikeCounts, J, mu, alpha, beta, NREMInts, interval,epoch_id);
        grad = @(x) compute_first_derivative(f, x);
        hess = @(x) compute_second_derivative(f, x);

        for n = 1:10
            gradient = grad(param_value);
            hessian = hess(param_value);

            % Update parameter
            param_value = param_value - gradient / hessian;

            % Check convergence
            if abs(gradient) < tol
                break;
            end
        end
        % Update parameter value
        eval([param ' = param_value;']);
    end

    logLikelihoods(iter) = compute_log_likelihood(alpha, 'alpha', gamma_t, spikeCounts, J, mu, alpha, beta, NREMInts, interval,epoch_id);
    alpha_iter(iter) =  alpha;
    mu_iter(iter) = mu;
    sprintf('Iter %i Loglikelihood = %i, alpha = %f, mu = %f',iter, logLikelihoods(iter) ,alpha,mu)
    toc
    % Check for convergence (if the change in log-likelihood is smaller than a threshold)
    if iter > 1 && abs(logLikelihoods(iter) - logLikelihoods(iter-1)) < 1e-5
        disp('Log likelihoods converged');
        break;
    end

    if iter > 1 && (abs(alpha_iter(iter) - alpha_iter(iter-1)) < 1e-5 | abs(mu_iter(iter) - mu_iter(iter-1)) < 1e-5)
        disp('Converged');
        break;
    end
end

% === Viterbi Algorithm: Decode Most Likely State Sequence ===
viterbiStates = zeros(numBins, 1);
[~, viterbiStates(1)] = max(gamma_t(1, :));
for t = 2:numBins
    [~, viterbiStates(t)] = max(gamma_t(t, :));
end
viterbiStates=viterbiStates-1;

% % Plot Results
% figure;
% plot( tvec_interp1(status),gamma_t.*15);
% xlabel('Time Bins');
% hold on;
% plot( tvec_interp1(status),spikeCounts(status), 'k');
% xlabel('Time Bins');
%
% subplot(2, 1, 2);
% stairs(viterbiStates, 'r', 'LineWidth', 2);
% xlabel('Time Bins'); ylabel('Inferred State'); title('Inferred State Sequence');
% ylim([-0.5, 1.5]); yticks([0, 1]); yticklabels({'DOWN', 'UP'});

viterbiStatesCleaned = viterbiStates; % Copy for cleaning

%%%%% Step 1 Remove noisy, rapid transitions
% Identify rapid transitions (states lasting < 3 bins)
rapidTransitionIndices = find(diff(viterbiStates) ~= 0);
rapidTransitionDurations = diff([0; rapidTransitionIndices; numBins]);

% Initialize variables to track noisy regions
isNoise = false(numBins, 1);
rapidAlternations = 0;

% Iterate through transitions to detect noisy regions
for i = 1:length(rapidTransitionIndices)
    if rapidTransitionDurations(i) <= 3
        rapidAlternations = rapidAlternations + 1;
    else
        rapidAlternations = 0; % Reset if a longer duration is found
    end

    % Mark as noise if more than 3 rapid alternations occur
    if rapidAlternations >= 3
        startIdx = max(1, rapidTransitionIndices(i - 3)); % Start from 4 transitions ago
        endIdx = min(numBins, rapidTransitionIndices(i + 1)); % End after the next state
        isNoise(startIdx:endIdx) = true; % Mark this region as noise
        rapidAlternations = 0; % Reset count after marking noise
    end
end

% Replace noisy regions with NaN
viterbiStatesCleaned(isNoise) = NaN;

%%%%% Step 2: Merge states spanning <= 2 bins into previous state
for t = 2:numBins - 1
    if viterbiStatesCleaned(t) ~= viterbiStatesCleaned(t - 1) && ...
            viterbiStatesCleaned(t) ~= viterbiStatesCleaned(t + 1)
        viterbiStatesCleaned(t) = viterbiStatesCleaned(t - 1);
    end
end

%%%%% Step 3 Find UP events with duration <= 30ms (3 time bins) and surrounded by DOWN
upOnsets = find(diff([0; viterbiStatesCleaned == 1]) == 1);
upOffsets = find(diff([viterbiStatesCleaned == 1; 0]) == -1);

for i = 1:length(upOnsets)
    if upOffsets(i) - upOnsets(i) <= 3
        viterbiStatesCleaned(upOnsets(i):upOffsets(i))=0;
    end
end

% i= tvec_interp1(status==1);
% 
% for i = 1:max(epoch_id)
%     subplot(8,7,i)
%     histogram(spikeCounts((epoch_id==i)),100)
% end
% 
% %%%%% Step 3: Modify long UP states (more than 2 seconnds) with embedded short DOWN states
% minUpDuration = 200; % 2 seconds (200 bins at 10 ms/bin)
% 
% for startIdx = find(diff([0; viterbiStatesCleaned]) == 1)' % Find UP start indices
%     endIdx = find(diff([viterbiStatesCleaned; 0]) == -1, 1, 'first') + startIdx - 1;
% 
%     if endIdx - startIdx + 1 >= minUpDuration
%         % Check for short DOWN states within this UP state
%         downStarts = find(diff([0; viterbiStatesCleaned(startIdx:endIdx)]) == -1);
%         downEnds = find(diff([viterbiStatesCleaned(startIdx:endIdx); 0]) == 1);
% 
%         if isempty(downEnds) | isempty(downStarts) % if no short states
%             continue
%         end
% 
%         if length(downStarts) > length(downEnds)
%             downStarts = downStarts(1:length(downEnds));
%         end
% 
%         for d = 1:length(downStarts)
%             if downEnds(d) - downStarts(d) + 1 >= 2 % Short DOWN state (with 20ms or more)
%                 if d < length(downStarts) && (downStarts(d + 1) - downEnds(d) <= 3)
%                     % Merge into DOWN state
%                     viterbiStatesCleaned(startIdx + downEnds(d):startIdx + downStarts(d + 1) - 1) = 0;
%                 end
%             end
%         end
%     end
% end

%%%%% Step 4: Remove extremely short (1 bin events) or evens transiting out
%%%%% of and/or into nan regions.

% Transition indices
upOnsets = find(diff([0; viterbiStatesCleaned == 1]) == 1);
upOffsets = find(diff([viterbiStatesCleaned == 1; 0]) == -1);
downOnsets = find(diff([0; viterbiStatesCleaned == 0]) == 1);
downOffsets = find(diff([viterbiStatesCleaned == 0; 0]) == -1);

% Remove transitions into or out of NaN regions 
upOnsets = upOnsets(~isnan(viterbiStatesCleaned(upOnsets)));
upOffsets = upOffsets(~isnan(viterbiStatesCleaned(upOffsets)));
downOnsets = downOnsets(~isnan(viterbiStatesCleaned(downOnsets)));
downOffsets = downOffsets(~isnan(viterbiStatesCleaned(downOffsets)));

% Remove transitions with duration equal or less than 30ms
upOnsets = upOnsets(~isnan(viterbiStatesCleaned(upOnsets)));
upOffsets = upOffsets(~isnan(viterbiStatesCleaned(upOffsets)));
upDurations = (upOffsets - upOnsets);
upOnsets(upDurations<=3)=[];
upOffsets(upDurations<=3)=[];

downOnsets = downOnsets(~isnan(viterbiStatesCleaned(downOnsets)));
downOffsets = downOffsets(~isnan(viterbiStatesCleaned(downOffsets)));
downDurations = (downOffsets - downOnsets);
downOnsets(downDurations<=2)=[];
downOffsets(downDurations<=2)=[];


upOnsets = tvec_interp1(upOnsets)';
upOffsets =  tvec_interp1(upOffsets)';
downOnsets =  tvec_interp1(downOnsets)';
downOffsets =  tvec_interp1(downOffsets)';

% Pair UP events with immediately following DOWN events
upDownPairs = [];
for i = 1:length(upOffsets)
    % Find the first DOWN onset that comes after the current UP offset
    if upOffsets(i)-upOnsets(i)<=0.05 % only selecting UP with at least 50ms duration
        continue
    end

    nextDownIdx = find(downOnsets > upOffsets(i), 1, 'first');
    if ~isempty(nextDownIdx) 
        if downOnsets(nextDownIdx)-upOffsets(i)<=0.03 
            % Record the UP-DOWN pair: [UP onset, UP offset, DOWN onset, DOWN offset]
            upDownPairs = [upDownPairs; i, nextDownIdx]; %#ok<AGROW>
        end
    end
end

% Pair DOWN events with immediately following UP events
downUpPairs = [];
for i = 1:length(downOffsets)
    % Find the first UP onset that comes after the current DOWN offset
    nextUpIdx = find(upOnsets > downOffsets(i), 1, 'first');
    if ~isempty(nextUpIdx) 
        if upOnsets(nextUpIdx)-downOffsets(i)<=0.03 % Next event should happen within at least 30ms
            if upOffsets(nextUpIdx)-upOnsets(nextUpIdx)>=0.05 % only selecting UP with at least 50ms duration
                % Record the DOWN-UP pair: [DOWN onset, DOWN offset, UP onset, UP offset]
                downUpPairs = [downUpPairs; i, nextUpIdx]; %#ok<AGROW>
            end
        end
    end
end

[status,interval,index] = InIntervals(tvec_interp1,NREMInts);

slow_waves_markov.UP_ints = [upOnsets upOffsets];
slow_waves_markov.DOWN_ints = [downOnsets downOffsets];
slow_waves_markov.UP_DOWN_index = upDownPairs;
slow_waves_markov.DOWN_UP_index = downUpPairs;

slow_waves_markov.NREM_t = tvec_interp1(status);
slow_waves_markov.spike_count = spikeCounts(status);
slow_waves_markov.alpha_t = alpha_t;
slow_waves_markov.beta_t = beta_t;
slow_waves_markov.gamma_t = gamma_t;
slow_waves_markov.xi_t = xi_t;
slow_waves_markov.viterbi_states = viterbiStatesCleaned;

slow_waves_markov.p = P;
slow_waves_markov.alpha = alpha;
slow_waves_markov.beta = beta;
slow_waves_markov.mu = mu;
slow_waves_markov.beta = beta;


disp('EM Algorithm completed!');
disp(['Transition Matrix: ', mat2str(P)]);
disp(['Estimated alpha: ', num2str(alpha)]);
disp(['Estimated beta: ', num2str(beta)]);

% === Helper Functions ===
    function prob = emission_prob(yk, mu, alpha, beta, hist_range, spikeCounts)
        % Compute the emission probability
        if isempty(hist_range)
            lambda(1) = exp(mu + alpha);
            lambda(2) = exp(alpha);
        else
            lambda(1) = exp(mu + alpha + beta * sum(spikeCounts(max(1, hist_range))));% lamda for DOWN
            lambda(2) = exp(alpha + beta * sum(spikeCounts(max(1, hist_range))));% lamda for UPs
        end
        prob = exp(-lambda) .* (lambda.^yk) ./ factorial(yk);
    end

    function logL = compute_log_likelihood(x, param, gamma_t, spikeCounts, J, mu, alpha, beta,NREMInts, interval,epoch_id)
        % Log-likelihood for a given parameter
        logL = 0;

        for nEpoch = 1:size(NREMInts,1)
            tidx = find(interval==nEpoch); % tidx out of the entire session
            NREM_tidx = find(epoch_id == nEpoch); % tidx during sleep windows
            for t = J+1:length(NREM_tidx)
                S_t = gamma_t(NREM_tidx(t), :);

                if strcmp(param, 'alpha')
                    lambda_t(1) = exp(mu + x + beta * sum(spikeCounts(max(tidx(t)-J,1):tidx(t)-1)));
                    lambda_t(2) = exp(x + beta * sum(spikeCounts(max(tidx(t)-J,1):tidx(t)-1)));

                    logL = logL + (S_t(1) * (spikeCounts(tidx(t)) * log(lambda_t(1))) - lambda_t(1)) + (S_t(2) * (spikeCounts(tidx(t)) * log(lambda_t(2))) - lambda_t(2));
                elseif strcmp(param, 'mu')
                    lambda_t(1) = exp(x + alpha + beta * sum(spikeCounts(max(tidx(t)-J,1):tidx(t)-1)));
                    %                 lambda_t(2) = exp(alpha + beta * sum(spikeCounts(max(tidx(t)-J,1):tidx(t)-1)));
                    logL = logL + (S_t(1) * (spikeCounts(tidx(t)) * log(lambda_t(1))) - lambda_t(1));
                    %             elseif strcmp(param, 'beta')
                    %                 logL = logL + S_t * (mu + alpha * S_t + x * sum(spikeCounts(max(tidx(t)-J,1):tidx(t)-1)));
                end
            end
        end

    end

    function grad = compute_first_derivative(f, x)
        % Compute first derivative using finite differences
        h = 1e-5;
        grad = (f(x + h) - f(x - h)) / (2 * h);
    end

    function hess = compute_second_derivative(f, x)
        % Compute second derivative using finite differences
        h = 1e-5;
        hess = (f(x + h) - 2 * f(x) + f(x - h)) / (h^2);
    end

end
