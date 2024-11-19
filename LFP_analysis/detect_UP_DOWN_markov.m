function detect_UP_DOWN_markov(tvec,slow_waves,spiketimes,LFP)


%% Buz detection pipeline (Down detection is ok, but UP is pretty bad)
temp_05 = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.5);
% Maybe

best_channel = find(LFP(probe_no).best_V1_channel==slow_waves(nprobe).best_channel);

if isempty(best_channel)
    [~,best_channel] = max(LFP(nprobe).best_V1_high_freq_power(:,7));
    temp = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1_high_freq(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.5);
    [spindles(probe_no)] = FindSpindles_masa(LFP(probe_no).best_V1_high_freq(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
else
    temp_05 = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.5);
    [spindles(probe_no)] = FindSpindles_masa(LFP(probe_no).best_V1(best_channel,:),LFP(probe_no).tvec','behaviour',Behaviour,'durations',[400 3000],'frequency',mean(1./diff(LFP(nprobe).tvec)),...
        'noise',[],'passband',[9 17],'thresholds',[1 3],'show','off');
end

% PSD slope quantification using fooof ()
tvec = LFP(1).tvec;
SR = round(1/mean(diff(tvec)));
nfft_seconds= 2;
nfft = 2^(nextpow2(SR*nfft_seconds));
win  = hanning(nfft);

clipDur = 10; % seconds
timebin_edges = tvec(1):10:tvec(end); % 10 seconds timebin edges for PSD slope
nClipSamps = round(SR*clipDur);
PSD_slope=[];

disp('PSD slope for slow waves started')


%%%%%%%%%%%%%%%%% Hidden Markov UP DOWN = ;
NREMInts = behavioural_state_merged.SWS;

% Filter the LFP: delta, high gamma and get MUA spike counts
filterparms.deltafilter = [0.5 8];%heuristically defined.  room for improvement here.
filterparms.gammafilter = [100 400];
filterparms.gammasmoothwin = 0.08; %window for smoothing gamma power (s)
filterparms.gammanormwin = 20; %window for gamma normalization (s)

SR = round(1/mean(diff(tvec)));
display('Filtering LFP')

lfp = [];
lfp.data = LFP(probe_no).best_V1(best_channel,:)';

if size(tvec,1)<size(tvec,2)
    lfp.timestamps = tvec';
else
    lfp.timestamps = tvec;
end

lfp.samplingRate = SR;

%%% Delta power
deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
deltaLFP.normamp = NormToInt(deltaLFP.data,'modZ',NREMInts,SR);

%%% Gamma power
gammaLFP= bz_Filter(lfp,'passband',filterparms.gammafilter,'filter','fir1','order',4);
gammaLFP.smoothamp = smooth(gammaLFP.amp,round(filterparms.gammasmoothwin.*SR),'moving');
gammaLFP.normamp = NormToInt(gammaLFP.smoothamp,'modZ',NREMInts,SR,'moving',filterparms.gammanormwin);

%%% Spikecount smoothed with σ = 30 ms
tvec_edges = [tvec(1)-1/(1/mean(diff(tvec))*2) tvec+1/(1/mean(diff(tvec))*2)];
if ~isempty(behavioural_state_merged.SWS)
    sleep_tvec = Restrict(tvec,behavioural_state_merged.SWS);
    sleep_index = ismember(tvec,sleep_tvec);
end

HPC_spike_counts=[];
V1_spike_counts=[];
w = gausswin(0.03*1/mean(diff(tvec))); % Smoothed with σ = 30 ms
w = w / sum(w);

V1_spike_count =[];
V1_spike_count = filtfilt(w,1,histcounts(V1_clusters(nprobe).spike_times,tvec_edges)')';
%             V1_spike_count = V1_spike_count-mean(V1_spike_count(sleep_index==1))./std(V1_spike_count(sleep_index==1));
%             V1_spike_count = NormToInt(V1_spike_count','modZ',NREMInts,SR)';
%             V1_MUA_sleep= V1_spike_count(sleep_index==1);

% Define observation sequences (spike rate, delta power, gamma power)
tvec_interp1 = tvec(1):0.001:tvec(end);

observations = [interp1(tvec,V1_spike_count,tvec_interp1)'./0.001,interp1(tvec,deltaLFP.normamp,tvec_interp1)',interp1(tvec,gammaLFP.normamp,tvec_interp1)'];

% Parameters for the detection algorithm
T = 1;                    % Bin size in ms
J = 20;                   % Number of history bins (corresponds to 20 ms memory)
beta = 0.1;               % History-dependence weight adjusted for 1 ms bins
P_DU = 0.9;               % Transition probability from DOWN to UP
P_UD = 0.9;               % Transition probability from UP to DOWN
P_DD = 0.1;               % Self-transition in DOWN state
P_UU = 0.1;               % Self-transition in UP state

[~, ~, ~, ~, ~, binnedArray] =psthAndBA(V1_clusters(nprobe).spike_times,slow_waves(nprobe).ints.UP(:,1),[0 0.1],0.01);
alpha =                 % Rate during UP periods (for spikes)
mu = -2;                  % Rate difference during DOWN and UP periods (for spikes)

% Observation model parameters for LFP (delta and gamma power)
% (Adjust means and variances based on observed data characteristics)
meanDelta_UP = 0.2;       % Expected delta power (low) during UP
stdDelta_UP = 0.05;       % Std of delta power during UP
meanDelta_DOWN = 1.0;     % Expected delta power (high) during DOWN
stdDelta_DOWN = 0.2;      % Std of delta power during DOWN
meanGamma_UP = 1.0;       % Expected gamma power (high) during UP
stdGamma_UP = 0.2;        % Std of gamma power during UP
meanGamma_DOWN = 0.2;     % Expected gamma power (low) during DOWN
stdGamma_DOWN = 0.05;     % Std of gamma power during DOWN

% Load or generate example data for spike counts and LFP power bands
% Assume spikeCounts, deltaPower, and gammaPower are vectors of the same length
% where each entry corresponds to a 1 ms bin.
spikeCounts = poissrnd(alpha, 10000, 1);  % Example Poisson-distributed spike counts
deltaPower = normrnd(meanDelta_DOWN, stdDelta_DOWN, 10000, 1);  % Example delta power
gammaPower = normrnd(meanGamma_DOWN, stdGamma_DOWN, 10000, 1);  % Example gamma power

% Initialize state estimate vector
numBins = length(spikeCounts);
stateEstimate = zeros(numBins, 1);  % 1 for UP, 0 for DOWN
logProbUP = log(P_UU);              % Initialize log probability for UP state
logProbDOWN = log(P_DD);            % Initialize log probability for DOWN state

% Iterate over each bin to estimate the state based on the history-dependent model
for t = 1:numBins
    % Calculate emission probabilities
    % Spike count probability
    rateUP = alpha;
    rateDOWN = alpha + mu;
    emissionProbSpikeUP = poisspdf(spikeCounts(t), rateUP);
    emissionProbSpikeDOWN = poisspdf(spikeCounts(t), rateDOWN);

    % Delta and gamma power probabilities (Gaussian for continuous LFP)
    emissionProbDeltaUP = normpdf(deltaPower(t), meanDelta_UP, stdDelta_UP);
    emissionProbDeltaDOWN = normpdf(deltaPower(t), meanDelta_DOWN, stdDelta_DOWN);
    emissionProbGammaUP = normpdf(gammaPower(t), meanGamma_UP, stdGamma_UP);
    emissionProbGammaDOWN = normpdf(gammaPower(t), meanGamma_DOWN, stdGamma_DOWN);

    % Combine emission probabilities for UP and DOWN states
    emissionProbUP = emissionProbSpikeUP * emissionProbDeltaUP * emissionProbGammaUP;
    emissionProbDOWN = emissionProbSpikeDOWN * emissionProbDeltaDOWN * emissionProbGammaDOWN;

    % Update log probabilities with transition and combined emission probabilities
    if t > 1
        % Compute log probabilities based on previous state and history
        historyFactor = exp(-beta * abs(t - J));  % History factor for smoothing

        if stateEstimate(t-1) == 1  % Previous state was UP
            logProbUP = log(P_UU * historyFactor) + log(emissionProbUP);
            logProbDOWN = log(P_UD * historyFactor) + log(emissionProbDOWN);
        else  % Previous state was DOWN
            logProbUP = log(P_DU * historyFactor) + log(emissionProbUP);
            logProbDOWN = log(P_DD * historyFactor) + log(emissionProbDOWN);
        end
    end

    % Choose state with higher probability
    if logProbUP > logProbDOWN
        stateEstimate(t) = 1;  % UP state
    else
        stateEstimate(t) = 0;  % DOWN state
    end
end

% Merge adjacent bins to get UP and DOWN periods
upPeriods = [];
downPeriods = [];
stateChangePoints = find(diff([0; stateEstimate; 0]));

for i = 1:2:length(stateChangePoints)-1
    tOn = stateChangePoints(i);       % Onset time of UP or DOWN state
    tOff = stateChangePoints(i+1) - 1; % Offset time of UP or DOWN state

    % Check if UP or DOWN period
    if stateEstimate(tOn) == 1
        upPeriods = [upPeriods; [tOn, tOff, tOff - tOn + 1]];
    else
        downPeriods = [downPeriods; [tOn, tOff, tOff - tOn + 1]];
    end
end

% Display UP and DOWN period information
disp('UP periods (start, end, duration):');
disp(upPeriods);
disp('DOWN periods (start, end, duration):');
disp(downPeriods);

% Plot state estimates
figure;
plot(stateEstimate, 'LineWidth', 1.5);
title('Estimated UP and DOWN States');
xlabel('Time (ms)');
ylabel('State (1 = UP, 0 = DOWN)');


