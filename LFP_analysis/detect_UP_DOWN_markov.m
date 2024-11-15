
temp_05


%% Buz detection pipeline (Down detection is ok, but UP is really loose)
temp_05 = DetectSlowWaves_masa('time',tvec,'lfp',LFP(probe_no).best_V1(best_channel,:),'spikes',V1_clusters(probe_no),'NREMInts',behavioural_state_merged.SWS,'sensitivity',0.5);
% Maybe 


%% Filter the LFP: delta, high gamma and get MUA spike counts
SR = round(1/mean(diff(time)));
display('Filtering LFP')

data = lfp;
lfp = [];
lfp.data = data';

if size(time,1)<size(time,2)
    lfp.timestamps = time';
else
    lfp.timestamps = time;
end

lfp.samplingRate = SR;

deltaLFP = bz_Filter(lfp,'passband',filterparms.deltafilter,'filter','fir1','order',1);
deltaLFP.normamp = NormToInt(deltaLFP.data,'modZ',NREMInts,SR);
% deltaLFP.timestamps = time;
% deltaLFP.samplingRate = SR;
% filter_type  = 'bandpass';
% passband = filterparms.deltafilter;
% filter_order = round(6*SR/(max(passband)-min(passband)));  % creates filter for ripple
% norm_freq_range = passband/(SR/2); % SR/2 = nyquist freq i.e. highest freq that can be resolved
% b_SW = fir1(filter_order, norm_freq_range,filter_type);
% signal = filtfilt(b_SW,1,lfp);
% zscored_SW = zscore(abs(hilbert(signal)));


%    case false
gammaLFP= bz_Filter(lfp,'passband',filterparms.gammafilter,'filter','fir1','order',4);
%    case true
%         [ MUA ] = MUAfromDat( basePath,'channels',SWChan);
%         gammaLFP.amp = MUA.data;
%         gammaLFP.samplingRate = MUA.samplingRate;
%         gammaLFP.timestamps = MUA.timestamps;
% end
gammaLFP.smoothamp = smooth(gammaLFP.amp,round(filterparms.gammasmoothwin.*SR),'moving' );
gammaLFP.normamp = NormToInt(gammaLFP.smoothamp,'modZ',NREMInts,SR,'moving',filterparms.gammanormwin);



% Define observation sequences (spike rate, delta power, gamma power)
observations = [SpikeRate(:), DeltaPower(:), GammaPower(:)];

% Set up parameters for the two states
numStates = 2;  % UP and DOWN
states = {'UP', 'DOWN'};

% Define duration distributions for UP and DOWN states
meanDuration_UP = 150; % mean duration of UP state in time steps
stdDuration_UP = 50;   % standard deviation of UP state duration
meanDuration_DOWN = 50; % mean duration of DOWN state in time steps
stdDuration_DOWN = 20;  % standard deviation of DOWN state duration

% Define Gaussian observation model parameters for UP and DOWN states
mu_UP = [mean(SpikeRate(SpikeRate > prctile(SpikeRate, 60))), ...
         mean(DeltaPower(DeltaPower < prctile(DeltaPower, 30))), ...
         mean(GammaPower(GammaPower > prctile(GammaPower, 60)))];
sigma_UP = cov(observations(SpikeRate > prctile(SpikeRate, 60), :));

mu_DOWN = [mean(SpikeRate(SpikeRate < prctile(SpikeRate, 40))), ...
           mean(DeltaPower(DeltaPower > prctile(DeltaPower, 60))), ...
           mean(GammaPower(GammaPower < prctile(GammaPower, 40)))];
sigma_DOWN = cov(observations(DeltaPower > prctile(DeltaPower, 60), :));

% Define duration PDFs for UP and DOWN states
durationPDF_UP = @(d) exp(-((d - meanDuration_UP).^2) / (2 * stdDuration_UP^2)) / (stdDuration_UP * sqrt(2 * pi));
durationPDF_DOWN = @(d) exp(-((d - meanDuration_DOWN).^2) / (2 * stdDuration_DOWN^2)) / (stdDuration_DOWN * sqrt(2 * pi));

% Run simplified Viterbi-like inference with explicit durations
numTimeSteps = size(observations, 1);
logProb = -inf(numStates, numTimeSteps);  % Log probability table
statePath = zeros(1, numTimeSteps);       % Best state path

% Initialize at t=1, start with both UP and DOWN
logProb(1, 1) = log(mvnpdf(observations(1, :), mu_UP, sigma_UP));
logProb(2, 1) = log(mvnpdf(observations(1, :), mu_DOWN, sigma_DOWN));

% Forward pass with explicit duration modeling
for t = 2:numTimeSteps
    % Calculate for UP state (if we enter UP now, we previously were in DOWN)
    for d = 1:min(t, 100)  % limit max duration to 100 steps for efficiency
        % Observation probability over last d steps in UP state
        obsProb = sum(log(mvnpdf(observations(t-d+1:t, :), mu_UP, sigma_UP)));
        
        % Add duration probability
        durationProb = log(durationPDF_UP(d));
        
        % Total log probability for this segment
        totalLogProb_UP = logProb(2, t-d) + durationProb + obsProb;
        
        % Update if better path found
        if totalLogProb_UP > logProb(1, t)
            logProb(1, t) = totalLogProb_UP;
            statePath(t-d+1:t) = 1;  % Mark this segment as UP
        end
    end
    
    % Calculate for DOWN state (if we enter DOWN now, we previously were in UP)
    for d = 1:min(t, 100)
        % Observation probability over last d steps in DOWN state
        obsProb = sum(log(mvnpdf(observations(t-d+1:t, :), mu_DOWN, sigma_DOWN)));
        
        % Add duration probability
        durationProb = log(durationPDF_DOWN(d));
        
        % Total log probability for this segment
        totalLogProb_DOWN = logProb(1, t-d) + durationProb + obsProb;
        
        % Update if better path found
        if totalLogProb_DOWN > logProb(2, t)
            logProb(2, t) = totalLogProb_DOWN;
            statePath(t-d+1:t) = 2;  % Mark this segment as DOWN
        end
    end
end

% Decode final path
decodedStates = statePath;

% Visualize the inferred UP and DOWN states
figure;
plot(decodedStates, 'LineWidth', 1.5);
yticks([1 2]);
yticklabels({'UP', 'DOWN'});
title('Inferred UP and DOWN States with Explicit Durations');
xlabel('Time');
ylabel('State');