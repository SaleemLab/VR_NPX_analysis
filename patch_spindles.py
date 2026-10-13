import re

file_path = r'c:\Users\masahiro.takigawa\Documents\GitHub\VR_NPX_analysis\LFP_analysis\DetectSpindles_masa.m'
with open(file_path, 'r') as f:
    code = f.read()

# 1. Input parser modifications
code = code.replace(
    "addParameter(p, 'savepath',      [],           @isstr);\nparse(p, varargin{:});",
    "addParameter(p, 'savepath',      [],           @isstr);\naddParameter(p, 'method',        'new',        @isstr);\nparse(p, varargin{:});"
)

code = code.replace(
    "best_channel    = p.Results.best_channel;",
    "best_channel    = p.Results.best_channel;\nmethod          = p.Results.method;\n\nif strcmp(method, 'default') && ismember('thresholds', p.UsingDefaults)\n    lowThresh = 2.5;\n    highThresh = 5;\nend"
)

# 2. Envelope modifications
env_old = """%% 3. Compute envelopes
rms_window = round(frequency / 5);  % ~200 ms window

sigma_envelope    = envelope(sigma_signal, rms_window, 'rms');
broadband_envelope = envelope(broadband_signal, rms_window, 'rms');

% Smooth the sigma envelope
sigma_envelope = smoothdata(sigma_envelope, 'gaussian', rms_window);

% Z-score the sigma envelope for thresholding
zscored_sigma = zscore(sigma_envelope);"""

env_new = """%% 3. Compute envelopes
rms_window = round(frequency / 5);  % ~200 ms window

if strcmp(method, 'new')
    sigma_envelope    = envelope(sigma_signal, rms_window, 'rms');
    % Smooth the sigma envelope
    sigma_envelope = smoothdata(sigma_envelope, 'gaussian', rms_window);
    % Z-score the sigma envelope for thresholding
    zscored_sigma = zscore(sigma_envelope);
else % 'default' (classic FMA FindSpindles logic)
    % Z-score the filtered signal
    z_sigma = zscore(sigma_signal);
    % Squared Hilbert envelope
    squared_envelope = abs(hilbert(z_sigma)).^2;
    % 320 ms moving average smooth
    windowLength = round(frequency * 0.320);
    smooth_win = ones(windowLength, 1) / windowLength;
    zscored_sigma = filtfilt(smooth_win, 1, squared_envelope);
    sigma_envelope = sqrt(zscored_sigma); % just for relative power calculation visualization
end

broadband_envelope = envelope(broadband_signal, rms_window, 'rms');"""

code = code.replace(env_old, env_new)

# 3. Relative power skip
rel_old = """valid_rel = mean_rel_power >= min_rel_power;
secondPass = secondPass(valid_rel, :);"""

rel_new = """if strcmp(method, 'new')
    valid_rel = mean_rel_power >= min_rel_power;
else
    valid_rel = true(size(secondPass, 1), 1); % Skip relative power check in default mode
end
secondPass = secondPass(valid_rel, :);"""

code = code.replace(rel_old, rel_new)

# 4. Zero crossings skip
cyc_old = """valid_cycles = n_cycles_vec >= min_cycles;
secondPass = secondPass(valid_cycles, :);"""

cyc_new = """if strcmp(method, 'new')
    valid_cycles = n_cycles_vec >= min_cycles;
else
    valid_cycles = true(size(secondPass, 1), 1); % Skip cycle check in default mode
end
secondPass = secondPass(valid_cycles, :);"""

code = code.replace(cyc_old, cyc_new)

with open(file_path, 'w') as f:
    f.write(code)

print("Patch applied successfully.")
