%% Continuous V1 SO phase and spindle-band amplitude per hemisphere, POST
% Recomputes the traces that main_LFP_analysis_sleep.m used for the ripple-onset
% SO phase / spindle amplitude (they were not saved), on the reference V1 shank of
% each hemisphere (best SO channel, as add_spindle_onset_peak_to_ripple_info):
%   SO phase      : bz_Filter 0.5-4 Hz (fir1, order 1), phase
%   spindle amp   : fir1 9-17 Hz (order 6*SR/8), filtfilt, |hilbert|, z-scored
% on extracted_LFP.mat average_V1 (both probes, truncated to the shorter
% best_HPC as in the original). Check: values at the nearest sample to every
% POST ripple onset are compared with ripples_all_POST (SO_phase_ripple_onset,
% spindle_amplitude_ripple_onset).
% Saved at 10 ms resolution (nearest sample; phase is not averaged):
%   nature_neuroscience_revision_analysis\SO_spindle_traces_POST.mat
%   traces(s).t (1 x n), .SO_phase (n x 2), .spindle_amp (n x 2), columns = hemisphere 1, 2
% sessions (optional): subset to run as a check; the file is saved only for all 22.
function extract_SO_spindle_traces_POST(sessions)
    rev_dir = addon_setup();
    analysis_folder = fileparts(rev_dir);
    out_file = fullfile(rev_dir, 'SO_spindle_traces_POST.mat');
    dt = 0.01;

    load(fullfile(analysis_folder, 'slow_waves_all_POST.mat'), 'slow_waves_all');
    load(fullfile(rev_dir, 'ripples_all_POST.mat'), 'ripples_all');
    SUBJECTS = {'M24016', 'M24017', 'M24018', 'M24062', 'M24064', 'M24065'};
    experiment_info = subject_session_stimuli_mapping(SUBJECTS, 'bilateral');
    experiment_info = experiment_info([4 5 6 17 18 19 21 33 34 35 44 45 46 47 56 58 59 60 70 71 72 73]);
    if nargin < 1, sessions = 1:numel(experiment_info); end

    traces = struct('t', {}, 'SO_phase', {}, 'spindle_amp', {}, 'ref_shank', {}, 'check', {});
    for nsession = sessions
        tic
        is_post = strcmp(experiment_info(nsession).StimulusName, 'SleepChronic');
        session_info = experiment_info(nsession).session(is_post);
        load(fullfile(session_info.probe(1).ANALYSIS_DATAPATH, 'extracted_LFP.mat'), 'LFP');

        % Reference shank per hemisphere (index into the concatenated V1 shanks)
        ref = NaN(1, 2);
        for probe_no = 1:2
            sw = slow_waves_all(probe_no);
            ref(probe_no) = find(sw.shank_id{nsession} == sw.shank{nsession}(sw.channel{nsession} == sw.best_channel(nsession)) ...
                & sw.probe_hemisphere{nsession} == probe_no);
        end

        n = min(length(LFP(1).best_HPC), length(LFP(2).best_HPC));
        tvec = LFP(1).tvec(1:n);
        data = [LFP(1).average_V1(:, 1:n)' LFP(2).average_V1(:, 1:n)'];
        lfp = struct('samplingRate', round(1 / mean(diff(LFP(1).tvec))), 'timestamps', tvec, 'data', data(:, ref));
        clear LFP data

        delta = bz_Filter(lfp, 'passband', [0.5 4], 'filter', 'fir1', 'order', 1);
        SO_phase = delta.phase;
        clear delta
        passband = [9 17];
        b_spindle = fir1(round(6 * lfp.samplingRate / (max(passband) - min(passband))), passband / (lfp.samplingRate / 2), 'bandpass');
        signal = zeros(size(lfp.data));
        for c = 1:2
            signal(:, c) = filtfilt(b_spindle, 1, lfp.data(:, c));
        end
        spindle_amp = zscore(abs(hilbert(signal)));
        clear signal

        % Check against the saved ripple-onset values
        chk = NaN(2, 2);   % rows hemisphere, cols [max |phase diff|, max |amp diff|]
        for probe_no = 1:2
            on = ripples_all(probe_no).onset(ripples_all(probe_no).session_count == nsession);
            tidx = interp1(tvec, 1:n, on, 'nearest', 'extrap');
            ph_saved = ripples_all(probe_no).SO_phase_ripple_onset{nsession}(ref, :)';
            amp_saved = ripples_all(probe_no).spindle_amplitude_ripple_onset{nsession}(ref, :)';
            dph = angle(exp(1i * (SO_phase(tidx, :) - ph_saved)));
            chk(probe_no, :) = [max(abs(dph(:)), [], 'omitnan') max(abs(spindle_amp(tidx, :) - amp_saved), [], 'all', 'omitnan')];
        end

        t = tvec(1):dt:tvec(end);
        k = interp1(tvec, 1:n, t, 'nearest');
        traces(nsession).t = t;
        traces(nsession).SO_phase = single(SO_phase(k, :));
        traces(nsession).spindle_amp = single(spindle_amp(k, :));
        traces(nsession).ref_shank = ref;
        traces(nsession).check = chk;
        fprintf('Session %d: ref shanks %s, max |diff| vs saved onset values: phase %.2g / %.2g rad, spindle %.2g / %.2g z (%.0f s)\n', ...
            nsession, mat2str(ref), chk(1, 1), chk(2, 1), chk(1, 2), chk(2, 2), toc);
        clear SO_phase spindle_amp lfp
    end
    if isequal(sessions, 1:numel(experiment_info))
        save(out_file, 'traces', 'dt', '-v7.3');
        fprintf('Saved %s\n', out_file);
    end
end
