%% V1-HC attractor remapping, step 1: lap tables and whole-session track decoding
% Builds the lap / block / transition tables (licks on the ephys clock) and
% decodes track in 100 ms non-overlapping bins (a whole-session grid plus bins
% locked to every lap end and onset) with block-held-out Bayesian (multinomial,
% Poisson) and PLS-KDE decoders (ar_decode_session.m), then summarises epochs
% around every lap transition (ar_epoch_table.m).
% PLS setups (good positions, bandwidth; the slow shuffle step) are reused from
% the first run (results\_v1_overlap50ms\decoding), whose PLS training bins are
% the same 100 ms bins.
% Output: results\data\lap_tables.mat, results\decoding\session_XX.mat,
%         results\data\epoch_table.mat
P = ar_setup();
build_tables = false;
sessions = 1:22;
overwrite = false;
old_dir = fullfile(P.results, '_v1_overlap50ms', 'decoding');

%% Lap tables
if build_tables || ~exist(fullfile(P.results, 'data', 'lap_tables.mat'), 'file')
    ar_build_lap_tables();
end

%% Decoding
S = load(fullfile(P.inputs, 'RUN1_decoding_input.mat'), 'D'); D = S.D; clear S
Lt = load(fullfile(P.results, 'data', 'lap_tables.mat'), 'L'); L = Lt.L;
out_dir = fullfile(P.results, 'decoding');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if isempty(gcp('nocreate')), parpool(6); end
for s = sessions
    f = fullfile(out_dir, sprintf('session_%02d.mat', s));
    if exist(f, 'file') && ~overwrite, fprintf('Session %d exists, skipping\n', s); continue; end
    fprintf('\n=== Session %d: %s %s ===\n', s, D(s).subject, D(s).session);
    old_setups = struct();
    fo = fullfile(old_dir, sprintf('session_%02d.mat', s));
    if exist(fo, 'file')
        O = load(fo, 'out');
        for r = 1:numel(O.out.region)
            if ~isempty(O.out.region(r).pls.setup), old_setups.(O.out.region(r).name) = O.out.region(r).pls.setup; end
        end
    end
    tic
    out = ar_decode_session(D(s), L(s).laps, s, P, old_setups);
    save(f, 'out', '-v7.3');
    fprintf('  saved (%.0f s)\n', toc);
end

%% Epoch table
ar_epoch_table(1:22);
