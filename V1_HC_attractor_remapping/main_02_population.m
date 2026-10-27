%% V1-HC attractor remapping, step 2: population track bias around lap transitions
% Needs results\data\epoch_table.mat (main_01_decode). For each decoder variant:
% time courses by lap in block, epoch summaries (grey bias toward the previous
% track, onset decoding on switch vs late laps), run-length dependence and the
% trial-by-trial chain; then a side-by-side comparison of decoder variants.
% Output: results\population\<decoder>\, results\population\decoder_comparison.*
P = ar_setup();

%% Session QC: cells and decoding accuracy at the track start
ar_session_qc();

%% All decoder variants
ar_population_analysis();

%% Decoder comparison
ar_decoder_comparison();

%% V1 with all animals (M24018 included)
ar_population_analysis('decoders', {'bayes_multinomial', 'pls'}, 'regions', {'V1'}, 'V1_all_animals', true, ...
    'out_name', 'population_V1_all_animals');

%% Raw log odds instead of z-scored (comparison)
ar_population_analysis('decoders', {'bayes_multinomial', 'bayes_multinomial_all', 'bayes_multinomial_lobo5', 'pls', 'pls_all'}, ...
    'measure', 'LO', 'out_name', 'population_rawLO', 'timecourses', false);
