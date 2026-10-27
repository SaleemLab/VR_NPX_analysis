%% V1-HC attractor remapping, step 3: single cells and V1 hemispheres
% Single cells (ar_single_cell.m): grey-screen and onset firing vs in-track
% preference, the drift-robust switch step (transition after the first lap of a
% new block minus the switch transition) vs a within-block baseline, and a split
% by field location. Then the population analysis for V1_L and V1_R separately.
% Output: results\single_cell\, results\population_V1_hemispheres\
P = ar_setup();

%% Single cells
ar_single_cell();

%% V1 hemispheres (population)
ar_population_analysis('decoders', {'bayes_multinomial', 'bayes_multinomial_all', 'pls'}, ...
    'regions', {'V1_L', 'V1_R'}, 'out_name', 'population_V1_hemispheres');
