# V1–HC attractor remapping

How quickly hippocampal (HPC) and V1 populations switch to the representation of the current track, and whether subjective expectation biases them. Mice run Track L (T1) and Track R (T2) in blocks of 7–13 laps, so the same track is expected within a block, and the first lap after a block switch is a surprise. Between laps the screen is grey for 4–6 s (no track information).

Data: `D:\corticohippocampal_replay\V1_HC_attractor_remapping` (`inputs\` holds copies of the revision-folder RUN1 data and is never edited; `results\` holds everything written here). The 22 sessions and 6 animals are the paper's familiar sessions (same order as the revision folder).

## Questions
1. Is decoding of the track at onset (0–0.5 s) worse on the first lap of a block than late in the block (laps 6+)?
2. Is there a spontaneous track bias in the grey screen (early: first 2 s after the lap ends; late: last 2 s before the next onset), toward the track of the block?
3. Does the grey bias grow with the number of laps of the same track just run?
4. Trial by trial, does end-of-track bias predict grey bias, and grey bias predict onset bias?
5. The same in V1. Main V1 analyses exclude M24018 (V1_L has 3–15 cells).
6. Later: theta-cycle "flicker" between maps at track onset (needs the server for LFP; theta phase will come from HPC spatial-cell MUA, and LFP only for validation).

## Pipeline
| Script | What it does | Output |
|---|---|---|
| `main_01_decode` | Lap / block tables (licks on ephys clock), whole-session decoding, epoch table | `results\data\lap_tables.mat`, `results\decoding\session_XX.mat`, `results\data\epoch_table.mat` |
| `main_02_population` | Session QC (HPC coverage of the track start), population analyses and figures for every decoder variant, decoder comparison, V1 with M24018 | `results\qc\`, `results\population\<decoder>\`, `results\population\decoder_comparison.*`, `results\population_V1_all_animals\` |
| `main_03_single_cell` | Single-cell grey and onset firing vs track preference, switch step, field-location split; V1_L / V1_R population | `results\single_cell\`, `results\population_V1_hemispheres\` |
| `main_04_trial_dynamics` | Late grey -> onset (0-200, 0-100, 100-200 ms) trial LME, also on raw log odds (`ar_trial_prediction`); continuity across onset vs within grey (`ar_onset_reset`); HPC-V1 grey coupling and its effect on the prediction (`ar_grey_coupling`); lead/lag inside the grey only (`ar_grey_lag`) | `results\trial_prediction*\`, `results\grey_coupling\`, `results\grey_lag\` |
| `main_05_theta_and_report` | Theta from HPC MUA: spectra (`ar_theta_mua`), phase and cycles (`ar_theta_phase`, Vollan et al. 2025 method), theta-cycle flicker around onset with per-session examples (`ar_theta_flicker`); report figures and HTML (`ar_report_figures`, `ar_build_report`, `report_template.html`) | `results\theta_mua\`, `results\theta_phase\`, `results\theta_flicker\`, `results\report\`, `results\V1_HC_remapping_report.html` |

Shared helpers: `ar_bin_session` (event-locked 100 ms bins), `ar_bayes_lap_stats`, `ar_rate_maps`, `ar_bayes_track`, `ar_balanced_train_laps` (Bayesian decoder pieces), `ar_trial_table`, `ar_trial_lme`, `ar_lme_coef`, `ar_stats`, `ar_save_figure`.

Results so far: [RESULTS_SUMMARY.md](RESULTS_SUMMARY.md).

## Decoders (`ar_decode_session.m`)
- Grid: 100 ms windows stepped by 50 ms over the whole session (laps and grey screens).
- Templates come only from in-lap data (track onset to offset). Grey-screen windows have no position and are never used for training.
- **Bayesian** (multinomial = main, Poisson = check): 2 tracks × 14 10 cm bins, running samples > 5 cm/s, flat prior. Track log odds = posterior mass on T1 vs T2, with position marginalised.
- **PLS-KDE** (paper decoder, location agnostic): PLS (3 components) on the good positions, then a 3-D KDE per track. Good positions and bandwidth are chosen once per session (`ar_pls_setup.m`, 500 shuffles). KDE log densities are computed in log space (`ar_log_kde.m`), so log odds are never infinite.
- **Cross-validation by held-out blocks.** The grey screen before lap k and lap k itself are decoded by a decoder trained without block b = block(k) and its neighbours. Training laps are equalised between tracks (subsampled).
  - `lobo3` (blocks b−1..b+1 held out): the nearest same-track training blocks are closer in time than the other track's, so slow drift would favour the block's track (liberal).
  - `lobo5` (b−2..b+2 held out): the other track's blocks are closer, so drift would favour the other track (conservative). Effects present in both cannot come from drift.
  - `_all`: trained on all laps. Only its grey-screen output is out of sample.
- Measures: zLO = log odds z-scored over the session, clipped at the 0.5/99.5 percentiles, then centred per epoch by the average of the previous-T1 and previous-T2 means, and signed toward the previous or upcoming track.

## Drift and the switch-locked test
Tracks alternate in blocks, so slow population drift is confounded with track: decoders trained on blocks far away in time lose much of their accuracy in some sessions (lobo5 vs lobo3), and grey-screen "bias toward the block's track" flips sign between lobo3 and lobo5. The drift-robust test (`switch_locked` in `ar_population_analysis.m`) compares, around every block switch, the transition after the first lap of the new track (i = 1) with the switch transition (i = 0): the two are about one lap apart and decoded by the same decoder, and the within-block change between consecutive transitions is the drift baseline. The step is also computed without lick transitions, after regressing out spike count and speed, and for spike count, speed and licks themselves.

## Statistics
Session means, then a signed-rank test across sessions. The number of sessions and animals in the predicted direction is reported. An LME on transitions with random intercepts for animal and session is also fitted. Animal means are shown in every figure.
