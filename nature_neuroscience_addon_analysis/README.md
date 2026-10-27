# Nature Neuroscience revision analyses

Code for the analyses added during the revision of the V1–hippocampus sleep
reactivation paper. Each `main_XX_*.m` script reruns one analysis end to end
(computation and figures). Each script is split into sections, so you can run
the whole file or one section at a time.

## How to run

1. Put `VR_NPX_analysis` on the MATLAB path. Each main script then calls
   `addon_setup`, which adds this folder and its subfolders to the path.
2. The data folder is `D:\` or `P:\corticohippocampal_replay\nature_neuroscience_revision_analysis`.
   Everything below is written there, in the output folder listed for each script.
3. Open a main script, set its flags at the top, and run it:
   - `run_analysis = false` only redraws the figures that can be made from saved
     results.
   - `run_decoding` and `overwrite` control the slow per-session decoders. Sessions
     that already have an output file are skipped unless `overwrite = true`.

## Main scripts

| Script | Analysis | Output folders |
|---|---|---|
| `main_00_data_preparation` | Shared inputs: PRE ripple info and PSTH, capped log odds, non-duplicated events, spindle at ripple onset/peak, RUN1 decoding input. Run once; all steps are off by default. | revision-folder root (`ripple_info_PRE_unique.mat`, `KDE_reactivation_ripples_PSTH_*_unique.mat`, `RUN1_decoding_input.mat`, ...) |
| `main_01_reactivation_AUC` | Pooled reactivation coherence AUC (ripple power, spindle, SO phase, dominant SO trough), PRE / POST, equal-n, hierarchical bootstrap | `PRE`, `PRE_equal_n`, `POST_matched`, `POST_matched_equal_n`, `POST`, `POST_original_AUC_hierarchical` |
| `main_02_session_level_stats` | Per-session and per-animal AUC and geomean coherence; Low vs High, PRE vs POST, vs shuffle; abs bias; per-session / per-animal bootstrapped static and temporal AUC of the original POST decoder | `PRE_vs_POST`, `session_AUC_vs_shuffle`, `session_bootAUC_PRE_vs_POST*`, `session_bootGeomean*`, `POST_animal`, `abs_bias`, `POST_original_bootAUC_per_session`, `POST_original_bootAUC_per_animal` |
| `main_03_sleep_controls` | PRE vs POST: interhemispheric SO synchrony, spindle and ripple power | `SO_bilateral_synchrony_PRE_vs_POST`, `spindle_ripple_PRE_vs_POST` |
| `main_04_MUA_equalised_event_control` | Firing-rate control: Low/High events matched for HC and dominant V1 MUA (`MUA_equalised_event_control`), and additionally non-dominant V1 MUA (`MUA_equalised_event_control_with_nondominant`); geomean, AUC and V1 MUA PSTH | `<folder>\geomean`, `<folder>\AUC`, `<folder>\MUA_PSTH` |
| `main_05_V1_split_half` | Split-half V1 decoders: dominance half vs selector half | `split_half_V1_AUC_POST`, `split_half_V1_AUC_POST_hierarchical` (decoders in `C:\Data\PLS_split_half`) |
| `main_06_low_rate_decoder` | PLS-KDE decoders from neurons < 10 Hz | `PLS_KDE_rate_lt10`, `POST_rate_lt10`, `POST_rate_lt10_bootAUC_per_animal` |
| `main_07_Bayesian_decoding` | Bayesian decoding: RUN1 cross-validation and POST replay (all cells, < 10 Hz) | `Bayesian_RUN_decoding*`, `Bayesian_replay_decoding*` |
| `main_08_Lseq_decoding` | Sequence-likelihood (L_seq) decoding: RUN1 cross-validation and POST replay | `Lseq_RUN_decoding`, `Lseq_replay_decoding` |
| `main_10_pairwise_coupling` | Pairwise V1-HC neuron coupling across POST ripples vs shared track preference; gain control by leave-one-out total-V1 or own-hemisphere MUA | `pairwise_V1_HC_coupling\counts`, `pairwise_V1_HC_coupling\gain_total_V1`, `pairwise_V1_HC_coupling\gain_own_hemisphere` |

Dependencies between analyses:
- Everything uses the `main_00` inputs.
- `main_06` reuses `reactivation_AUC_coherence_PRE` from `main_01`.
- `main_08` replay scoring reads the events of `main_07`.

## Folder layout

| Folder | Contents |
|---|---|
| `shared\` | Helpers used by more than one analysis: AUC (`reactivation_coherence_auc`, `hemisphere_conditions`), log-odds capping, PLS-KDE decoder, RUN1 behaviour, permutation test, low-vs-high bar plots (`plot_AUC_low_vs_high*`, `AUC_low_high_*`) |
| `00_data_preparation\` ... `08_Lseq_decoding\` | The functions of the main script with the same number |
