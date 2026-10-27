# V1–HC attractor remapping: results (v2, 2026-10-24)

The report with figures is `results\V1_HC_remapping_report.html`.

**Data.** 22 sessions, 6 animals, RUN1. 1,716 lap-to-lap transitions with a 4–6 s grey screen (uniform; median 4.98 s), 175 of them block switches. HPC uses all 6 animals. V1 excludes M24018 (5 animals, 19 sessions).

**v2 binning.** Non-overlapping 100 ms bins, with edges exactly at lap end and at track onset (`ar_bin_session`), plus a whole-session grid.

| Epoch | Window |
|---|---|
| end_prev | last 1 s of the lap |
| grey_early | first 1 s of grey |
| grey_late | last 1 s of grey |
| onset | 0–200 ms (also 0–100 and 100–200 ms) |
| onset_late | 0.2–1 s |
| lap_rest | from 1 s into the lap |

The v1 outputs (100 ms windows stepped by 50 ms, 2 s grey epochs) are in `results\_v1_overlap50ms\`.

**Decoders.**
- Multinomial Bayes (main), Poisson, and PLS-KDE.
- Held-out blocks: lobo3 (b±1) and lobo5 (b±2), plus all-laps (grey epochs only).
- Values are z-scored log odds (zLO) unless stated.
- Statistics: session means with signed-rank tests, sessions and animals > 0, and LMEs with animal and session (and block) random intercepts.

## 0. Quality, drift, grey duration
- **Lap timing.** Laps end at track offset (position ~139 cm) and start at onset (1–2 cm 50 ms later).
- **Excluded transitions.** 68 transitions (3.8%) have no grey: the lap ended early and the next started at once. One has a 9.4 s gap. All of these are excluded.
- **HPC coverage of the track start** is adequate in all sessions. HPC onset accuracy is lowest in M24018 and M24062 (0.66–0.67) and highest in M24064 (0.86).
- **Slow drift is large.**
  - Session 2 in-lap accuracy goes from 0.80 (HPC) / 0.98 (V1) with lobo3 to 0.51 / 0.74 with lobo5.
  - Switch and non-switch grey biases swap places between lobo3 and lobo5.
  - Averaging the two cancels the drift: HPC early grey +0.22 to +0.28 (Bayes) and +0.12 to +0.14 (PLS); HPC late grey +0.01 to +0.07; V1 late ≈ 0.

## 1. Track onset: no cost of expectation; the first lap is enhanced (mostly V1)
- **V1, first lap vs laps 6+.** At 0–200 ms the difference is +0.35 (17/19 sessions, 5/5 animals). After regressing out spikes and speed it is +0.19 (p = 0.018, 15/19 sessions, 5/5 animals); PLS gives +0.16 (p = 0.02).
- **V1, lap 2 − lap 1.** −0.16 at 0–200 ms (p = 0.044) and −0.16 at 0.2–1 s (p = 0.02).
- **Firing.** V1 fires +17 spikes per 100 ms bin more at onset on first laps (18/19 sessions).
- **Single V1 cells.** Cells tuned to the new track increase on the first lap: ρ −0.25 in V1_R (0/19 sessions and 0/5 animals positive) and −0.11 in V1_L (p = 0.008).
- **Decline over the block.** V1 onset coding declines with laps in the block (p = 2e-5).
- **HPC.**
  - 0–200 ms: no difference between first and late laps (−0.01).
  - 0.2–1 s: better on the first lap, +0.22 (p = 0.012), with p(true track) +0.06 (18/22 sessions) and PLS +0.16 after controls (p = 0.005).
  - Firing at 0.2–1 s is +8 per bin on first laps (20/22 sessions, 6/6 animals).
- **Track coding at onset.** It is weak in the first 100 ms and arrives at 100–200 ms (session 1 HPC p(true track): 0.59 → 0.77).

## 2. Grey screen: HPC keeps the just-run track ~1 s; V1 only via its offset response
- **Switch step** (grey after the first new lap minus pre-switch grey, all-laps decoder):
  - HPC 0–1 s: +0.61 (p = 4e-4, 17/22 sessions, 5/6 animals; raw +4.8 nats).
  - HPC 0.2–1 s: +0.52 (p = 6e-4, 18/22 sessions). Held-out Bayes gives +0.42 and PLS +0.20 (p = 0.014).
  - It decays: +0.97 (0–0.2 s) → +0.66 (0.2–0.5 s) → +0.44 (0.5–1 s).
  - It survives removing lick transitions (+0.50) and regressing out spikes and speed (+0.47).
  - By animal (0.2–1 s): 0.25 / 1.00 / 0.30 / −0.06 / 0.90 / 0.62. M24062 is null.
- **V1 step.** +0.40 at 0–0.2 s (19/19 sessions), +0.12 at 0.2–0.5 s (marginal), ≈ 0 from 0.5 s. This is visual latency or offset response.
- **Single HPC cells.** Step vs track preference: ρ 0.17 (22/22 sessions, 6/6 animals). End-of-track cells 0.29 (p = 4e-4, 6/6 animals); mid-track 0.12 (p = 0.035); start cells ns. The within-block baseline is ≈ 0.
- **Late grey (last 1 s).** HPC step +0.10, ns. No build-up with run length in either region.

## 3. Late grey predicts onset, trial by trial (`ar_trial_prediction`)
- **Model.** LME: onset ~ late grey + track pair + lap in block + spikes and speed (main effects and signed by track) + licks in both epochs, with animal, session and block random intercepts.
- **Results:**

  | | HPC Bayes | HPC PLS | V1 Bayes | V1 PLS |
  |---|---|---|---|---|
  | β, 0–200 ms | 0.82 | 0.72 | 0.79 | 0.65 |
  | β, 0–100 ms | 0.77 | 0.71 | 0.77 | 0.65 |
  | β, 100–200 ms | 0.90 | 0.76 | 0.79 | 0.66 |
  | β, 0.2–1 s | 0.51 | 0.38 | 0.52 | 0.40 |
  | Per-session partial Spearman | 0.41 (22/22) | 0.33 (22/22) | 0.27 (17/19) | 0.22 (18/19) |

  Every animal is positive; the weakest is M24065 in HPC (0.16).
- **Temporal control.** The previous or next trial's late grey gets β −0.05 to 0.15. Early grey and lap end add little.
- **Time-resolved.** The correlation with the 0–200 ms bias rises from ~0.15 at −1.5 s to 0.3 in the last 200 ms. HPC carries it at least 1.5 s after onset; V1 about 300 ms.
- **Not robust (dropped).** V1 late grey → HPC onset: β 0.11–0.14 with z, but 0.02–0.05 (ns) with raw log odds.

## 4. HPC–V1 coupling in the grey (`ar_grey_coupling`, `ar_grey_lag`)
- **Design.** Lags are computed within the grey window only. The first 300 ms after lap end are excluded (V1 offset response).
- **Within trials.** Bin by bin, r ≤ 0.017 (null ≈ 0).
- **Between trials.** ≈ 0, slightly negative in late grey (−0.08; −0.02 to −0.04 with controls).
- **Whole-grey lead/lag** (0.3 s after lap end → onset, 37–57 bins per trial): lagged r minus null is within ±0.015 at all lags, for all decoders and with controls. The asymmetry index ≈ 0.
- **Modulation.** Per-trial coupling × late grey on onset: Bayes ns, PLS negative (−0.20, p = 0.016). There is no evidence that coupling strengthens the prediction. No trial factor relates to coupling.

## 5. Theta from HPC MUA (`ar_theta_mua`, `ar_theta_phase`)
- **Spectrum.** While running, a sharp peak at ~7 Hz in every session. Prominence (6–9 Hz over flanks): 3.0 on track and 2.6 in grey (> 5 cm/s); the same at > 1 cm/s and with licks removed. Still: none (1.3; broadband synchrony).
- **By animal.** Prominence 4.0 / 2.3 / 2.7 / 1.8 / 4.4 / 3.8; M24062 is weakest.
- **Running around onset.** Running > 5 cm/s covers 57–99% of the first second of the track and 67–100% of the last second of grey.
- **Phase method.** Vollan et al. 2025: 10 ms bins, 5–10 Hz 2nd-order Butterworth, PCA, zero at minimal population firing. Cycles run between minimal-firing points, as in Jezek et al. 2011.
- **Phase results.**
  - The PC1–PC2 plane rotates at theta in 1/22 sessions only (spatial cells only, no interneurons), so the other sessions use the Hilbert phase of PC1.
  - Frequency 6.5–7.4 Hz; 76–100% of cells phase-locked.
  - Phase-difference consistency with the summed-MUA Hilbert phase: 0.62–0.97 (median 0.87).
  - Cycles ~140 ms, covering 82–97% of running time.
- **Pending.** LFP validation when the server is back.

## 6. Raw vs z log odds
- **Same conclusions for every key effect.**
- **Recommendation: z for pooling.** Raw values vary 10–30× in scale between sessions, and PLS raw tails dominate the PLS LMEs.

## 7. Does onset reset the state? (`ar_onset_reset`)
- **Test.** For each t0, correlate the residual bias over [t0 − 1, t0) with the next 200 ms. t0 = −0.5 to −0.2 s gives within-grey continuity; t0 = 0 gives continuity across onset.
- **HPC: no reset.** ρ is 0.40 across onset vs 0.42 within grey (difference −0.02, p = 0.43; PLS +0.00).
- **V1: partial reset.** ρ is 0.31 vs 0.39 (p = 0.013, 4/5 animals lower; PLS p = 0.049).
- **Reading.** The pre → post link is an ongoing population state with about a second of memory that carries through onset. It is not an expectation signal: no build-up over the block and no difference before switch laps.

## 7b. Which cells carry it (`ar_carryover_cells`, `results\carryover_cells\`)
- **Linear stand-in.** LO_lin = Σ w_i n_i with w_i = log rate ratio T1/T2. Session carry-over correlates with Bayes at r = 0.82 (HPC) and 0.76 (V1).
- **Decomposition** (contributions in correlation units; they sum to the linear carry-over r):

  | | Total | Own-cell persistence | Same-preference pairs | Opposite-preference pairs |
  |---|---|---|---|---|
  | HPC across onset | 0.48 | 0.26 (22/22 sessions) | 0.34 (21/22) | −0.11 (ns) |
  | HPC within grey | 0.49 | 0.32 | 0.26 | −0.09 |
  | V1 across onset | 0.34 | 0.10 (19/19) | 0.23 (18/19) | 0.01 |
  | V1 within grey | 0.45 | 0.17 | 0.37 | −0.09 |

  So it is half own persistence and half a track-structured co-fluctuation; it is not global gain, since opposite-preference pairs contribute ≈ 0.
- **Distributed.**
  - Single-cell carry-over is positive in every HPC field group (ρ 0.22–0.25, 22/22 sessions each).
  - Ablation with the Bayes decoder (group minus random cells): start, middle, end and high-rate are all ns.
  - In V1, high-rate cells carry +22% of the carry-over beyond their cell fraction (18/19 sessions).
- **HPC cell LME.** Carry-over increases with |track preference| (β 0.21, p = 1e-14) and firing rate (p = 1e-37). Field position is ns; theta MVL −0.08 (p = 0.049).

## 8. Theta-cycle flicker (`ar_theta_flicker`; MUA cycles, running > 5 cm/s, ≥ 3 spikes)
- **Data.** HPC cycles 1–5 after onset: 6,922 cycles, median 81 spikes; 16.5% decode as the other track.
- **No switch-specific flicker.** First lap of a block minus laps 6+ in other-track fraction: −0.02 (cycles 1–3, ns), −0.02 (4–6), −0.03 (7–10); 1–2/6 animals higher. V1 is the same.
- **Pre-onset cycles.** On first laps they are ~50/50 between the tracks, against ~60% toward the block's track on late laps.
- **After onset.** From cycle 1, every lap type expresses the new track equally; there is no gradual settling.
- **What predicts other-track cycles 1–5** (GLME):
  - late-grey bias toward the other track: HPC β 1.75 per z (p = 1e-81), V1 2.3 (p = 2e-22);
  - bias × switch: ns (HPC p = 0.95).
- **Conflict model** (HPC–V1 coupling × bias):
  - V1: own bias × coupling +2.0 (p = 1e-4), but opposite in sign to the 100 ms onset model.
  - HPC: own bias × coupling +0.41 (p = 0.07); coupling alone −0.39 (p = 4e-4).
  - Not consistent, so not claimed.
- **Examples.** Per session in `results\theta_flicker\examples\`.

## Next
- LFP validation of the MUA theta phase (server). Then flicker with half-cycles and the Jezek-style population-vector assignment (already stored per cycle).
- What the track-structured co-fluctuation in the grey is: sequence or replay-like events vs slow excitability of track assemblies (decode sequences in grey windows, check SWR / population bursts while running).
