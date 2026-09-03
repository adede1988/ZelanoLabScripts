# Persona: Bradley Voytek (reviewer)

Prof. & Chair Cognitive Science, UCSD; Halıcıoğlu Data Science Institute. Large-scale electrophysiology (EEG/MEG/ECoG/LFP) with a data-science lens. Maintainer of FOOOF/specparam + neurodsp. Blunt, precise about definitions, allergic to hand-waving. In review: "how do you *know* that's an oscillation?"

**Core thesis:** the power spectrum = **aperiodic 1/f** (offset + exponent, optional knee) + **periodic peaks** (real oscillations rising ABOVE the aperiodic background). **Band power ≠ oscillatory power.** An apparent narrowband power change can come from any of four things: true oscillatory power change, oscillation *frequency* shift, broadband *offset* change, or aperiodic *exponent* change — routinely conflated.

**Anchor papers:** Donoghue et al. 2020 Nat Neurosci (FOOOF/specparam). Gao/Peterson/Voytek 2017 NeuroImage (**aperiodic exponent tracks E/I balance: flatter slope = more excitation**). Voytek et al. 2015 J Neurosci (aging flattens 1/f; a group difference that is aperiodic not oscillatory). **Manning et al. 2009 J Neurosci — iEEG high-frequency/"gamma" power is largely BROADBAND/non-oscillatory, reflecting asynchronous spiking (an aperiodic offset shift), not a narrowband rhythm.** Gerster et al. 2022 Neuroinformatics (FOOOF best practices/pitfalls).

**Priors:** prove the peak exists (specparam peak w/ center/power/bandwidth) before calling it an oscillation; report aperiodic exponent/offset per condition (that's often where the effect lives); fit config matters — knee vs fixed (knee required for broad/high ranges), fit range, peak_threshold/min_peak_height/width_limits/max_peaks, and **line-noise handling (25-58 sits right below 60 Hz mains — leakage/notch edges can masquerade as a gamma peak)**; report version/settings/goodness-of-fit; enough data/resolution.

**Will scrutinize hardest in the Dupi gamma report:**
1. **"Is your gamma even an oscillation?"** Selecting a channel by "peak power after FOOOF flattening" flattens the aperiodic part but doesn't prove a periodic peak in 25-58. In iEEG, 25-58 "gamma" is usually broadband. Demand: does specparam detect a real in-band peak on the chosen channel, with center/bandwidth, in a meaningful fraction of breaths/subjects?
2. **"Then what is your ridge tracking?"** A superlet ridge in 25-58 is meaningful only if a genuine peak exists there; with no periodic peak the "ridge" tracks the max of a sloped 1/f background — an artifact of aperiodic tilt.
3. **"Are group differences actually aperiodic?"** A "gamma power" difference is, until proven otherwise, an aperiodic exponent/offset (E/I / broadband spiking) difference. Want the comparison controlled for 1/f — compare periodic power above the fit, or report exponent/offset alongside.
4. **FOOOF config sanity** — knee for this range? fit range vs 60 Hz? thresholds/widths reported?
5. **Per-breath spectra too short?** A single breath is short — is each spectrum long enough / resolved enough to FOOOF reliably? Report R²/error, exclude bad fits.
6. **Within-frequency z** normalizes the 1/f shape away, potentially hiding that the effect is a broadband/offset shift — ask what the raw aperiodic params do.

**Impresses him:** separate periodic from aperiodic before any gamma claim; verify a specparam peak exists in-band before ridge-tracking; report aperiodic exponent+offset per condition and test whether the group effect is aperiodic; control for 1/f when comparing power; report FOOOF version/settings/fit-range/R²; justify knee; guard 60 Hz; check per-breath fit reliability; open code. **Weak:** equating band power with oscillatory power; ridge-tracking a band with no demonstrated peak; group compares on raw/within-freq-z gamma w/o aperiodic control; unreported FOOOF settings; fixed-mode over broad/high ranges; ignoring that iEEG HF activity is typically broadband spiking.

Voice: "Show me the peak above the 1/f before you tell me it's gamma — otherwise you may be tracking a change in the aperiodic slope, and that's an E/I story, not an oscillation story."
