import base64, os
D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, "figures")

def img(name):
    with open(os.path.join(R, name), "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

HEAD = """<title>Depth Gain and Smell Recovery</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Zilla+Slab:wght@400;500;600;700&family=Source+Sans+3:ital,wght@0,300;0,400;0,600;0,700;1,400&family=IBM+Plex+Mono:wght@400;500;600&display=swap">
<style>
  :root{--ground:#F6F8F7;--surface:#FFF;--surface-alt:#EDF1EF;--line:#D3DCD8;--line-soft:#E3E9E6;
    --ink:#151E1C;--ink-mid:#47554F;--ink-soft:#6D7C76;--resp:#16776A;--resp-soft:#DCEDE9;
    --card:#A32F47;--card-soft:#F6DFE3;--warn:#8A6410;--warn-soft:#F6EBD3;--ok:#2C6B3F;--ok-soft:#DDECE1;
    --stop:#93341F;--stop-soft:#F5E1DA;}
  @media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
    --ground:#0E1513;--surface:#161F1D;--surface-alt:#1C2725;--line:#2C3937;--line-soft:#232F2D;
    --ink:#E3EAE7;--ink-mid:#A9B7B2;--ink-soft:#7E8D88;--resp:#4FBBA6;--resp-soft:#12332E;
    --card:#E4788C;--card-soft:#391A22;--warn:#D9A748;--warn-soft:#33280F;--ok:#6FBE86;--ok-soft:#16301F;
    --stop:#E08767;--stop-soft:#351D14;}}
  :root[data-theme="dark"]{
    --ground:#0E1513;--surface:#161F1D;--surface-alt:#1C2725;--line:#2C3937;--line-soft:#232F2D;
    --ink:#E3EAE7;--ink-mid:#A9B7B2;--ink-soft:#7E8D88;--resp:#4FBBA6;--resp-soft:#12332E;
    --card:#E4788C;--card-soft:#391A22;--warn:#D9A748;--warn-soft:#33280F;--ok:#6FBE86;--ok-soft:#16301F;
    --stop:#E08767;--stop-soft:#351D14;}
  *{box-sizing:border-box}
  body{background:var(--ground);color:var(--ink);font-family:"Source Sans 3",ui-sans-serif,system-ui,sans-serif;
    font-size:17px;line-height:1.62;margin:0;padding:0 clamp(16px,5vw,40px) 96px;-webkit-font-smoothing:antialiased}
  .wrap{max-width:1120px;margin:0 auto}.col{max-width:68ch}
  h1,h2,h3{font-family:"Zilla Slab",Georgia,serif;text-wrap:balance;margin:0;line-height:1.18}
  .mast{padding:clamp(40px,7vw,72px) 0 30px;border-bottom:1px solid var(--line)}
  .eyebrow{font-family:"IBM Plex Mono",monospace;font-size:11.5px;font-weight:500;letter-spacing:.16em;
    text-transform:uppercase;color:var(--ink-soft);margin:0 0 18px}
  h1{font-size:clamp(32px,5.5vw,52px);font-weight:600;letter-spacing:-.015em}
  h1 .a{color:var(--card)}
  .standfirst{font-size:clamp(17px,2vw,20px);color:var(--ink-mid);max-width:62ch;margin:20px 0 0}
  section{padding-top:56px}
  h2{font-size:clamp(23px,3.1vw,30px);font-weight:600;letter-spacing:-.01em;padding-bottom:12px;
    border-bottom:2px solid var(--ink);margin-bottom:24px}
  h3{font-size:18.5px;font-weight:600;margin:30px 0 10px}
  p{margin:0 0 17px}
  code{font-family:"IBM Plex Mono",monospace;font-size:.875em;background:var(--surface-alt);
    border:1px solid var(--line-soft);border-radius:3px;padding:.08em .34em;white-space:nowrap}
  .mono{font-family:"IBM Plex Mono",monospace;font-size:.9em}
  ul{margin:0 0 17px;padding-left:1.15em}li{margin-bottom:8px}
  .stats{display:grid;gap:1px;grid-template-columns:repeat(auto-fit,minmax(185px,1fr));background:var(--line);
    border:1px solid var(--line);border-radius:3px;overflow:hidden;margin:36px 0 0}
  .stat{background:var(--surface);padding:18px 20px 20px}
  .stat .k{font-family:"IBM Plex Mono",monospace;font-size:10.5px;letter-spacing:.13em;text-transform:uppercase;
    color:var(--ink-soft);margin:0 0 9px}
  .stat .v{font-family:"Zilla Slab",Georgia,serif;font-size:37px;font-weight:600;line-height:1;font-variant-numeric:tabular-nums}
  .stat .v.pos{color:var(--resp)}.stat .v.neu{color:var(--card)}
  .stat .n{font-size:13.5px;color:var(--ink-mid);margin:8px 0 0;line-height:1.45}
  figure{margin:26px 0 30px;background:var(--surface);border:1px solid var(--line);border-radius:4px;padding:16px}
  figure img{width:100%;height:auto;display:block;border-radius:2px}
  figcaption{font-size:14px;color:var(--ink-mid);margin-top:12px;line-height:1.5}
  figcaption b{color:var(--ink)}
  .note{background:var(--surface);border:1px solid var(--line);border-left:3px solid var(--resp);
    border-radius:3px;padding:18px 22px;margin:0 0 22px}
  .note.warn{border-left-color:var(--warn)}.note.stop{border-left-color:var(--card)}
  .note p:last-child{margin-bottom:0}
  .note .lab{font-family:"IBM Plex Mono",monospace;font-size:10.5px;letter-spacing:.13em;text-transform:uppercase;
    color:var(--ink-soft);display:block;margin-bottom:8px}
  .scroll{overflow-x:auto;border:1px solid var(--line);border-radius:3px;background:var(--surface);margin-bottom:22px}
  table{border-collapse:collapse;width:100%;font-size:14.5px}
  th,td{text-align:left;padding:8px 13px;border-bottom:1px solid var(--line-soft);white-space:nowrap}
  thead th{font-family:"IBM Plex Mono",monospace;font-size:10px;font-weight:600;letter-spacing:.11em;
    text-transform:uppercase;color:var(--ink-soft);background:var(--surface-alt);border-bottom:1px solid var(--line)}
  tbody tr:last-child td{border-bottom:0}
  td.num{font-variant-numeric:tabular-nums;text-align:right}
  td.sub{font-family:"IBM Plex Mono",monospace;font-weight:700;font-size:13px}
  td.w{white-space:normal;min-width:280px;font-size:14px;line-height:1.5}
  tr.hi td{background:var(--ok-soft)}
  tr.dim td{color:var(--ink-soft)}
  .pill{display:inline-block;font-family:"IBM Plex Mono",monospace;font-size:10px;font-weight:600;
    letter-spacing:.06em;padding:2px 7px;border-radius:2px;text-transform:uppercase}
  .p-keep{background:var(--ok-soft);color:var(--ok)}.p-drop{background:var(--stop-soft);color:var(--stop)}
  .p-warn{background:var(--warn-soft);color:var(--warn)}
  caption{caption-side:bottom;text-align:left;font-size:13px;color:var(--ink-soft);padding:10px 13px;
    border-top:1px solid var(--line-soft)}
  footer{margin-top:60px;padding-top:20px;border-top:1px solid var(--line);font-size:14px;color:var(--ink-soft)}
</style>"""

BODY = f"""
<div class="wrap">
<header class="mast">
  <p class="eyebrow">Dupilumab cohort · respiration → HRV · 2 September 2026</p>
  <h1>Depth Gain <span class="a">&amp;</span> Smell Recovery</h1>
  <p class="standfirst">Patients whose sense of smell improved most over the first month of treatment
  also showed the largest increase in how strongly a deeper breath drives their heart rate. The effect
  is specific to respiratory–cardiac coupling, absent from general heart-rate variability, and does not
  survive multiplicity correction at this sample size.</p>
  <div class="stats">
    <div class="stat"><p class="k">Primary correlation</p><div class="v pos">+0.79</div>
      <p class="n">Δ depth gain vs Δ olfaction, session 1→2, n = 7 subjects</p></div>
    <div class="stat"><p class="k">Family-corrected</p><div class="v neu">0.159</div>
      <p class="n">Max-statistic permutation across 7 metrics — not significant</p></div>
    <div class="stat"><p class="k">Sessions analysed</p><div class="v">25</div>
      <p class="n">Of 29 breathing finals; 13 change intervals from 9 subjects</p></div>
    <div class="stat"><p class="k">Duration-slope null</p><div class="v">+0.69</div>
      <p class="n">Mean artifact in the original estimand. Depth-slope null: ≈ 0</p></div>
  </div>
</header>

<section>
  <div class="col">
    <h2>What was asked, and what changed</h2>
    <p>The question was whether the respiration → HRV transfer function steepens across sessions, and
    whether that steepening tracks olfactory recovery. The honest answer required first replacing the
    estimand, because the original one could not have answered it.</p>
  </div>
  <div class="note stop">
    <span class="lab">The original slope was measuring an artifact</span>
    <p>Regressing within-breath heart-period range on breath <b>duration</b> has a built-in positive
    slope that has nothing to do with physiology: a longer breath is a longer window, a longer window
    contains more heartbeats, and the range of more samples is simply larger. Across the 25 sessions the
    surrogate null for the duration slope ran <span class="mono">+0.44 to +1.03</span> — the same size as
    the raw effect being measured.</p>
  </div>
  <div class="note">
    <span class="lab">Breath depth has no such artifact</span>
    <p>Depth does not set the length of the analysis window, so the order-statistics effect does not
    apply. Measured directly: the depth-slope null across all 25 sessions falls between
    <span class="mono">−0.011 and +0.022</span> — indistinguishable from zero. The depth slope
    (θ<sub>depth</sub>) is therefore the well-posed version of the question: <b>how much cardiac
    modulation does a deeper breath recruit?</b></p>
  </div>
  <figure><img src="{img('fig1_nulls.png')}" alt="Left: bar chart of surrogate null slopes for duration shift, duration i.i.d., and depth. Duration nulls near 0.69, depth null near zero. Right: scatter of observed depth slope against its null, all points above the identity line.">
    <figcaption><b>The artifact, measured.</b> Both duration nulls sit near +0.69 while the depth null
    sits at zero. On the right, every session's observed depth slope exceeds its own null — real signal
    against a genuinely flat baseline.</figcaption></figure>
</section>

<section>
  <div class="col">
    <h2>The primary result</h2>
    <p>The pre-specified contrast is session 1 → 2, the interval every subject shares and which is
    tightly matched at 28–42 days. Seven subjects have it: AB, DB, GH, JH, JL, KS and PC.</p>
  </div>
  <figure><img src="{img('fig2_primary.png')}" alt="Scatter of change in depth slope against change in olfactory composite for seven subjects with error bars, showing a positive trend from PC at lower left to JH at upper right.">
    <figcaption><b>Δ depth gain vs Δ olfaction, session 1→2.</b> Spearman ρ = +0.79. PC, the one subject
    whose olfaction declined, is also the one whose depth gain fell most. JH, who improved most, gained
    most. Error bars are surrogate-derived standard errors on each subject's change.</figcaption></figure>
  <div class="scroll"><table>
    <thead><tr><th>Measure</th><th>What it is</th><th>ρ</th><th>Jackknife range</th><th>P(ρ&gt;0)</th></tr></thead>
    <tbody>
      <tr class="hi"><td class="sub">adjLogRSA</td><td class="w">RSA level at a fixed reference breath</td><td class="num">+0.886</td><td class="num">+0.80 to +0.90</td><td class="num">0.992</td></tr>
      <tr class="hi"><td class="sub">thetaVol</td><td class="w">Depth slope, surrogate-corrected</td><td class="num">+0.786</td><td class="num">+0.66 to +1.00</td><td class="num">1.000</td></tr>
      <tr class="hi"><td class="sub">logMedPV</td><td class="w">Median peak-valley RSA per breath</td><td class="num">+0.771</td><td class="num">+0.60 to +0.90</td><td class="num">0.999</td></tr>
      <tr class="hi"><td class="sub">b_vol</td><td class="w">Depth slope, uncorrected</td><td class="num">+0.750</td><td class="num">+0.60 to +0.94</td><td class="num">0.999</td></tr>
      <tr><td class="sub">logRMSSD</td><td class="w">Beat-to-beat variability</td><td class="num">+0.486</td><td class="num">+0.10 to +0.80</td><td class="num">0.989</td></tr>
      <tr class="dim"><td class="sub">logMASD</td><td class="w">Outlier-robust beat-to-beat variability</td><td class="num">−0.257</td><td class="num">−0.60 to +0.10</td><td class="num">0.203</td></tr>
      <tr class="dim"><td class="sub">pNN20</td><td class="w">Proportion of successive differences &gt; 20 ms</td><td class="num">−0.257</td><td class="num">−0.60 to +0.10</td><td class="num">0.123</td></tr>
    </tbody>
    <caption>P(ρ&gt;0) is the fraction of 20,000 draws in which the rank correlation stays positive when
    each subject's change score is perturbed by its own standard error. n = 7 for the depth slopes,
    n = 6 for panel metrics (JL failed the beat-quality gate).</caption>
  </table></div>
</section>

<section>
  <div class="col">
    <h2>The pattern is construct-specific</h2>
    <p>This is the part that makes the result more than a lucky ordering. The metrics do not scatter
    randomly around zero: they split cleanly along a physiological line.</p>
  </div>
  <figure><img src="{img('fig3_specificity.png')}" alt="Horizontal bar chart of Spearman correlations by metric. Four RSA and depth-response measures are positive between 0.75 and 0.89 in teal; three beat-to-beat variability measures are near zero or negative in red.">
    <figcaption><b>Four coupling measures positive, three general-HRV measures not.</b> Everything that
    indexes respiratory–cardiac coupling — RSA level, depth response, peak-valley RSA — lands between
    +0.75 and +0.89. Everything that indexes generic beat-to-beat variability lands between −0.26 and
    +0.49.</figcaption></figure>
  <div class="col">
    <p>If this were noise, the seven metrics would not partition by construct. The interpretation it
    supports is narrow and specific: what tracks olfactory recovery is <b>how strongly respiration drives
    the heart</b>, not how variable the heart is overall. That is a coupling story, not a vagal-tone
    story, and it is a more constrained claim than the one we set out to test.</p>
    <p>Two confound checks came back clean. Change in mean heart rate correlates with olfactory change at
    ρ = −0.21, and change in breathing rate at ρ = +0.09 — both far below the effect, so neither the
    heart-rate mechanism nor a breathing-rate shift explains it.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>What does not hold</h2>
    <p>Three things temper this, and they belong in the same breath as the result.</p>
    <h3>It is not significant after multiplicity correction</h3>
    <p>Screening seven metrics against one outcome inflates the best one. The honest inferential device
    is a max-statistic permutation that shuffles subject labels and records the largest |ρ| anywhere in
    the family. Observed max is 0.886; the permutation p is <b>0.159</b>. With 5,040 orderings the
    resolution floor is 0.0002, so this is a real failure to reach significance, not a resolution limit.</p>
    <h3>It weakens when the six-month intervals are pooled in</h3>
    <p>Across all 13 adjacent intervals the depth-slope correlation attenuates to +0.335 and the
    beat-to-beat measures turn negative. The max statistic falls to 0.335 with a permutation
    <b>p = 0.729</b> over 362,880 orderings — no signal at all once the intervals are pooled. The
    heart-rate confound also grows, from ρ = −0.21 in the primary contrast to <b>+0.39</b> here, so the
    pooled analysis is both weaker and less clean.</p>
    <p>The session 2 → 3 intervals span ~140 days rather than ~31, and four of six subjects <em>lose</em>
    olfactory ground over them. Whether that reflects waning response, variable adherence, or something
    procedural is a question for you — session means show no systematic artifact (composite −0.51, +0.28,
    +0.26 for sessions 1/2/3), so it looks like genuine heterogeneity in trajectory rather than a
    measurement problem.</p>
  </div>
  <figure><img src="{img('fig4_traj.png')}" alt="Two line plots across sessions 1 to 3. Left: olfactory composite trajectories per subject, some rising throughout, others peaking at session 2. Right: depth-gain trajectories per subject.">
    <figcaption><b>Trajectories diverge after the first month.</b> JH and AB keep improving through six
    months; GH, JN, KS and TB peak at one month and fall back. Pooling a 31-day change with a 140-day
    change assumes a linearity these trajectories do not support.</figcaption></figure>
  <div class="col">
    <h3>The sample is still small</h3>
    <p>Seven subjects for the primary contrast. A perfect rank ordering at n = 7 reaches p ≈ 0.024
    one-sided for a single pre-declared metric — but log RMSSD was the pre-declared primary, and it came
    in at +0.486. The depth slope, which produced the headline, was pre-specified as a better-posed
    alternative but not as <em>the</em> primary, so claiming its p-value uncorrected would be
    retrospective.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>Reliability gate</h2>
    <p>Before any outcome was touched, every candidate metric was screened on κ = mean SE² / var(Δ),
    which asks whether the between-subject spread exceeds the measurement noise. Metrics failing κ ≤ 0.5
    were dropped without ever being correlated with olfaction. This is what the previous round of this
    analysis failed outright, when the duration slope came in at κ ≈ 1.</p>
  </div>
  <div class="scroll"><table>
    <thead><tr><th>Metric</th><th>var(Δ)</th><th>mean SE</th><th>κ</th><th>Verdict</th></tr></thead>
    <tbody>
      <tr><td class="sub">pNN20</td><td class="num">108.63</td><td class="num">3.202</td><td class="num">0.10</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr><td class="sub">logMASD</td><td class="num">0.0706</td><td class="num">0.092</td><td class="num">0.13</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr><td class="sub">logMedPV</td><td class="num">0.0616</td><td class="num">0.105</td><td class="num">0.19</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr><td class="sub">logRMSSD</td><td class="num">0.0468</td><td class="num">0.096</td><td class="num">0.22</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr class="hi"><td class="sub">thetaVol</td><td class="num">0.0201</td><td class="num">0.065</td><td class="num">0.23</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr class="hi"><td class="sub">b_vol</td><td class="num">0.0179</td><td class="num">0.065</td><td class="num">0.26</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr class="hi"><td class="sub">adjLogRSA</td><td class="num">0.0440</td><td class="num">0.113</td><td class="num">0.30</td><td><span class="pill p-keep">retain</span></td></tr>
      <tr class="dim"><td class="sub">b_len</td><td class="num">0.3284</td><td class="num">0.423</td><td class="num">0.60</td><td><span class="pill p-drop">drop</span></td></tr>
      <tr class="dim"><td class="sub">theta (duration)</td><td class="num">0.1837</td><td class="num">0.423</td><td class="num">1.08</td><td><span class="pill p-drop">drop</span></td></tr>
      <tr class="dim"><td class="sub">depthSlope (panel)</td><td class="num">0.0862</td><td class="num">0.305</td><td class="num">1.14</td><td><span class="pill p-drop">drop</span></td></tr>
      <tr class="dim"><td class="sub">logHF</td><td class="num">0.0692</td><td class="num">0.398</td><td class="num">2.46</td><td><span class="pill p-drop">drop</span></td></tr>
      <tr class="dim"><td class="sub">evokedRSA</td><td class="num">0.0000</td><td class="num">0.009</td><td class="num">3.09</td><td><span class="pill p-drop">drop</span></td></tr>
    </tbody>
    <caption>The duration slope fails its own reliability gate. The depth slope passes comfortably.</caption>
  </table></div>
</section>

<section>
  <div class="col">
    <h2>Data problems found</h2>
  </div>
  <div class="note warn">
    <span class="lab">Metadata gap — cost 11 sessions on the first pass</span>
    <p>The <code>noseMouth</code> field is empty for 31–46% of breaths in 11 sessions, concentrated in the
    <code>audio</code> and <code>focus</code> blocks. No session anywhere carries an explicit
    <code>"mouth"</code>; values are <code>"nose"</code>, <code>"NA"</code> or empty. Since respiration is
    recorded with a nasal cannula, a detected breath with a valid RSA value was necessarily nasal, so
    empty means unrecorded metadata. Filtering on <code>== "nose"</code> silently discarded most of the
    newly preprocessed data and collapsed those sessions to a single block. Worth populating the field,
    or documenting that empty means nose.</p>
  </div>
  <div class="note warn">
    <span class="lab">Corrupt file</span>
    <p><span class="mono">250811_Dupi_NMH_TB_2_PEA_threshold_preproc.mat</span> fails to open with an HDF5
    <code>inflate()</code> error. It has now failed on two attempts four days apart, so this is genuine
    corruption rather than a network read failure. TB_2 has no threshold component until it is rebuilt.</p>
  </div>
  <div class="note warn">
    <span class="lab">Beat quality</span>
    <p>JL_1 rejects 42% of its NN intervals and fails the beat-quality gate. PD_1, JA_1 and JA_2 yield only
    52–59 spontaneous breaths against a 60-breath minimum, so they contribute to the slope analysis but not
    the vagal panel. An ectopic filter was added this round after the previous one — which only checked
    each interval against a local median — let through short-long ectopic pairs that inflated RMSSD
    six-fold in one session and manufactured a spurious perfect correlation.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>What I would do next</h2>
    <ul>
      <li><b>Pre-register the depth slope as the primary metric and re-test.</b> It is now well motivated
      on artifact grounds, it passes the reliability gate at κ = 0.23, and it produced the headline. A
      single pre-declared metric at n = 10 would reach significance for ρ ≈ 0.65, which is below what was
      observed here.</li>
      <li><b>Add subjects with a session 1→2 pair.</b> PD, JA and BS each need one more usable session to
      enter the primary contrast; PD_1, JA_1 and JA_2 are three of the four sessions currently just below
      the breath-count threshold, so a longer spontaneous block in the protocol would recover them.</li>
      <li><b>Decide the six-month question separately.</b> Whether coupling tracks olfaction over 140 days
      is a different question from whether it does over 31, and the trajectories say they should not be
      pooled.</li>
      <li><b>Add a 5-minute seated spontaneous-breathing block</b> to every session. All the vagal metrics
      rest on ~9 minutes of <code>audio</code> + <code>focus</code> data; doubling that is the cheapest
      available reduction in κ, and κ is what governs whether any n helps.</li>
    </ul>
  </div>
</section>

<footer>
  25 of 29 breathing finals · 13 adjacent intervals · 9 subjects · primary contrast n = 7.
  Beat series from hand-tuned per-session <span class="mono">beatSpec</span>; spontaneous blocks only for
  vagal metrics; block-centred within session; 200 circular-shift and 200 permuted-NN surrogates per
  session; 300-draw block bootstrap for standard errors; 20,000-draw perturbation for rank stability;
  exact subject-level permutation for inference. Analysis run 2 September 2026 against preprocessed
  files stable since 1 September 23:14.
</footer>
</div>
"""

os.makedirs(os.path.join(D, "reports"), exist_ok=True)
with open(os.path.join(D, "reports", "rsa_report.html"), "w", encoding="utf-8") as f:
    f.write(HEAD + BODY)
print("wrote reports/rsa_report.html", os.path.getsize(os.path.join(D, "reports", "rsa_report.html")) // 1024, "KB")
