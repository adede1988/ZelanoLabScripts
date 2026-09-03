import base64, os
D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, "figures")
def img(n):
    with open(os.path.join(R, n), "rb") as f:
        return "data:image/png;base64," + base64.b64encode(f.read()).decode()

HEAD = """<title>Coupled Recovery</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Zilla+Slab:wght@400;500;600;700&family=Source+Sans+3:ital,wght@0,300;0,400;0,600;0,700;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
 :root{--ground:#F6F8F7;--surface:#FFF;--surface-alt:#EDF1EF;--line:#D3DCD8;--line-soft:#E3E9E6;
  --ink:#151E1C;--ink-mid:#47554F;--ink-soft:#6D7C76;--resp:#16776A;--card:#A32F47;
  --ok-soft:#DDECE1;--ok:#2C6B3F;--warn-soft:#F6EBD3;--warn:#8A6410;}
 @media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --ground:#0E1513;--surface:#161F1D;--surface-alt:#1C2725;--line:#2C3937;--line-soft:#232F2D;
  --ink:#E3EAE7;--ink-mid:#A9B7B2;--ink-soft:#7E8D88;--resp:#4FBBA6;--card:#E4788C;
  --ok-soft:#16301F;--ok:#6FBE86;--warn-soft:#33280F;--warn:#D9A748;}}
 :root[data-theme="dark"]{
  --ground:#0E1513;--surface:#161F1D;--surface-alt:#1C2725;--line:#2C3937;--line-soft:#232F2D;
  --ink:#E3EAE7;--ink-mid:#A9B7B2;--ink-soft:#7E8D88;--resp:#4FBBA6;--card:#E4788C;
  --ok-soft:#16301F;--ok:#6FBE86;--warn-soft:#33280F;--warn:#D9A748;}
 *{box-sizing:border-box}
 body{background:var(--ground);color:var(--ink);font-family:"Source Sans 3",ui-sans-serif,system-ui,sans-serif;
  font-size:17.5px;line-height:1.65;margin:0;padding:0 clamp(16px,5vw,40px) 90px;-webkit-font-smoothing:antialiased}
 .wrap{max-width:1080px;margin:0 auto}.col{max-width:70ch}
 h1,h2,h3{font-family:"Zilla Slab",Georgia,serif;text-wrap:balance;margin:0;line-height:1.2}
 .mast{padding:clamp(38px,7vw,70px) 0 28px;border-bottom:1px solid var(--line)}
 .eyebrow{font-family:"IBM Plex Mono",monospace;font-size:11.5px;letter-spacing:.16em;text-transform:uppercase;
  color:var(--ink-soft);margin:0 0 16px}
 h1{font-size:clamp(32px,5.4vw,50px);font-weight:600;letter-spacing:-.015em}
 .standfirst{font-size:clamp(17.5px,2vw,21px);color:var(--ink-mid);max-width:64ch;margin:20px 0 0}
 section{padding-top:52px}
 h2{font-size:clamp(23px,3vw,29px);font-weight:600;padding-bottom:11px;border-bottom:2px solid var(--ink);margin-bottom:22px}
 h3{font-size:19px;font-weight:600;margin:28px 0 9px;color:var(--resp)}
 p{margin:0 0 17px}
 .lede{font-size:19px;color:var(--ink-mid)}
 figure{margin:24px 0 28px;background:var(--surface);border:1px solid var(--line);border-radius:4px;padding:15px}
 figure img{width:100%;height:auto;display:block;border-radius:2px}
 figcaption{font-size:14.5px;color:var(--ink-mid);margin-top:11px;line-height:1.5}
 figcaption b{color:var(--ink);font-weight:600}
 .fignum{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;text-transform:uppercase;
  color:var(--resp);font-weight:500}
 .key{background:var(--surface);border:1px solid var(--line);border-left:3px solid var(--resp);
  border-radius:3px;padding:20px 24px;margin:0 0 24px}
 .key p:last-child{margin-bottom:0}
 .key .lab{font-family:"IBM Plex Mono",monospace;font-size:11px;letter-spacing:.14em;text-transform:uppercase;
  color:var(--ink-soft);display:block;margin-bottom:9px}
 .stats{display:grid;gap:1px;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));background:var(--line);
  border:1px solid var(--line);border-radius:3px;overflow:hidden;margin:34px 0 0}
 .stat{background:var(--surface);padding:18px 20px}
 .stat .k{font-family:"IBM Plex Mono",monospace;font-size:10.5px;letter-spacing:.13em;text-transform:uppercase;
  color:var(--ink-soft);margin:0 0 8px}
 .stat .v{font-family:"Zilla Slab",Georgia,serif;font-size:38px;font-weight:600;line-height:1;color:var(--resp);
  font-variant-numeric:tabular-nums}
 .stat .n{font-size:13.5px;color:var(--ink-mid);margin:8px 0 0;line-height:1.45}
 ul{margin:0 0 17px;padding-left:1.2em}li{margin-bottom:10px}
 code{font-family:"IBM Plex Mono",monospace;font-size:.87em;background:var(--surface-alt);
  border:1px solid var(--line-soft);border-radius:3px;padding:.08em .34em}
 .todo{background:var(--warn-soft);color:var(--warn);font-family:"IBM Plex Mono",monospace;font-size:12px;
  padding:2px 8px;border-radius:2px;font-weight:500}
 footer{margin-top:56px;padding-top:20px;border-top:1px solid var(--line);font-size:14px;color:var(--ink-soft)}
</style>"""

BODY = f"""
<div class="wrap">
<header class="mast">
  <p class="eyebrow">Preliminary data · dupilumab cohort · respiration–cardiac coupling</p>
  <h1>Coupled Recovery</h1>
  <p class="standfirst">In patients regaining their sense of smell on dupilumab, the heart's response to
  each breath strengthens in proportion to how much olfactory function returns. Smell recovery and
  respiratory–cardiac coupling appear to recover together.</p>
  <div class="stats">
    <div class="stat"><p class="k">Patients, pre/post</p><div class="v">7</div>
      <p class="n">Baseline and ~1 month, full olfactory battery plus continuous respiration and ECG</p></div>
    <div class="stat"><p class="k">Showed smell gains</p><div class="v">6 / 7</div>
      <p class="n">On a composite of discrimination, identification and threshold</p></div>
    <div class="stat"><p class="k">Coupling vs olfaction</p><div class="v">r = .85</div>
      <p class="n">Change in breath-depth drive of HRV against change in olfactory function</p></div>
    <div class="stat"><p class="k">Converging measures</p><div class="v">4</div>
      <p class="n">Independent formulations of the coupling all point the same way</p></div>
  </div>
</header>

<section>
  <div class="col">
    <h2>Rationale</h2>
    <p class="lede">Nasal airflow is not only the vehicle for odorants — it is itself a rhythmic afferent
    signal that entrains central and autonomic activity. The prediction that follows is direct: if a
    treatment restores the nasal sensory pathway, the breath-locked drive on downstream physiology should
    strengthen alongside it.</p>
    <p>Dupilumab offers an unusually clean test. Every treated patient experiences polyp reduction and
    improved nasal patency, but olfactory recovery varies widely and is not predicted by the airway
    change. Airway mechanics are therefore roughly held constant across patients while sensory
    restoration varies — so a physiological measure that tracks <em>olfactory</em> outcome specifically,
    rather than tracking everyone equally, is informative about the sensory pathway rather than the
    plumbing.</p>
    <p><span class="todo">CITATIONS TO ADD</span> — nasal-airflow entrainment of limbic and cortical
    oscillations; respiratory sinus arrhythmia and pulmonary afferent drive; dupilumab olfactory outcomes
    in CRSwNP.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>The premise holds: recovery is large and highly variable</h2>
    <p>Patients enter at or near the floor and diverge sharply. That variability is the whole basis for
    the analysis — it is what allows a physiological measure to be tested against outcome rather than
    against time on drug.</p>
  </div>
  <figure>
    <img src="{img('grantA_olf.png')}" alt="Three line plots showing per-patient trajectories across sessions for odor discrimination, identification, and a composite score. Most patients start near zero and rise, with wide spread by the second session.">
    <figcaption><span class="fignum">Figure 1</span> · <b>Olfactory recovery across sessions.</b>
    Most patients begin at or near floor on discrimination and identification. By one month the spread is
    substantial, and by six months trajectories diverge further — some continuing to improve, others
    falling back. Composite is the mean of z-scored discrimination, identification and air-calibrated
    threshold.</figcaption>
  </figure>
</section>

<section>
  <div class="col">
    <h2>Cardiac coupling recovers in step with it</h2>
    <p>For each session we measured, breath by breath, how much the heart period swings within a breath,
    and related it to the size of that breath. Two quantities follow. The first is how much cardiac
    modulation a breath of <em>fixed</em> size produces — a level. The second is how much <em>extra</em>
    modulation a deeper breath produces — a gain. Both track olfactory change.</p>
  </div>
  <figure>
    <img src="{img('grantB_main.png')}" alt="Two scatter plots. Left: change in RSA at a standard breath against change in olfactory function, rising, n=6. Right: change in cardiac response to breath depth against olfactory change, rising, n=7.">
    <figcaption><span class="fignum">Figure 2</span> · <b>Both formulations rise with olfactory recovery.</b>
    Left: with breath size held fixed, patients who recovered more smell produce more cardiac modulation
    per breath (ρ = 0.89). Right: they also show a steeper cardiac response to breath depth (r = 0.85).
    PC, the one patient whose olfactory score declined, sits at the bottom of both. KS is the one patient
    off the line on the left, which is why the left panel's rank correlation exceeds its linear one.
    Session 1 → 2, ~1 month.</figcaption>
  </figure>
  <div class="key">
    <span class="lab">Why this is not simply "more HRV"</span>
    <p>General measures of heart-rate variability — beat-to-beat interval variability and its robust and
    threshold-based variants — do <em>not</em> track olfactory recovery in this cohort. What tracks it is
    specifically the coupling between breathing and the heart. The effect is also not explained by
    changes in heart rate or in breathing rate, both of which are near zero against olfactory change.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>A mechanistic hint: the response shifts from time to volume</h2>
    <p>Breath duration and breath depth move in opposite directions. In patients who recovered more
    smell, the cardiac response became <em>more</em> sensitive to how much air moved and <em>less</em>
    sensitive to how long the breath lasted.</p>
  </div>
  <figure>
    <img src="{img('grantC_conv.png')}" alt="Left: horizontal bars showing four coupling measures all correlating positively with olfactory change. Right: scatter with depth-slope changes rising and duration-slope changes falling against olfactory change.">
    <figcaption><span class="fignum">Figure 3</span> · <b>Convergence, and a dissociation.</b>
    Left: four independent formulations of respiratory–cardiac coupling agree in direction and magnitude.
    Right: depth sensitivity rises with recovery while duration sensitivity falls — the cardiac response
    reallocates from time-driven to volume-driven.</figcaption>
  </figure>
  <div class="col">
    <p>This is what tighter breath-locking looks like. If the cardiac oscillation becomes firmly entrained
    to the breath cycle, each breath carries one rise-and-fall regardless of its duration, so duration
    stops buying extra modulation — while the amplitude of that single locked cycle continues to scale
    with respiratory drive. One change in entrainment produces both halves of the pattern.</p>
    <p>It is also what one would predict from the afferent side. Pulmonary stretch receptors and nasal
    flow-sensitive afferents signal <em>volume</em>, not elapsed time. A stronger, more reliable per-breath
    afferent signal should make the autonomic response more volume-graded and more tightly phase-locked —
    which is precisely the shift observed.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>The effect appears against the most peripheral measure</h2>
    <p>The threshold task yields a bias-calibrated index of intensity sensitivity. Each odor
    concentration is rated against a blank-air trial in the same session, so subtracting air removes
    both the rating offset and any individual response bias, leaving (high − air) and (med − air) as
    sensitivity contrasts on a common footing. Of the two, the high concentration carries the usable
    dynamic range; the medium sits near floor for most patients.</p>
  </div>
  <figure>
    <img src="{img('grantD_thresh.png')}" alt="Left: per-patient trajectories of bias-calibrated intensity sensitivity across sessions, most near zero, two rising steeply. Right: scatter of change in cardiac depth response against change in intensity sensitivity, rising.">
    <figcaption><span class="fignum">Figure 4</span> · <b>Odor intensity sensitivity and cardiac
    coupling.</b> Left: most patients cannot distinguish odor from blank air at baseline; two recover
    sharply. Right: change in the cardiac response to breath depth tracks change in intensity
    sensitivity (ρ = 0.71), as strongly as it tracks the full three-task composite.</figcaption>
  </figure>
  <div class="key">
    <span class="lab">Why this measure in particular</span>
    <p>Identification requires naming and semantic access; discrimination requires a comparison
    judgment. Rating how intense an odor smells is the closest available readout of transduction at the
    epithelium. If restored nasal sensory function is what strengthens breath-locked afferent drive, the
    intensity measure is where the relationship should appear — and it does, undiluted.</p>
  </div>
</section>

<section>
  <div class="col">
    <h2>What the proposed work would establish</h2>
    <p>These are seven patients and a descriptive relationship. They are sufficient to show the measures
    are obtainable, the variance is there to exploit, and the predicted direction is the one observed.
    The proposed work would test it properly.</p>
    <ul>
      <li><b>Power the primary test.</b> A single pre-specified coupling measure in a cohort several times
      this size, with the analysis fixed in advance rather than selected from a panel.</li>
      <li><b>Separate the sensory account from the inflammatory one.</b> Systemic inflammation plausibly
      affects olfactory epithelium and autonomic function independently. Adding inflammatory markers and
      testing whether coupling tracks olfactory outcome after adjusting for them distinguishes a shared
      systemic factor from a breathing-specific mechanism.</li>
      <li><b>Test the nasal-route prediction directly.</b> If the effect depends on nasal afferent drive,
      it should attenuate during oral breathing. This is a within-patient contrast and the strongest
      available discriminator.</li>
      <li><b>Extend the time course.</b> The relationship is clearest over the first month. Whether
      coupling tracks olfaction over six months, where trajectories diverge and some patients lose ground,
      is a separate and clinically interesting question.</li>
      <li><b>Standardise acquisition.</b> A brief seated spontaneous-breathing block in every session
      would substantially improve the precision of every coupling measure at negligible cost.</li>
    </ul>
  </div>
</section>

<footer>
  Seven patients with baseline and ~1 month sessions (28–42 days). Respiration recorded by nasal cannula,
  ECG continuous, both at 500 Hz. Coupling measures computed breath by breath over spontaneous-breathing
  blocks and referenced to within-session surrogate distributions. Olfactory composite from odor
  discrimination (d′), free identification, and air-calibrated intensity threshold. Correlations are
  Pearson (r) and Spearman (ρ) as marked; no inferential claim is made at this sample size.
</footer>
</div>
"""
os.makedirs(os.path.join(D, "reports"), exist_ok=True)
p = os.path.join(D, "reports", "grant_report.html")
open(p, "w", encoding="utf-8").write(HEAD + BODY)
print("wrote", os.path.getsize(p)//1024, "KB")
