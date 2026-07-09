# G25 → Unsupervised Ancestry Decomposition

A browser-based tool that estimates ancestry mixtures from [Global25 (G25)](https://eurogenes.blogspot.com/2025/02/g25-available-again.html) PCA coordinates, producing the same kind of stacked bar charts you'd get from ADMIXTURE — but without telling it which populations to look for ahead of time.

No installation, no server, nothing to configure. Open the HTML file in any modern browser and it works.

![K=4 decomposition of 32 world populations showing genuine soft ancestry mixtures across European, Near Eastern, East Asian, and Sub-Saharan components](screenshot.png)

---

## What this actually does

If you've used a tool like ADMIXTURE before, you already know the output: a chart where each person is a bar, and the bar is divided into colored segments — "40% this, 30% that, 30% the other thing." What you might not know is how much heavy machinery normally goes into producing that chart. ADMIXTURE and STRUCTURE work directly on raw genotypes, hundreds of thousands of SNPs, and they lean on real population-genetics assumptions like Hardy-Weinberg equilibrium to do it.

This tool skips all of that. It starts from G25 coordinates instead — 25 numbers per person that are already a compressed, PCA-based summary of genetic variation — and asks whether ancestry structure can be recovered from that compressed representation alone, with no reference populations chosen in advance. The components just have to emerge from the data.

That's useful if you want to:

- Poke at population structure in a batch of G25 samples without picking reference groups first
- Sanity-check what an ADMIXTURE run might show, without running the actual genotype pipeline
- See how populations relate to each other as mixtures rather than fixed categories
- Try different K values (number of components) interactively and watch the structure shift

### What it isn't

This is not a replacement for SNP-based ADMIXTURE, and it won't pretend to be. It works in a reduced 25-dimensional space, not on allele frequencies, so there's no population-genetic model underneath it — no HWE, no linkage assumptions. The "components" it finds are directions in PCA space, not literal ancestral populations, and any fine-scale structure that needs hundreds of thousands of SNPs to resolve simply isn't there to find. Treat the results as analogous to a real ADMIXTURE run, not identical to one.

---

## The algorithm, in two passes

**Plain version:** the tool is trying to explain each person's G25 coordinates as a mix of a small number of "ancestral profiles," where the mixing percentages for each person have to be non-negative and add up to 100%. It alternates between two questions — "given everyone's mixing percentages, what do the ancestral profiles look like?" and "given the profiles, what's each person's best mix?" — and keeps refining both until the numbers stop moving.

**Technical version**, for anyone who wants the actual method:

It's **Alternating Least Squares (ALS) with simplex-projected mixture proportions**.

Given an input matrix **X** (N samples × D = 25 coordinates), the tool solves for:

- **Q** (N × K) — per-sample mixture proportions, each row constrained to the probability simplex (non-negative, sums to 1)
- **P** (K × D) — ancestral component centroids in G25 space, left unconstrained since PCA coordinates can be negative

such that **X ≈ Q · P**, minimizing Frobenius reconstruction error.

The simplex constraint on Q is doing the real work here — it's what mirrors ADMIXTURE's requirement that ancestry fractions be non-negative and sum to one, and it's why you get genuine soft mixtures instead of the hard 0%/100% calls a Gaussian Mixture Model would hand you in high-dimensional space with few samples.

**The two update steps:**

1. **Fix Q, solve for P.** With Q held constant, this is ordinary least squares: P = (QᵀQ)⁻¹QᵀX, with a small ridge term (λ = 10⁻⁷) added to QᵀQ for stability. The K × K matrix is inverted by Gauss-Jordan elimination, which is fine since K stays small (typically 2–20).

2. **Fix P, solve for Q.** Each row of Q updates independently by minimizing ‖xᵢ − qᵢP‖² subject to qᵢ living on the simplex. This runs as projected gradient descent: compute the gradient (PᵀP · qᵢ − Pᵀxᵢ), take a step, then project back onto the simplex using Duchi et al.'s (2008) algorithm. Step size is 0.9 / tr(PᵀP), and each row gets 80 inner iterations.

**Initialization:** component centroids P start from **K-means++ seeding** — the first centroid picked at random, each next one chosen with probability proportional to its squared distance from the nearest existing centroid. This spreads the starting points out and makes the result far less sensitive to which random seed you happened to get, compared to plain random initialization. Q starts uniform (1/K per component) plus a little noise, then gets projected onto the simplex.

**Multiple restarts:** the objective isn't convex, so it can get stuck in local minima. The tool runs several independent restarts with different seeds and keeps whichever one reconstructs the data best.

**Convergence:** the algorithm stops once the change in mean squared reconstruction error between iterations drops below 10⁻¹², or when it hits the iteration cap, whichever comes first.

---

## Parameters

| Parameter | Default | Range | Effect |
|-----------|---------|-------|--------|
| **Mode** | Single K | Single K / Sweep K range | Run one K value, or test a range and compare them |
| **K** | 4 | 2–20 | Number of components to discover. Low K gives you a broad split; high K resolves finer structure but starts overfitting past a point. K=2 usually shows a continental-scale split; K=4–6 starts revealing sub-continental structure |
| **K max** | 10 | 3–20 | Upper end of the K range in sweep mode (the lower end is fixed at 2) |
| **Selection Criterion** | BIC | BIC / AIC / Raw MSE / Cross-Validation | *(Sweep mode only)* How the "best" K is picked once the sweep finishes. See below |
| **CV Folds** | 5 | 2–10 | *(Cross-Validation only)* Number of random train/test splits averaged per K |
| **Holdout %** | 20 | 10–40 | *(Cross-Validation only)* Fraction of samples set aside as the test set in each fold |
| **Iterations** | 500 | 50–5000 | Max ALS iterations per run. 500 is normally enough to converge; bump it up for large or messy datasets |
| **Restarts** | 10 | 1–30 | Independent random restarts per K. More restarts lower the odds of landing in a bad local minimum, at the cost of runtime scaling linearly. 10 is a reasonable default; go to 20–30 if you want publication-quality stability |
| **Palette** | Classic | Classic / Earth Tones / Vibrant | Just the colors. No effect on the analysis |

### Picking a K

There's no "correct" K, same as with ADMIXTURE — different values just show you different levels of structure:

- **K=2:** the coarsest split (e.g., West Eurasian vs. everything else, in the sample data)
- **K=3–4:** major continental groupings show up
- **K=5–8:** sub-continental patterns emerge — Northern vs. Southern European, Near Eastern vs. South Asian, that sort of thing
- **K>10:** increasingly fine-grained, and increasingly likely to be fitting noise if your sample size is small

In sweep mode, the tool runs every K in the range and tries to find the **elbow** — the last K where adding one more component still helped by a meaningful amount, and after which the curve settles down and stays settled for the rest of the range you tested. That second part matters: a single small step surrounded by bigger ones isn't enough to count. Early versions of this tool stopped at the *first* small step it saw, which turned out to be too easy to fool — a curve can dip briefly by chance and then keep making real improvements for several more K values, especially with a smaller sample count where different restarts land on slightly different local optima and the curve gets noisier step to step. The current version instead finds the *last* step that's still big, and only calls it an elbow if nothing big happens again afterward. If the curve never actually settles within your tested range — still jumping around, or still trending down, all the way to K max — the tool says so plainly ("no elbow, showing K=X, likely just the top of your range") rather than quietly pretending it found one. The chart tags the result **ELBOW** or **NO ELBOW** so you always know which one you're looking at.

Which curve the elbow gets found on depends on the **Selection Criterion**:

- **BIC (default).** Reconstruction error alone always keeps dropping as K grows, so picking on raw error would just point you at the largest K in your range. BIC adds a penalty for extra free parameters — treating the model as having K·D of them (just the component profiles, P) — and picks the K with the best penalized score. It penalizes harder than AIC, so it tends to land on a more conservative, parsimonious K.
- **AIC.** Same idea, lighter penalty. Tends to allow a somewhat higher K than BIC would.
- **Raw MSE.** No penalty at all — mostly useful for seeing the unpenalized error curve, not for picking K on its own.
- **Cross-Validation.** Doesn't rely on a penalty formula at all — see below.

Whichever criterion is active, the chart plots that score against K and highlights the elbow (or the fallback bar, tagged accordingly); click any bar to jump to that K's decomposition.

One caveat worth knowing about BIC/AIC specifically: only P's K·D component-profile values count as free parameters here — Q's per-sample mixture proportions don't, since they're treated as latent per-sample estimates rather than fixed model parameters (the same convention classical factor analysis uses for factor scores). An earlier version of this tool counted Q's N·(K−1) proportions too, which seems more "honest" at first glance, but backfires badly in practice: with a few hundred samples or more, that term dominates the penalty and makes it grow by roughly the same fixed amount at every step, regardless of K, which reliably swamps the shrinking fit improvement from any real additional structure and pushes the "best" K straight to the bottom of the range — not because the data lacks structure, but because the formula was structurally biased toward the smallest K available. Counting only P avoids that. The residual risk runs the other way: with very few samples relative to K, the penalty can still be too weak to catch real over-fitting, since Q is free to fit almost anything once K gets close to N. The tool detects this (when K·D exceeds half the available data points) and prints a warning under the chart. When you see it, treat the visual elbow and whether the resulting groups make biological sense as more trustworthy than the highlighted bar.

### Cross-Validation (a sturdier, slower alternative)

BIC, AIC, and raw MSE all share a blind spot: they only ever measure how well a K fits the samples it was trained on. Past a certain K — especially once components start approaching the number of coordinate dimensions (25) or the number of samples — a model can start reconstructing individual samples' idiosyncrasies almost perfectly, and none of those three criteria can tell the difference between that and real structure, because they never check against anything the model didn't see.

**Cross-Validation** does check against something it didn't see. For each K, the tool runs several folds; each fold randomly holds out a slice of your samples (the **Holdout %**), fits the K components on everyone else, then freezes those components and projects the held-out samples onto them — the same simplex-constrained least-squares step ALS already uses internally, just applied to samples that had no influence on the components being tested. The reconstruction error on that held-out slice, averaged across folds (**CV Folds**), is the score for that K. A K that only looked good by memorizing quirks of its training samples won't transfer to samples it never saw, so held-out error stops improving — or gets worse — right around where real structure runs out, without needing to count parameters or guess at a penalty at all. The same elbow-detection logic described above applies to this curve too.

The cost is speed: CV effectively refits the model from scratch for every fold at every K, so it's roughly `folds ×` as slow as a normal sweep. To keep that bounded, each fold's fit uses at most 5 restarts internally regardless of the main **Restarts** setting. Once CV picks a best K, the tool re-fits that K one more time on your *entire* dataset (using your normal **Restarts** setting) to produce the actual decomposition you see — the CV folds themselves are only ever used to pick K, never to produce the displayed result. Clicking a different bar in CV mode triggers this same full re-fit (with live restart progress in the status line), which is slower than clicking a bar under BIC/AIC/MSE — those already have a complete decomposition cached for every K, since their fits used the whole dataset from the start.

If a fold's held-out set ends up small (under 15 samples), the tool flags that the CV estimate may be noisy — worth raising the fold count, the holdout percentage, or your overall sample size before trusting fine differences between bars.

---

## Input format

Standard G25/Vahaduo format — one sample per line, comma-separated:

```
Population_Label,coord1,coord2,coord3,...,coord25
```

For example:

```
English,0.0389,0.0596,0.0497,0.0167,0.003,-0.003,0.005,0.0014,...
Greek_Crete,0.0609,0.0845,0.0285,-0.0065,-0.0135,0.003,-0.0095,...
Yoruba,0.082,-0.012,-0.068,-0.059,0.01,-0.018,0.035,0.025,...
```

- Labels can use letters, numbers, underscores, and hyphens
- Coordinates can be positive or negative — G25 space is centered
- Lines starting with `#` are comments and get skipped
- The tool figures out the number of dimensions on its own (doesn't have to be exactly 25)
- A sample dataset of 32 world populations is built in — hit **Load Sample Data**

---

## Output

### Stacked bar chart

The classic ADMIXTURE-style plot: one bar per sample, colored segments showing each component's proportion. Three sort options:

- **Input** — original order from the pasted data
- **Cluster** — grouped by dominant component, then sorted by proportion within each group (usually the clearest view)
- **Name** — alphabetical

### Component proportions table

Every sample, its exact percentage on each component, and which component is dominant. Follows whatever sort order the chart is using.

### CSV export

**Export CSV** downloads the full proportions matrix, one row per sample, percentages from 0–100 for each component — ready to drop into R, Python, Excel, or wherever else you need it.

### Model selection chart (sweep mode)

A bar chart scored by the active **Selection Criterion** (BIC, AIC, raw MSE, or Cross-Validation) — one bar per K, with the elbow (or fallback) bar highlighted and tagged **ELBOW** / **NO ELBOW**. Click any bar to view that K — for BIC/AIC/MSE this is instant, since every K's full decomposition is already cached; for Cross-Validation it triggers a quick re-fit on your whole dataset, with restart progress shown in the status line. Whichever bar you're currently viewing gets a white outline, separate from the highlight color — so if you click away from the elbow to compare a different K, it stays clear which one the algorithm picked versus which one you're looking at. Lower is always better, regardless of which criterion is showing.

If you change **Restarts**, **K max**, or the criterion and re-run, the chart dims with a "Recalculating…" label until the new sweep finishes, so it's never ambiguous whether you're looking at a fresh result or a leftover one from before.

Sweep mode runs off the main thread in a Web Worker, so the page stays responsive while it works — the status line shows live progress ("Calculating K=6 of 10 (4/10 restarts)...", or fold-by-fold progress for Cross-Validation) instead of freezing, and a **Cancel** button appears if you want to stop partway through.

---

## Interpreting results with AI

![Interpret K=10 components to understand which populations might represent the clusters. Creates a LLM prompt.](decomposition_interpret_with_ai.jpg)

The decomposition tells you the components exist. It doesn't tell you what they mean. If a component peaks at 45% in Yoruba and Dinka samples, that's obviously capturing Sub-Saharan African ancestry — but "obviously" only holds if you already know something about population genetics. The algorithm doesn't.

This is where an LLM — [Claude](https://claude.ai), [ChatGPT](https://chat.openai.com), [Gemini](https://gemini.google.com), whichever you have handy — earns its keep. It can look at which populations load highest on each component and match that pattern against what's known from the ancient DNA literature.

### Built-in prompt generator

Once you've run a decomposition, an **Interpret with AI** panel shows up below the results. Click **Generate Prompt** and it builds a structured summary — population averages and top-scoring samples per component — formatted for any LLM to read. Copy it, paste it into your model of choice, done.

### Doing it by hand

If you'd rather write your own, export the CSV and adapt something like this:

```
I ran an unsupervised ancestry decomposition on G25 scaled PCA coordinates
with K=6 components. G25 coordinates are 25-dimensional scaled PCA projections
that capture the major axes of human genetic variation.

No source populations were predefined — the components emerged unsupervised.
Below are the component proportions for each sample (percentages, rows sum
to 100%).

[paste CSV contents here]

Please interpret what each of the 6 components most likely represents in
population-genetic terms (e.g., "Western Hunter-Gatherer," "Near Eastern
Farmer," "East Asian," etc.) based on which samples load highest on each
component.
```

### Getting better interpretations

- **Use informative population labels.** "Yoruba" and "Han_Chinese" give the model something to work with. "Sample_001" doesn't.
- **Say what K you used and give context on your dataset.** If everything's European, mention that — it changes how the model should read the components.
- **Follow up.** Once the components are labeled, you can ask things like "which of my samples show the most Near Eastern Farmer ancestry?" or "what's the historical story behind the split between Component 3 and Component 5?"
- **Any capable model works.** Free tiers of Claude, ChatGPT, and Gemini all know enough population genetics and ancient DNA to be useful here.

---

## Technical notes

- **Why not a Gaussian Mixture Model?** In 25 dimensions with relatively few samples, a GMM will plant a Gaussian right on top of each cluster and drive posterior responsibilities to 100%/0% — hard assignments, not the soft mixtures you actually want. The simplex-constrained factorization avoids this by construction, not by tuning.

- **Negative coordinates.** G25 coordinates can go negative, which rules out standard NMF. Here, only Q is forced non-negative; P stays unconstrained and can sit anywhere in coordinate space.

- **Simplex projection.** Projecting an arbitrary vector onto the probability simplex uses the O(n log n) method from Duchi, Shalev-Shwartz, Singer, and Chandra (2008).

- **Cost.** Runtime scales as O(restarts × iterations × N × K × D) for the Q-update. With defaults (10 restarts, 500 iterations, D=25), expect well under a second for datasets up to a few hundred samples, and a few seconds for datasets in the thousands.

- **Threading.** Sweep mode runs in a Web Worker built from an inline `Blob`, so it's still a single self-contained HTML file — no separate `.js` file to ship alongside it. The worker gets its own copy of the ALS code, since workers don't share scope with the page; single-K mode still runs on the main thread, since one K with a handful of restarts finishes fast enough that it isn't worth the extra plumbing.

- **BIC/AIC parameter counting.** Both are computed as `n·ln(RSS/n) + penalty·numParams`, with `numParams = K·D` (just P's component profiles) and `n = N·D`. Q's mixture proportions are deliberately *not* counted, even though there are N·(K−1) of them, because they're per-sample latent estimates rather than fixed model parameters — the same convention used for factor scores in classical factor analysis. Counting them (as an earlier version of this tool did) makes the parameter total scale with sample size N, and for datasets much beyond a couple hundred samples that term dominates the penalty and forces the criterion toward the smallest K in the range regardless of the actual error curve. Leaving Q out avoids that, at the cost of a milder residual risk in the other direction — with very few samples relative to K, the penalty can still be too weak to catch real over-fitting, which is what the in-app warning watches for (see "Picking a K" above).

- **Elbow detection.** Given a sequence of scores across K, the tool computes the step-to-step change and finds the *last* step whose size is still at least 15% of the first (largest) step. The elbow is the K that step landed on — and it only counts as an elbow if no later step is big again, i.e. every step after it stays under that 15% threshold. If the very last step tested is still big, or the curve never once drops below the threshold, there's no elbow, and the tool falls back to the single lowest-scoring K with a note that it's likely just the top of the range. An earlier version stopped at the *first* small step instead of requiring everything after it to also stay small; that turned out to be unreliable on noisy curves, where a single small step can appear early by chance while real, substantial improvement is still ahead. Requiring the flattening to be sustained through the rest of the tested range fixed that, at the cost of reporting "no elbow" more often — which, on the noisy curves this tool tends to produce, is usually the more honest answer.

- **Cross-Validation implementation.** Each fold: shuffle sample indices with the tool's seeded RNG, hold out `Holdout %` of them, fit P via the normal ALS restarts (capped at 5 per fold regardless of the main Restarts setting) on the remaining training samples, then freeze P and run a single projected-gradient Q-update — the same math `runALS` uses internally for its per-sample step — for each held-out sample against that frozen P. Held-out reconstruction error is averaged across all folds and all K values tested. This intentionally never touches Q for training samples across folds or reuses a fit between K values; every (K, fold) pair is an independent fit from a fresh k-means++ initialization, which is why CV costs roughly `folds ×` a normal sweep.

- **Determinism.** Same parameters, same result, every time — the tool uses a seeded LCG rather than the browser's built-in randomness. Change K or the restart count and you get a different seed, so results shift accordingly.

---

## Running it

Open `g25-admixture-tool.html` in Chrome, Firefox, Safari, or Edge. It's all client-side JavaScript — nothing you paste in ever leaves your machine.

---

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE) for the full text.

---

## References

- Alexander, D.H., Novembre, J., & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.
- Duchi, J., Shalev-Shwartz, S., Singer, Y., & Chandra, T. (2008). Efficient Projections onto the ℓ₁-Ball for Learning in High Dimensions. *ICML 2008*.
- Pritchard, J.K., Stephens, M., & Donnelly, P. (2000). Inference of population structure using multilocus genotype data. *Genetics*, 155(2), 945–959.
- Lazaridis, I. et al. (2014). Ancient human genomes suggest three ancestral populations for present-day Europeans. *Nature*, 513(7518), 409–413.
- Global25 (G25) coordinates by Davidski: [Eurogenes Blog](https://eurogenes.blogspot.com/2025/02/g25-available-again.html)
