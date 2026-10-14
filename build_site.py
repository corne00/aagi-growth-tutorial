"""Build the static documentation pages. Numerical results come from R only."""
from pathlib import Path
import html
import json

ROOT = Path(__file__).resolve().parent
PAGES = [
 ("index", "Introduction", "Why model growth dynamically?"),
 ("dataset", "APSIM data", "The APSIM wheat dataset"),
 ("response", "Response functions", "Approximating the response"),
 ("hierarchy", "Genotype responses", "Shared and genotype-specific responses"),
 ("design", "Design matrix", "From trajectories to a design matrix"),
 ("gibbs", "Gibbs sampling", "Bayesian fitting with Gibbs sampling"),
 ("prediction", "Prediction", "Predicting biomass in held-out years"),
 ("uncertainty", "Uncertainty", "What does an uncertainty interval describe?"),
 ("interpretation", "ALE interpretation", "Interpreting environmental responses"),
 ("gxe", "G × E", "Where do G, E and G × E appear?"),
 ("reduction", "SVD extension", "Reducing correlated features"),
 ("reproduce", "Reproduce", "Run the analysis yourself"),
]

def equation(text):
    return '<div class="equation">\\[' + text + '\\]</div>'

def code(name, title="Show R code"):
    return f'<details class="code"><summary>{title}</summary><button class="copy" type="button">Copy</button><pre><code data-function="{name}"></code></pre></details>'

def plot(name, caption="", height=""):
    return f'<figure><div class="plot" id="{name}" {height} role="img" aria-label="{html.escape(caption or name)}"></div><figcaption>{caption}</figcaption></figure>'

def controls(content):
    return '<div class="controls">' + content + '</div>'

def select(name, label, options=None):
    opts = ''.join(f'<option value="{v}">{l}</option>' for v,l in (options or []))
    return f'<label>{label}<select id="{name}">{opts}</select></label>'

def slider(name, label, lo, hi, value, step=1):
    return f'<label class="range">{label} <output for="{name}" id="{name}-value">{value}</output><input type="range" id="{name}" min="{lo}" max="{hi}" step="{step}" value="{value}"></label>'

common = controls(select("genotype", "Genotype") + select("location", "Location") + select("year", "Year"))
content = {}
content["index"] = '''<p>Final biomass does not describe when growth occurred or how growth responded to the environment. Different growth trajectories can reach similar endpoints.</p>'''
content["index"] += plot("endpoints", "Actual simulation trajectories, selected for similar final biomass. Values and timing are unchanged; endpoints are not forced to coincide.")
content["index"] += '''<p>Repeated measurements let us study a response function relating environmental conditions to the rate of biomass accumulation.</p>'''
content["index"] += equation(r"\frac{dB_g(t)}{dt}=f_g(\mathbf u(t)),\qquad B_g(t)=B_g(0)+\int_0^t f_g(\mathbf u(s))\,ds.")
content["index"] += '''<div class="process" aria-label="Environment to response to biomass"><span>Environment<br><small>driver sequence</small></span><span class="arrow">→</span><span>Response function<br><small>biomass growth rate</small></span><span class="arrow">→</span><span>Growth trajectory<br><small>integrated growth</small></span></div>'''
content["index"] += controls(slider("time", "Inspect the growing season", 0, 140, 70))
content["index"] += plot("flow", "One real weather sequence, its held-out fitted response, and accumulated biomass. The vertical line follows the same day in all panels.")
content["dataset"] = '''<p>We use APSIM-wheat simulations from Bustos-Korts and colleagues: 199 genotypes, four Australian locations and 31 years. Biomass is available daily.</p><p>APSIM simulates crop growth from physiological parameters, weather, soil and management. Our reduced model uses biomass observations and weather-derived predictors.</p>
<div class="apsim-flow"><div>Genotype physiology<br>Weather<br>Soil and management</div><span>→</span><div>APSIM-wheat</div><span>→</span><div>Biomass over time<br>Other crop traits</div></div><p id="dataset-counts" class="caption"></p>'''
content["dataset"] += controls(select("year", "Year") + select("genotype", "Genotype 1") + select("genotype2", "Genotype 2") + select("genotype3", "Genotype 3") + select("driver", "Weather", [("temperature","Mean temperature"),("Radn","Radiation"),("Rain","Rainfall")]))
content["dataset"] += '<div class="plot-grid">' + ''.join(plot("data-"+l, l) for l in ["Emerald","Merredin","Narrabri","Yanco"]) + '</div>'
content["dataset"] += plot("weather", "Observed weather at the four locations. Genotypes share the weather sequence within each environment.")
content["dataset"] += '''<p>We do not use markers, APSIM soil-water status, leaf area or other simulated crop states as predictors. Thermal time and trailing rainfall are calculated from weather.</p><h2>Data tables</h2><div id="table-preview"></div>'''
content["dataset"] += code("read_apsim") + '''<p class="reference">Bustos-Korts, D. et al. (2019). From QTLs to Adaptation Landscapes: Using Genotype-To-Phenotype Models to Characterize G×E Over Time. <i>Frontiers in Plant Science</i> 10:1540. <a href="https://doi.org/10.3389/fpls.2019.01540">Article</a> · <a href="https://doi.org/10.6084/m9.figshare.24843705.v1">Dataset, CC BY 4.0</a>.</p>'''
content["response"] = '''<p>We approximate the response using candidate functions of the environmental drivers. Splines describe curved responses; product terms allow interactions. The response is nonlinear in its drivers and linear in its coefficients.</p>'''
content["response"] += equation(r"f(\mathbf u)=\sum_{k=1}^{K}\beta_k\phi_k(\mathbf u).")
content["response"] += controls(select("feature", "Inspect a fitted component"))
content["response"] += plot("basis", "A basis function evaluated along real Yanco 2010 weather; its fitted contribution; the complete population response. Coefficients are posterior means from the full-data fit.")
content["response"] += '''<p>The compact library contains 16 features: an intercept, three natural-spline functions each for temperature, thermal time, radiation and 14-day rainfall, a day-length term, and two thermal-time interactions.</p>'''
content["response"] += code("make_basis", "Show R code: training basis") + code("evaluate_basis", "Show R code: evaluate features")
content["hierarchy"] = '''<p>Genotypes share a population response and have their own deviations. Estimating these responses jointly allows well-observed genotypes to inform genotypes with fewer observations.</p>'''
content["hierarchy"] += equation(r"f_g(\mathbf u)=\boldsymbol\phi(\mathbf u)^\top\boldsymbol\beta+\boldsymbol\phi(\mathbf u)^\top\mathbf b_g,\qquad b_{gk}\mid\lambda_{Gk}\sim N(0,\lambda_{Gk}^{-1}).")
content["hierarchy"] += controls(select("genotype", "Genotype"))
content["hierarchy"] += plot("hierarchy", "Population and genotype responses along Yanco 2010 weather. These are full-fit response summaries, not an out-of-sample evaluation.")
content["hierarchy"] += '''<p>The zero mean is a prior on genotype deviations; it is not a constraint forcing the fitted deviations to sum to zero. We interpret complete genotype responses as well as their shared component.</p><h2>What does pooling change?</h2><p>In a real-data exercise, g008, g052 and g185 have only one or three training environments. Nine other genotypes retain their training data. Both methods are evaluated on the same held-out years for the three sparsely observed genotypes.</p>'''
content["hierarchy"] += controls(select("density", "Training environments for sparse genotypes", [("one","One"),("three","Three")]))
content["hierarchy"] += plot("pooling", "Held-out RMSE for independent and hierarchical fits. This exercise uses 12 real genotypes; it is separate from the full-data analysis.")
content["hierarchy"] += '''<details><summary>Location extension</summary>'''+equation(r"f_{g\ell}(\mathbf u)=\boldsymbol\phi(\mathbf u)^\top(\boldsymbol\beta+\mathbf b_g+\mathbf c_\ell).")+'''<p>Location deviations use sum-to-zero contrasts and a fixed Gaussian shrinkage prior. They describe persistent differences among known locations. This compact extension omits the separate genotype–location offset in the research slides.</p><div id="location-results"></div></details>'''
content["hierarchy"] += code("joint_genotype_update", "Show R code: joint shared/genotype update")
content["design"] = '''<p>We model changes in biomass over observation intervals. Integrating each feature along the weather sequence produces a design-matrix row, without numerically differentiating the biomass observations.</p>'''
content["design"] += equation(r"\Delta B_{gej}=\sum_k\underbrace{\int_{t_j}^{t_{j+1}}\phi_k(\mathbf u_e(s))\,ds}_{X_{ejk}}(\beta_k+b_{gk})+\varepsilon_{gej}.")
content["design"] += controls(slider("interval", "Observation interval", 0, 10, 4))
content["design"] += plot("interval-growth", "The selected biomass increment and its matching integrated-feature row, for g008 in Yanco 2010.")
content["design"] += plot("matrix", "Integrated and scaled feature matrix. Rows are observation intervals; columns are candidate functions.")
content["design"] += equation(r"\mathbf y=\mathbf Z\boldsymbol\theta+\boldsymbol\varepsilon.")
content["design"] += '''<p>Each shared feature is accompanied by genotype-specific coefficient blocks. We compute sufficient statistics once, so the sampler does not repeatedly multiply the complete regression matrix.</p><p>Observation spacing is 14 days, retaining the final shorter interval. We use trapezoidal quadrature on the daily weather grid and training-only column scaling.</p>'''
content["design"] += code("observation_intervals") + code("integrated_features") + code("prepare_design")
content["gibbs"] = '''<p>Gibbs sampling alternates between coefficients, shrinkage precisions and residual variance. Every retained joint draw defines a response function.</p>'''
content["gibbs"] += equation(r"\mathbf y\mid\boldsymbol\theta,\sigma^2\sim N(\mathbf Z\boldsymbol\theta,\sigma^2I),\quad\theta_k\mid\lambda_{q(k)}\sim N(0,\lambda_{q(k)}^{-1}).")
content["gibbs"] += '''<div class="process"><span>Coefficients<br><small>Gaussian update</small></span><span class="arrow">→</span><span>Precisions<br><small>Gamma update</small></span><span class="arrow">→</span><span>Residual variance<br><small>inverse-Gamma update</small></span></div>'''
content["gibbs"] += controls(select("parameter", "Parameter", [("beta","Population intercept coefficient"),("precision","Genotype intercept precision"),("sigma","Residual standard deviation")]) + slider("iteration", "Retained draws per chain", 100, 1200, 1200, 100))
content["gibbs"] += plot("traces", "Four independently initialised chains from the full-data hierarchy. Warmup draws are excluded.")
content["gibbs"] += '<div id="diagnostic-table"></div>'
content["gibbs"] += '''<p>Coefficients and genotype deviations are sampled jointly using block Gaussian algebra. ARD continuously shrinks coefficients; its precision is not an inclusion probability.</p><details><summary>Priors and conditional equations</summary>'''
content["gibbs"] += equation(r"\lambda_{0k}\sim\operatorname{Gamma}(2,1),\quad\lambda_{Gk}\sim\operatorname{Gamma}(2,0.02),\quad\sigma^2\sim\operatorname{InvGamma}(2,0.1).")
content["gibbs"] += '''<p>Gamma parameters are shape and rate; inverse-Gamma parameters are shape and scale. Coefficient priors are independent of residual variance. All coefficients use the documented integrated-feature scaling.</p>'''
content["gibbs"] += equation(r"V=(\mathbf Z^\top\mathbf Z/\sigma^2+\Lambda)^{-1},\quad m=V\mathbf Z^\top\mathbf y/\sigma^2,\quad\boldsymbol\theta\mid\cdots\sim N(m,V).")
content["gibbs"] += '</details>' + code("draw_gaussian") + code("fit_gibbs", "Show complete R sampler")
content["prediction"] = '''<p>We evaluate each posterior response along a held-out weather sequence and integrate it to predict biomass. The genotype and location have been observed in other years.</p>'''
content["prediction"] += equation(r"B^{(s)}_{ge}(t)=B_{ge}(0)+\int_0^t f^{(s)}_g(\mathbf u_e(v))\,dv.")
content["prediction"] += common + controls(select("model", "Model", [("hierarchical","Shared + genotype responses"),("shared","Shared response only")]))
content["prediction"] += '<p id="prediction-split" class="caption"></p>' + plot("prediction", "Weather, predicted growth rate and biomass. APSIM biomass is shown at its native daily times; model predictions are evaluated at observation checkpoints.")
content["prediction"] += '<h2>Held-out performance</h2><div id="metrics-table"></div>'
content["prediction"] += '''<p>Entire year blocks are withheld. The basis, knots and scaling are learned from training years only. This tests transfer across years at known locations and genotypes; it does not test new genotypes or new locations.</p>'''
content["prediction"] += code("posterior_draws") + code("evaluate_fit")
content["uncertainty"] = '''<p>Credible intervals describe uncertainty in the fitted biomass trajectory. Predictive intervals additionally include the interval discrepancy specified by our likelihood.</p>'''
content["uncertainty"] += common + controls(select("level", "Interval", [("0.5","50%"),("0.8","80%"),("0.95","95%"),("0.99","99%")]) + select("kind", "Type", [("credible","Credible"),("predictive","Predictive")]))
content["uncertainty"] += plot("uncertainty", "Joint posterior trajectories and intervals in held-out years. Initial biomass is fixed to the observed zero at sowing.")
content["uncertainty"] += plot("coverage", "Terminal-biomass coverage across all held-out trajectories, separately by year block. The diagonal denotes nominal coverage.")
content["uncertainty"] += '''<p>The simulator data contain no added measurement noise. Residual variance represents reduced-model discrepancy. Predictive discrepancies are independent between fitted intervals and accumulate along each path; credible intervals omit those discrepancies.</p><p>Coverage must be assessed, not assumed. Narrow intervals can miss APSIM trajectories when the response model is incomplete.</p>'''
content["uncertainty"] += code("posterior_draws")
content["interpretation"] = '''<p>Weather drivers are correlated. Accumulated local effects (ALE) summarise local changes in the fitted growth-rate response within the observed driver distribution.</p>'''
content["interpretation"] += controls(select("ale-genotype", "Response") + select("driver", "Driver", [("temperature","Mean temperature (°C)"),("thermal_time","Thermal time (°C days)"),("Radn","Radiation (MJ/m²/day)"),("rain14","14-day rainfall (mm)")]) + select("bins", "Bins", [("6","6"),("10","10"),("15","15")]))
content["interpretation"] += plot("ale", "Full-fit ALE of biomass growth rate, with a 95% posterior interval. Bin counts show the reference data supporting each local comparison.")
content["interpretation"] += equation(r"\delta_k=\frac{1}{n_k}\sum_{i:\,u_{ij}\in I_k}\left[f(z_k,\mathbf u_{i,-j})-f(z_{k-1},\mathbf u_{i,-j})\right].")
content["interpretation"] += '''<p>We accumulate the local differences and centre the curve using empirical bin weights. The reference sample consists of 6,000 evenly indexed environment-days, each weighted equally. One-dimensional ALE does not display every interaction and is not a causal effect.</p>'''
content["interpretation"] += code("ale_curve")
content["gxe"] = '''<p>Environmental sequences generate different growth conditions. Genotypes differ in their responses, so genotype contrasts can change across environments and during the season.</p>'''
content["gxe"] += controls(select("genotype", "Genotype 1") + select("genotype2", "Genotype 2") + select("year", "Year"))
content["gxe"] += plot("contrasts", "Predicted genotype differences at four locations, with the corresponding APSIM differences. All model curves come from held-out year fits.")
content["gxe"] += equation(r"M_{ge}(t)=\mu(t)+G_g(t)+E_e(t)+I_{ge}(t).")
content["gxe"] += controls(slider("time", "Common days after sowing", 14, 112, 84, 14))
content["gxe"] += plot("gxe-components", "Observed and predicted mean-square fractions, and component correlations, across the balanced 199-genotype × 124-environment panel.")
content["gxe"] += '''<p>At each common day, we subtract genotype and environment means to obtain the double-centred interaction. These are descriptive orthogonal sums-of-squares fractions, not REML variance-component estimates. No APSIM maturity information is used to align predictions.</p>'''
content["gxe"] += '<details class="code"><summary>Show R code: G/E/G×E decomposition</summary><pre><code>M &lt;- biomass_matrix\nmu &lt;- mean(M)\nG &lt;- rowMeans(M) - mu\nE &lt;- colMeans(M) - mu\nI &lt;- sweep(sweep(M, 1, G + mu, "-"), 2, E, "-")\nfraction &lt;- c(mean(G^2), mean(E^2), mean(I^2))\nfraction &lt;- fraction / sum(fraction)</code></pre></details>'
content["reduction"] = '''<p>Correlated basis functions can produce a wide, redundant design. Singular value decomposition replaces the nonconstant feature family with fewer coordinates, estimated from training data.</p>'''
content["reduction"] += equation(r"X=U\Sigma V^\top,\qquad\widetilde X=XV_r.")
content["reduction"] += controls(slider("rank", "Inspect retained feature rank", 1, 15, 8))
content["reduction"] += plot("spectrum", "Cumulative squared singular-value energy of the training design. Moving the slider inspects the projection; it does not refit the model.")
content["reduction"] += '<div id="svd-results"></div>'
content["reduction"] += '''<p>The fitted extension retains 99% of training spectral energy and is evaluated on the first held-out year block. The prior acts on the reduced coordinates, so this is a distinct model variant. Feature compression does not imply low-rank genetic effects.</p>'''
content["reduction"] += code("prepare_design")
content["reproduce"] = '''<p>All displayed analysis results are generated in R from the APSIM wheat dataset. The website serves saved summaries; it does not run a new Gibbs fit when you move a slider.</p><h2>Small example</h2><p>The repository includes daily data for 12 actual genotypes at four locations in three years. Run it first to inspect the workflow.</p><pre><code>install.packages(c("data.table", "jsonlite"))\n# From the repository directory:\n# Rscript analysis/run_example.R</code></pre><h2>Full analysis</h2><pre><code>Rscript analysis/download_data.R\nRscript analysis/run.R full\nRscript analysis/extensions.R\nRscript analysis/pooling.R\nRscript analysis/refit_failed.R\nRscript analysis/refit_location.R\nRscript analysis/check.R\nRscript analysis/export.R\nRscript analysis/finalise_export.R\npython3 build_site.py</code></pre><p>R, data.table and jsonlite are free, open-source software. Natural splines and linear algebra use packages bundled with R. ASReml and other proprietary software are not required. Python is used only to assemble the static HTML pages; all numerical analysis is R.</p><h2>Files and provenance</h2><ul><li><a href="https://github.com/corne00/aagi-growth-tutorial">Browse or download the repository</a></li><li><a href="results/manifest.json">Analysis manifest</a></li><li><a href="results/heldout_metrics.csv">Held-out metrics</a></li><li><a href="results/performance.json">Profiling and optimisation measurements</a></li><li><a href="analysis/03_gibbs.R">Complete R sampler</a></li><li><a href="analysis/export.R">Prediction, interval and figure exports</a></li><li><a href="data/example.tsv">Small real-data exercise</a></li></ul><div id="run-manifest"></div><h2>Model limitations</h2><p>This is a compact approximation of simulated biomass accumulation, not a replacement for APSIM. It has no explicit biomass feedback or positivity constraint. Its Gaussian interval discrepancies need not capture persistent model errors. The optional location and SVD extensions have first-block evaluation only.</p><h2>References</h2><p>Bustos-Korts, D. et al. (2019). From QTLs to Adaptation Landscapes. <i>Frontiers in Plant Science</i>, 10:1540. <a href="https://doi.org/10.3389/fpls.2019.01540">doi:10.3389/fpls.2019.01540</a>.</p><p>APSIM wheat simulation dataset, version 1. <a href="https://doi.org/10.6084/m9.figshare.24843705.v1">doi:10.6084/m9.figshare.24843705.v1</a>. CC BY 4.0. Biomass/weather file unchanged before the documented processing.</p><p>Verburg et al., supplied APSIM project slides: motivation for the response model, integrated regression and Bayesian hierarchy. The compact workshop fits are regenerated analyses; previous slide performance numbers are not reused.</p>'''

def split_results():
    path = ROOT / "assets/results.js"
    if not path.exists():
        return ['assets/results.js']
    data = json.loads(path.read_text().removeprefix('window.WORKSHOP_DATA=').rstrip(';\n'))
    # Deduplicate environment weather, which is shared across genotype/model curves.
    weather = {}
    for curve in data['trajectories'].values():
        key = f"{curve['location']}_{curve['year']}"
        weather[key] = curve.pop('weather')
        curve['weather_key'] = key
    data['weather'] = weather
    generated = ROOT / 'assets/generated'
    generated.mkdir(exist_ok=True)
    for p in generated.glob('*.js'): p.unlink()
    scripts = []
    base = {k:v for k,v in data.items() if k not in ['trajectories','ale','code','weather','traces']}
    base.update(trajectories={}, ale={}, code={}, weather={})
    files = [('000-base.js', 'window.WORKSHOP_DATA=', base)]
    for field in ['trajectories','ale','code','weather']:
        chunk = {}; number=0
        for key,value in data[field].items():
            if chunk and len(json.dumps(chunk,separators=(',',':'))) + len(json.dumps(value)) > 60000:
                files.append((f'{field}-{number:03}.js', f'Object.assign(window.WORKSHOP_DATA.{field},', chunk)); chunk={};number+=1
            chunk[key]=value
        if chunk: files.append((f'{field}-{number:03}.js', f'Object.assign(window.WORKSHOP_DATA.{field},',chunk))
    files.append(('traces.js','window.WORKSHOP_DATA.traces=',data['traces']))
    for filename,prefix,obj in files:
        suffix = ');' if prefix.startswith('Object.assign') else ';'
        (generated / filename).write_text(prefix+json.dumps(obj,separators=(',',':'),ensure_ascii=False)+suffix)
        scripts.append('assets/generated/'+filename)
    return scripts

scripts = split_results()
for number,(slug,short,title) in enumerate(PAGES):
    navigation = ''.join(f'<a href="{s}.html"'+(' aria-current="page"' if s==slug else '')+f'><span>{i+1:02}</span>{html.escape(label)}</a>' for i,(s,label,_) in enumerate(PAGES))
    previous = f'<a href="{PAGES[number-1][0]}.html">← {PAGES[number-1][1]}</a>' if number else '<span></span>'
    following = f'<a href="{PAGES[number+1][0]}.html">{PAGES[number+1][1]} →</a>' if number<len(PAGES)-1 else '<a href="index.html">Introduction →</a>'
    document = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="robots" content="noindex,nofollow"><title>{html.escape(title)} · Dynamic crop growth</title><link rel="stylesheet" href="assets/workshop.css"><link rel="stylesheet" href="assets/vendor/katex/katex.min.css"></head>
<body data-page="{slug}"><a class="skip" href="#main">Skip to content</a><header><a class="brand" href="index.html">Dynamic crop growth</a><span>Bayesian SINDy · AAGI</span><a href="https://github.com/corne00/aagi-growth-tutorial">R code on GitHub</a></header><div class="layout"><nav aria-label="Workshop pages">{navigation}</nav><main id="main"><p class="chapter">{number+1:02} / {len(PAGES):02}</p><h1>{html.escape(title)}</h1>{content[slug]}<div class="page-links">{previous}{following}</div><footer>APSIM wheat · Reproducible Bayesian response modelling · <a href="reproduce.html">Data, methods and references</a></footer></main></div>
<script src="assets/vendor/plotly.min.js"></script><script src="assets/vendor/katex/katex.min.js"></script><script src="assets/vendor/katex/contrib/auto-render.min.js"></script>{''.join('<script src="'+s+'"></script>' for s in scripts)}<script src="assets/workshop.js"></script></body></html>'''
    (ROOT / (slug+'.html')).write_text(document)
(ROOT/'.nojekyll').write_text('')
print(f'Built {len(PAGES)} pages and {len(scripts)} reproducible result shards.')
