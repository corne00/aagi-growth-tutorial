# APSIM introduction page

This package implements only the new Introduction page. It uses a restrained documentation layout: top navigation, main article, right-hand contents list, ordinary explanatory paragraphs, an interactive four-location biomass figure and an APSIM inputs–outputs schematic.

## View locally

Open `index.html` in a modern browser. All JavaScript, equation fonts, plotting code and data are bundled. If your browser restricts local scripts, run `python -m http.server 8000` in this directory and visit `http://localhost:8000`.

## Update the existing GitHub Pages site

Copy `index.html`, `.nojekyll`, and the contents of `assets/` into your existing website repository. Merge the assets directories; do not delete the repository's other pages. The new page uses its own stylesheet and JavaScript files, so the old chapter files are not rewritten. Their redesign will happen separately, page by page.

From your repository directory:

```bash
git add .
git commit -m "Introduce APSIM biomass dataset"
git push
```

GitHub Pages will deploy the updated introduction. The new header deliberately links only to Introduction and the references on this page. Old pages remain at their original URLs but are not presented as part of the revised tutorial.

## Data

The trajectories are extracted directly from `SimulatedData124Envs_Biomass.txt` in the supplied APSIM dataset. Locations: Emerald, Merredin, Narrabri, Yanco. Available years in this introduction: 1993, 2008, 2010. Six genotypes: g008, g052, g185, g127, g009, g105. The first three are displayed by default.

Illustrative genotypes were selected by a farthest-point procedure on standardized biomass profiles at days 60, 100 and the final available day across the four locations in 2008. The selection illustrates contrast; it is not a random or representative population sample. Exact selection description is stored in the embedded data file. The plotted points themselves are unmodified daily APSIM values. No trajectories are padded or carried forward. The figure uses the same axis scales across locations.

The full dataset comprises 199 genotypes, 124 site-years, 1983–2013. It contains model-generated phenotypes, not observed field biomass. The analysis planned for later pages uses biomass and weather drivers, excluding markers and APSIM water-status outputs as predictors. No fitted model or SVD methodology is introduced on this page.

## References

Bustos-Korts et al. (2019). From QTLs to Adaptation Landscapes: Using Genotype-To-Phenotype Models to Characterize G×E Over Time. Frontiers in Plant Science 10:1540. https://doi.org/10.3389/fpls.2019.01540

Holzworth et al. (2014). APSIM — Evolution towards a new generation of agricultural systems simulation. Environmental Modelling & Software 62:327–350. https://doi.org/10.1016/j.envsoft.2014.07.009

## Checks and dependencies

Programmatic checks covered genotype and year selection, shared plot scales, finite ordered trajectory data, equations, local links and empty selections. Full visual browser review was unavailable in the creation environment; review the page in your browser before presenting.

Plotly.js 3.1.0 and KaTeX 0.16.22 are bundled with their licenses under `assets/vendor/`.

Every page is public when deployed to GitHub Pages. The noindex setting discourages search indexing and is not access control.
