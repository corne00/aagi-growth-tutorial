"""Validate exported numerical artefacts; does not compute model results."""
import json, math
from pathlib import Path
root=Path(__file__).resolve().parents[1]
generated=root/'assets/generated'
d=json.loads((generated/'000-base.js').read_text().removeprefix('window.WORKSHOP_DATA=').rstrip(';\n'))
for p in sorted(generated.glob('*.js')):
    if p.name in ['000-base.js','traces.js']:continue
    text=p.read_text()
    if text.startswith('Object.assign(window.WORKSHOP_DATA.'):
        field=text.removeprefix('Object.assign(window.WORKSHOP_DATA.').split(',',1)[0]
        d[field].update(json.loads(text.split(',',1)[1].removesuffix(');')))
assert d['counts']['trajectories']==199*124
assert d['counts']['intervals']==278539
assert d['source']['md5']=='b208ae77ec2c9f9e18dc605d62d2812e'
assert all(x['passed'] for x in d['diagnostics'].values())
assert len(d['trajectories'])==288
for c in d['trajectories'].values():
    c['weather']=d['weather'][c['weather_key']]
    assert len(c['time'])==len(c['mean_biomass'])==len(c['observed'])
    assert c['time'][0]==0 and c['observed'][0]==0
    assert c['rate']['time']==c['daily_time']
    assert len(c['daily_time'])==len(set(c['daily_time']))
    # Mean response integration equals mean biomass at the fitted checkpoints.
    rates=c['rate']['mean']; cumulative=[0.0]
    for a,b in zip(rates,rates[1:]):cumulative.append(cumulative[-1]+(a+b)/2)
    for t,b in zip(c['time'],c['mean_biomass']):assert abs(cumulative[t]-b)<2e-5
    for interval in c['intervals']:
        for kind in ['credible','predictive']:
            s=interval[kind]
            assert len(s['lower'])==len(c['time'])
            assert all(a<=b<=z for a,b,z in zip(s['lower'],s['median'],s['upper']))
for curve in d['ale'].values():
    assert sum(curve['counts'])==6000
    assert all(n>0 for n in curve['counts'])
    assert len(curve['x'])==len(curve['median'])
for v in d['gxe'].values():
    assert abs(sum(v['observed'])-1)<1e-6
    assert abs(sum(v['predicted'])-1)<1e-6
for scope in ['trajectory checkpoints','terminal biomass']:
    metrics=[m for m in d['metrics'] if m['model']=='hierarchical' and m['scope']==scope]
    assert len(metrics)==4
    if scope=='terminal biomass':assert sum(m['n'] for m in metrics)==199*124
print('Source, convergence, integration, intervals, ALE and balanced held-out panel checks passed.')
