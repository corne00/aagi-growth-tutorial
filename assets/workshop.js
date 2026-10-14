/* All scientific values are saved outputs of analysis/export.R. */
(() => {
  'use strict';
  const D = window.WORKSHOP_DATA;
  const page = document.body.dataset.page;
  const $ = id => document.getElementById(id);
  if (!D || !D.trajectories) {
    const p = document.createElement('p'); p.className = 'error';
    p.textContent = 'Analysis results could not be loaded. Reload the page or run the reproducible export workflow.';
    $('main').prepend(p); return;
  }
  Object.values(D.trajectories).forEach(c => { c.weather = D.weather[c.weather_key]; });
  const colors = ['#235f83','#b76435','#518365','#875b88'];
  const units = {temperature:'Mean temperature (°C)',Radn:'Radiation (MJ/m²/day)',Rain:'Rainfall (mm/day)',
    thermal_time:'Thermal time (°C days)',rain14:'14-day rainfall (mm)'};
  const num = (x, digits=3) => Number(x).toFixed(digits);
  const trace = (x,y,name,color=colors[0],extra={}) => ({x,y,name,type:'scatter',mode:'lines',line:{color,width:2},...extra});
  const config = {responsive:true,displayModeBar:false};
  function plot(id,traces,layout={}) {
    if (!$(id)) return;
    const base={paper_bgcolor:'#fff',plot_bgcolor:'#fff',height:360,
      font:{family:'Arial, sans-serif',color:'#34434b',size:12},
      margin:{l:68,r:18,t:24,b:52},hovermode:'x unified',
      legend:{orientation:'h',y:-.2,x:0,font:{size:11}},
      xaxis:{title:{text:'Days after sowing',font:{size:12}},showgrid:false,zeroline:false,linecolor:'#ccd4d9'},
      yaxis:{title:{text:'Biomass (t/ha)',font:{size:12}},gridcolor:'#edf0f2',zeroline:false},
      ...layout};
    return Plotly.react($(id),traces,base,config);
  }
  function bands(x,s,color,name,extra={}) {
    const rgb = color === colors[0] ? '35,95,131' : color === colors[1] ? '183,100,53' : '81,131,101';
    return [trace(x,s.lower,'',color,{line:{width:0},showlegend:false,hoverinfo:'skip',...extra}),
      trace(x,s.upper,'',color,{line:{width:0},fill:'tonexty',fillcolor:`rgba(${rgb},0.14)`,showlegend:false,hoverinfo:'skip',...extra}),
      trace(x,s.median,name,color,extra)];
  }
  function fill(id,values,defaultValue) {
    if (!$(id)) return;
    $(id).replaceChildren(...values.map(v=>{
      const option=document.createElement('option'); option.value=Array.isArray(v)?v[0]:v;
      option.textContent=Array.isArray(v)?v[1]:v; return option;
    }));
    if (defaultValue !== undefined) $(id).value=String(defaultValue);
  }
  function controls() {
    fill('genotype',D.genotypes,'g008');fill('genotype2',D.genotypes,'g052');fill('genotype3',D.genotypes,'g185');
    fill('location',D.locations,'Yanco'); fill('year',D.years,2010);
  }
  const val = (id,fallback) => $(id)?.value || fallback;
  const curve = (model='hierarchical',location=val('location','Yanco'),year=val('year',2010),genotype=val('genotype','g008')) =>
    D.trajectories[[model,location,year,genotype].join('|')];
  const table = (id,heads,rows) => {
    const container=$(id); if(!container) return;
    const wrapper=document.createElement('div');wrapper.className='table-scroll';
    const table=document.createElement('table');const head=document.createElement('thead');const hr=document.createElement('tr');
    heads.forEach(h=>{const cell=document.createElement('th');cell.textContent=h;hr.append(cell);});head.append(hr);table.append(head);
    const body=document.createElement('tbody');rows.forEach(row=>{const tr=document.createElement('tr');row.forEach(v=>{
      const td=document.createElement('td');td.textContent=String(v);tr.append(td);});body.append(tr);});table.append(body);
    wrapper.append(table);container.replaceChildren(wrapper);
  };
  function diagnostics(id='diagnostic-table',names=Object.keys(D.diagnostics)) {
    table(id,['Fit','Maximum R̂','Minimum bulk ESS','Minimum tail ESS','Check'], names.map(name=>{
      const d=D.diagnostics[name];return [name.replaceAll('_',' '),num(d.max_rhat,4),Math.round(d.min_bulk_ess),Math.round(d.min_tail_ess),d.passed?'Passed':'Needs longer sampling'];
    }));
  }
  function metricRows(metrics,label='') {return metrics.map(m=>[label || m.model,m.fold || '',m.scope,num(m.rmse),num(m.r2),num(m.correlation)]);}
  function metrics(id,rows){table(id,['Model','Held-out years','Response','RMSE (t/ha)','R²','Correlation'],rows);}
  function interval(c,level=.95,kind='credible') {return c.intervals.find(i=>Math.abs(i.level-level)<1e-8)[kind];}
  function gridLayout(rows,titles,height=660) {
    const layout={height,grid:{rows,columns:1,pattern:'independent',roworder:'top to bottom'},
      margin:{l:78,r:20,t:26,b:62},legend:{orientation:'h',y:-.09,x:0}};
    for(let i=1;i<=rows;i++){
      const suffix=i===1?'':i;
      layout['xaxis'+suffix]={title:{text:i===rows?'Days after sowing':''},showgrid:false,zeroline:false};
      layout['yaxis'+suffix]={title:{text:titles[i-1],font:{size:12}},gridcolor:'#edf0f2',zeroline:false};
    }return layout;
  }
  function flowPlot(id,c,time) {
    const w=c.weather;const traces=[trace(w.time,w.temperature,'Temperature',colors[1],{xaxis:'x',yaxis:'y',showlegend:false}),
      trace(c.rate.time,c.rate.mean,'Mean predicted growth rate',colors[2],{xaxis:'x2',yaxis:'y2',showlegend:false}),
      trace(c.daily_time,c.daily_biomass,'APSIM', '#6a7378',{xaxis:'x3',yaxis:'y3',line:{color:'#6a7378',width:1.5}}),
      trace(c.time,c.mean_biomass,'Mean predicted biomass',colors[0],{xaxis:'x3',yaxis:'y3'})];
    const layout=gridLayout(3,['Temperature (°C)','Growth (t/ha/day)','Biomass (t/ha)']);
    for(const axis of ['xaxis','xaxis2','xaxis3'])layout[axis].range=[0,Math.max(...c.daily_time)];
    if(time!==undefined){
      layout.shapes=[1,2,3].map(i=>({type:'line',xref:i===1?'x':'x'+i,yref:i===1?'y domain':'y'+i+' domain',
        x0:time,x1:time,y0:0,y1:1,line:{color:'#7b858b',width:1,dash:'dot'}}));
      const rates=c.rate.time.filter(t=>t<=time);
      traces.splice(1,0,trace(rates,c.rate.mean.slice(0,rates.length),'Accumulated rate area',colors[2],{
        xaxis:'x2',yaxis:'y2',fill:'tozeroy',fillcolor:'rgba(81,131,101,.12)',showlegend:false,line:{width:0},hoverinfo:'skip'}));
    }plot(id,traces,layout);
  }
  function introduction() {
    plot('endpoints',D.introduction.map((c,i)=>trace(c.time,c.biomass,`${c.id}: ${num(c.endpoint,2)} t/ha`,colors[i])),
      {height:340,xaxis:{title:{text:'Time'},showgrid:false,zeroline:false},yaxis:{title:{text:'Biomass'},showgrid:false,zeroline:false},legend:{y:-.22}});
    const c=curve('hierarchical','Yanco',2010,'g008');
    $('time').max=Math.max(...c.daily_time);flowPlot('flow',c,+val('time',70));
  }
  function dataset() {
    $('dataset-counts').textContent=`${D.counts.records.toLocaleString()} daily records · ${D.counts.trajectories.toLocaleString()} trajectories · biomass in t/ha.`;
    const chosen=['genotype','genotype2','genotype3'].map(id=>val(id));
    const traces=[];
    D.locations.forEach((location,i)=>{
      plot('data-'+location,chosen.map((g,j)=>{
        const c=curve('hierarchical',location,val('year'),g);return trace(c.daily_time,c.daily_biomass,g,colors[j],{showlegend:true});
      }),{height:280,yaxis:{title:{text:'Biomass (t/ha)'},range:[0,18],gridcolor:'#edf0f2',zeroline:false},
        xaxis:{title:{text:'Days after sowing'},range:[0,190],showgrid:false},margin:{l:55,r:10,t:10,b:55},legend:{y:-.28}});
      const c=curve('hierarchical',location,val('year'),chosen[0]);const driver=val('driver','temperature');
      traces.push(trace(c.weather.time,c.weather[driver],location,colors[i]));
    });
    const driver=val('driver','temperature');plot('weather',traces,{yaxis:{title:{text:units[driver]},gridcolor:'#edf0f2',zeroline:false},height:320});
    const c=curve('hierarchical','Yanco',val('year'),chosen[0]);
    table('table-preview',['Table','Identifier','Time','Example value'],[
      ['Observations',`${c.location}_${c.year} / ${c.genotype}`,c.daily_time[30],num(c.daily_biomass[30])+' t/ha'],
      ['Drivers',`${c.location}_${c.year}`,c.weather.time[30],num(c.weather.temperature[30],1)+' °C'],
      ['Metadata',`${c.location}_${c.year} / ${c.genotype}`,'—',`genotype=${c.genotype}; location=${c.location}; year=${c.year}`]]);
  }
  function response() {
    const b=D.basis; const j=+val('feature',0);const name=b.features[j];
    plot('basis',[
      trace(b.time,b.values[j],name,colors[1],{xaxis:'x',yaxis:'y',showlegend:false}),
      trace(b.time,b.contributions[j],'Selected fitted contribution',colors[1],{xaxis:'x2',yaxis:'y2'}),
      trace(b.time,b.rate,'Complete population response',colors[0],{xaxis:'x2',yaxis:'y2'})],
      gridLayout(2,['Scaled feature','Growth rate (t/ha/day)'],520));
  }
  function hierarchy() {
    const h=D.hierarchy;const g=val('genotype','g008');const selected=h.genotypes[g];
    plot('hierarchy',[
      trace(h.time,h.shared,'Population response', '#68747c',{xaxis:'x',yaxis:'y'}),
      ...bands(h.time,selected.response,colors[0],g,{xaxis:'x',yaxis:'y'}),
      ...bands(h.time,selected.deviation,colors[1],'Genotype deviation',{xaxis:'x2',yaxis:'y2'})],
      gridLayout(2,['Growth rate (t/ha/day)','Deviation (t/ha/day)'],520));
    const density=val('density','one');const names=['independent','hierarchical'];
    const scope='terminal biomass';
    const records=names.map(m=>D.pooling[`pooling_${density}_${m}`]);
    plot('pooling',[{type:'bar',x:['Independent','Hierarchical'],y:records.map(r=>r.metrics.find(m=>m.scope===scope).rmse),
      marker:{color:[colors[1],colors[0]]},hovertemplate:'%{x}<br>RMSE %{y:.3f} t/ha<extra></extra>'}],
      {height:300,xaxis:{title:{text:''},showgrid:false},yaxis:{title:{text:'Terminal RMSE (t/ha)'},gridcolor:'#edf0f2'}});
    metrics('location-results',metricRows(D.location.metrics,'Location extension').map(r=>{r[1]=D.location.fold;return r;}));
  }
  function design() {
    const d=D.design;const i=+val('interval',4);const c=curve('hierarchical','Yanco',2010,'g008');
    $('interval').max=d.time0.length-1;
    plot('interval-growth',[trace(c.daily_time,c.daily_biomass,'APSIM',colors[0]),
      {x:[d.time0[i],d.time1[i]],y:[c.observed[i],c.observed[i+1]],name:`ΔB = ${num(d.dy[i])} t/ha`,
       type:'scatter',mode:'markers',marker:{size:8,color:colors[1]}}],
      {height:300,shapes:[{type:'rect',xref:'x',yref:'paper',x0:d.time0[i],x1:d.time1[i],y0:0,y1:1,
        fillcolor:'rgba(183,100,53,.10)',line:{width:0}}]});
    plot('matrix',[{type:'heatmap',z:d.X,x:d.columns,y:d.time1,colorscale:[[0,'#b76435'],[.5,'#fff'],[1,'#235f83']],
      zmid:0,colorbar:{title:{text:'Integrated feature'},thickness:12},hovertemplate:'%{x}<br>Interval ending day %{y}<br>%{z:.3f}<extra></extra>'}],
      {height:380,margin:{l:65,r:80,t:15,b:120},xaxis:{tickangle:-40,title:{text:''}},
       yaxis:{title:{text:'Interval end (DAS)'},autorange:'reversed'},
       shapes:[{type:'rect',xref:'paper',yref:'y',x0:0,x1:1,y0:d.time1[i]-5,y1:d.time1[i]+5,
        line:{color:'#202b32',width:1},fillcolor:'rgba(0,0,0,0)'}]});
  }
  function gibbs() {
    const parameter=val('parameter','beta');const maximum=D.traces[0][parameter].length;
    $('iteration').max=maximum;if(+$('iteration').value>maximum)$('iteration').value=maximum;
    const count=+val('iteration',maximum);
    const x=Array.from({length:count},(_,i)=>i+1);
    plot('traces',D.traces.map((c,i)=>trace(x,c[parameter].slice(0,count),`Chain ${c.chain}`,colors[i],{line:{color:colors[i],width:1}})),
      {height:330,xaxis:{title:{text:'Retained iteration'},showgrid:false},yaxis:{title:{text:parameter==='sigma'?'Residual SD (t/ha)':parameter==='precision'?'Genotype precision':'Population coefficient'},gridcolor:'#edf0f2'}});
    diagnostics();
  }
  function prediction() {
    const model=val('model','hierarchical');const c=curve(model);
    $('prediction-split').textContent=`Held out: ${c.fold} · ${c.location} ${c.year} · ${c.genotype}. All observations in this year block were excluded from fitting.`;
    flowPlot('prediction',c);
    metrics('metrics-table',metricRows(D.metrics));
  }
  function uncertainty() {
    const c=curve();const kind=val('kind','credible');const level=+val('level',.95);
    const s=interval(c,level,kind);
    plot('uncertainty',[
      ...(kind==='credible'?c.draws.map(draw=>trace(c.time,draw,'',colors[0],{opacity:.18,showlegend:false,line:{color:colors[0],width:.7},hoverinfo:'skip'})):[]),
      ...bands(c.time,s,colors[0],`${Math.round(level*100)}% ${kind} interval`),
      trace(c.daily_time,c.daily_biomass,'APSIM','#343e44',{line:{color:'#343e44',width:1.5}})],{height:360});
    plot('coverage',[trace([0,1],[0,1],'Nominal','#7b858b',{line:{color:'#7b858b',width:1,dash:'dash'}}),
      ...Object.values(D.coverage).map((c,i)=>trace(c.nominal,c[kind],c.fold,colors[i],{mode:'lines+markers'}))],
      {height:350,xaxis:{title:{text:'Nominal coverage'},range:[0,1],tickformat:'.0%',showgrid:false},
       yaxis:{title:{text:'Empirical terminal coverage'},range:[0,1],tickformat:'.0%',gridcolor:'#edf0f2'}});
  }
  function interpretation() {
    const driver=val('driver','temperature');const g=val('ale-genotype','population');const bins=val('bins',10);
    const a=D.ale[[g,driver,bins].join('|')];
    const traces=[...bands(a.x,a,colors[0],g,{xaxis:'x',yaxis:'y'}),
      {type:'bar',x:a.x,y:a.counts,name:'Reference days',xaxis:'x2',yaxis:'y2',marker:{color:'#a7b6bf'},showlegend:false}];
    const layout=gridLayout(2,['ALE (t/ha/day)','Reference days'],490);
    layout.xaxis.title={text:units[driver]};layout.xaxis2.title={text:units[driver]};
    layout.yaxis.domain=[.36,1];layout.yaxis2.domain=[0,.22];
    layout.shapes=[{type:'line',xref:'paper',yref:'y',x0:0,x1:1,y0:0,y1:0,line:{width:1,color:'#b5bec4',dash:'dot'}}];
    plot('ale',traces,layout);
  }
  function gxe() {
    const first=val('genotype','g008');const second=val('genotype2','g052');
    const traces=[];D.locations.forEach((l,i)=>{
      const a=curve('hierarchical',l,val('year'),first),b=curve('hierarchical',l,val('year'),second);
      const last=Math.min(Math.max(...a.time),Math.max(...b.time));
      const common=a.time.filter(t=>t<=last && b.time.includes(t));
      const aMean=a.mean_biomass,bMean=b.mean_biomass;
      traces.push(trace(common,common.map(t=>aMean[a.time.indexOf(t)]-bMean[b.time.indexOf(t)]),l,colors[i]));
      traces.push(trace(common,common.map(t=>a.observed[a.time.indexOf(t)]-b.observed[b.time.indexOf(t)]),l+' APSIM',colors[i],
        {line:{color:colors[i],width:1,dash:'dot'},showlegend:false}));
    });
    plot('contrasts',traces,{height:340,yaxis:{title:{text:`${first} − ${second} (t/ha)`},gridcolor:'#edf0f2',zeroline:true,zerolinecolor:'#b8c2c9'}});
    const c=D.gxe[val('time',84)];
    plot('gxe-components',[{type:'bar',name:'APSIM',x:['G','E','G×E'],y:c.observed,marker:{color:'#889aa6'},xaxis:'x',yaxis:'y'},
      {type:'bar',name:'Predicted',x:['G','E','G×E'],y:c.predicted,marker:{color:colors[0]},xaxis:'x',yaxis:'y'},
      {type:'bar',name:'Correlation',x:['G','E','G×E'],y:c.correlation,marker:{color:colors[2]},xaxis:'x2',yaxis:'y2',showlegend:false}],
      {height:340,grid:{rows:1,columns:2,pattern:'independent'},margin:{l:65,r:25,t:20,b:65},barmode:'group',
       xaxis:{title:{text:'Component'}},xaxis2:{title:{text:'Component'}},yaxis:{title:{text:'Mean-square fraction'},range:[0,1],tickformat:'.0%',gridcolor:'#edf0f2'},
       yaxis2:{title:{text:'Observed–predicted correlation'},range:[-1,1],gridcolor:'#edf0f2'}});
  }
  function reduction() {
    const rank=+val('rank',8);const energy=D.svd.full_energy;
    $('rank').max=energy.length;
    plot('spectrum',[trace(energy.map((_,i)=>i+1),energy,'Training spectral energy',colors[0],{mode:'lines+markers'})],
      {height:330,xaxis:{title:{text:'Retained nonconstant feature rank'},showgrid:false},yaxis:{title:{text:'Cumulative squared singular-value energy'},range:[0,1.02],tickformat:'.0%',gridcolor:'#edf0f2'},
       shapes:[{type:'line',xref:'x',yref:'paper',x0:rank,x1:rank,y0:0,y1:1,line:{color:colors[1],width:1,dash:'dot'}}]});
    const baseline=D.metrics.filter(m=>m.model==='hierarchical' && m.fold===D.svd.fold);
    metrics('svd-results',[...metricRows(baseline,'Original features'),...metricRows(D.svd.metrics,`SVD: ${D.svd.features} total features`).map(r=>{r[1]=D.svd.fold;return r;})]);
  }
  function reproduce() {
    const s=D.specification;
    table('run-manifest',['Item','Value'],[
      ['Full-data records',D.counts.records.toLocaleString()],['Regression intervals',D.counts.intervals.toLocaleString()],
      ['Main response library',s.features.length+' features'],['Observation spacing',s.spacing_days+' days, plus native final interval'],
      ['Quadrature',s.integration],['Full-fit chains',s.chains],['Full-fit warmup / retained per chain',s.warmup+' / '+s.samples_per_chain],
      ['Seed',s.seed],['Source checksum',D.source.md5],['R',s.R],['data.table / jsonlite',s.packages['data.table']+' / '+s.packages.jsonlite]]);
  }
  controls();
  if(page==='response')fill('feature',D.basis.features.map((f,i)=>[i,f]),0);
  if(page==='interpretation')fill('ale-genotype',[['population','Population'],...D.genotypes], 'population');
  if($('level'))$('level').value='0.95';
  if($('iteration')){$('iteration').max=D.traces[0].beta.length;$('iteration').value=D.traces[0].beta.length;}
  const render={index:introduction,dataset,response,hierarchy,design,gibbs,prediction,uncertainty,interpretation,gxe,reduction,reproduce}[page];
  const update=()=>{document.querySelectorAll('input[type=range]').forEach(input=>{const out=$(input.id+'-value');if(out)out.value=input.value;});render?.();};
  document.querySelectorAll('select,input[type=range]').forEach(input=>input.addEventListener(input.tagName==='INPUT'?'input':'change',update));
  document.querySelectorAll('[data-function]').forEach(code=>{code.textContent=D.code[code.dataset.function] || 'See the complete analysis script in the repository.';});
  document.querySelectorAll('.copy').forEach(button=>button.addEventListener('click',async()=>{
    const code=button.parentElement.querySelector('code');
    try{await navigator.clipboard.writeText(code.textContent);button.textContent='Copied';setTimeout(()=>button.textContent='Copy',1200);}
    catch{button.textContent='Select code to copy';}
  }));
  if(window.renderMathInElement)renderMathInElement(document.querySelector('main'),{
    delimiters:[{left:'\\[',right:'\\]',display:true},{left:'\\(',right:'\\)',display:false}],throwOnError:false});
  update();
})();
