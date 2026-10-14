/* Programmatic DOM/interaction tests; no browser rendering claim. */
const fs=require('fs');const path=require('path');const assert=require('assert');
const {JSDOM}=require('jsdom');
const root=path.resolve(__dirname,'..');
const pages=fs.readdirSync(root).filter(p=>p.endsWith('.html'));
let updates=0,controls=0;
for(const page of pages){
  const html=fs.readFileSync(path.join(root,page),'utf8');
  const dom=new JSDOM(html,{runScripts:'outside-only',url:'https://corne00.github.io/aagi-growth-tutorial/'+page});
  const w=dom.window;const errors=[];w.addEventListener('error',e=>errors.push(e.error||e.message));
  w.Plotly={react(element,traces,layout){
    assert(element,'Plot element exists');
    for(const t of traces){
      if(t.type==='heatmap'){assert(Array.isArray(t.z));assert(t.z.length>0);continue;}
      assert(Array.isArray(t.x),page+': x array');assert(Array.isArray(t.y),page+': y array');
      assert.equal(t.x.length,t.y.length,page+': matching coordinates');
      assert(t.y.every(Number.isFinite),page+': finite values');
    }updates++;return Promise.resolve();
  }};
  Object.defineProperty(w.navigator,'clipboard',{value:{writeText:async()=>{}}});
  for(const script of [...w.document.querySelectorAll('script[src]')]){
    const src=script.getAttribute('src');
    if(src.includes('plotly.min.js'))continue;
    assert(fs.existsSync(path.join(root,src)),page+': script exists '+src);
    w.eval(fs.readFileSync(path.join(root,src),'utf8'));
  }
  const plotCount=w.document.querySelectorAll('.plot').length;
  for(const control of [...w.document.querySelectorAll('select,input[type=range]')]){
    if(control.tagName==='SELECT'){
      if(control.options.length>1)control.value=control.options[1].value;
      control.dispatchEvent(new w.Event('change',{bubbles:true}));
    }else{
      control.value=control.min;control.dispatchEvent(new w.Event('input',{bubbles:true}));
      control.value=control.max;control.dispatchEvent(new w.Event('input',{bubbles:true}));
    }controls++;
  }
  assert.equal(w.document.querySelectorAll('.plot').length,plotCount,'Controls do not create new plot elements');
  assert.equal(errors.length,0,page+': no JS errors: '+errors.map(String).join(';'));
  assert.equal(w.document.querySelectorAll('.katex-error').length,0,page+': valid equations');
  for(const code of w.document.querySelectorAll('[data-function]'))assert(code.textContent.length>20,'Real code inserted');
  for(const anchor of w.document.querySelectorAll('a[href]')){
    const href=anchor.getAttribute('href');if(href.startsWith('http')||href.startsWith('#'))continue;
    assert(fs.existsSync(path.join(root,href)),page+': local link '+href);
  }
  dom.window.close();
}
console.log(`${pages.length} pages, ${controls} controls and ${updates} plot updates passed.`);
