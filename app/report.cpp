#include "report.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace vmd::cli {
namespace {

struct Point { double x; double y; };

std::string html_escape(const std::string& value) {
    std::string result;
    for (char character : value) {
        switch (character) {
        case '&': result += "&amp;"; break;
        case '<': result += "&lt;"; break;
        case '>': result += "&gt;"; break;
        case '"': result += "&quot;"; break;
        case '\'': result += "&#39;"; break;
        default: result += character;
        }
    }
    return result;
}

std::vector<Point> envelope(const std::vector<Point>& input, std::size_t max_intervals = 2000) {
    if (input.size() <= max_intervals) return input;
    const std::size_t buckets = std::max<std::size_t>(1, max_intervals / 2);
    const double width = static_cast<double>(input.size()) / static_cast<double>(buckets);
    std::vector<Point> output;
    output.reserve(buckets * 2);
    for (std::size_t bucket = 0; bucket < buckets; ++bucket) {
        const std::size_t begin = static_cast<std::size_t>(std::floor(static_cast<double>(bucket) * width));
        const std::size_t end = std::min(input.size(), static_cast<std::size_t>(std::ceil(static_cast<double>(bucket + 1) * width)));
        if (begin >= end) continue;
        std::size_t minimum = begin;
        std::size_t maximum = begin;
        for (std::size_t i = begin + 1; i < end; ++i) {
            if (input[i].y < input[minimum].y) minimum = i;
            if (input[i].y > input[maximum].y) maximum = i;
        }
        if (minimum <= maximum) {
            output.push_back(input[minimum]);
            if (maximum != minimum) output.push_back(input[maximum]);
        } else {
            output.push_back(input[maximum]);
            output.push_back(input[minimum]);
        }
    }
    return output;
}

void write_points(std::ostream& output, const std::vector<Point>& points) {
    output << '[';
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (i != 0) output << ',';
        output << '[' << points[i].x << ',' << points[i].y << ']';
    }
    output << ']';
}

std::vector<Point> vector_series(const Eigen::VectorXd& values, double sample_rate_hz) {
    std::vector<Point> points(static_cast<std::size_t>(values.size()));
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        points[static_cast<std::size_t>(i)] = {static_cast<double>(i) / sample_rate_hz, values(i)};
    }
    return envelope(points);
}

std::vector<Point> mode_series(const Result& result, Eigen::Index mode, double sample_rate_hz) {
    std::vector<Point> points(static_cast<std::size_t>(result.modes.cols()));
    for (Eigen::Index i = 0; i < result.modes.cols(); ++i) {
        points[static_cast<std::size_t>(i)] = {static_cast<double>(i) / sample_rate_hz, result.modes(mode, i)};
    }
    return envelope(points);
}

std::vector<Point> spectrum_series(const Result& result, Eigen::Index mode, double sample_rate_hz) {
    const Eigen::Index n = result.spectra.rows();
    const Eigen::Index center = n / 2;
    std::vector<Point> points;
    points.reserve(static_cast<std::size_t>(n - center));
    for (Eigen::Index i = center; i < n; ++i) {
        points.push_back({static_cast<double>(i - center) * sample_rate_hz / static_cast<double>(n),
                          std::abs(result.spectra(i, mode))});
    }
    return envelope(points);
}

std::vector<Point> omega_series(const Result& result, Eigen::Index mode) {
    std::vector<Point> points(static_cast<std::size_t>(result.omega_hz.rows()));
    for (Eigen::Index i = 0; i < result.omega_hz.rows(); ++i) {
        points[static_cast<std::size_t>(i)] = {static_cast<double>(i), result.omega_hz(i, mode)};
    }
    return envelope(points);
}

const char* color(Eigen::Index index) {
    static constexpr const char* palette[] = {
        "#1769aa", "#c43c39", "#27845f", "#8a4fa3", "#d17b0f",
        "#007f86", "#b13c75", "#5d6d2c", "#6b5fca", "#8b5a2b"
    };
    return palette[static_cast<std::size_t>(index) % (sizeof(palette) / sizeof(palette[0]))];
}

}  // namespace

void write_html_report(
    const std::filesystem::path& path,
    const InputSignal& input,
    const Options& options,
    const Result& result) {
    std::ofstream output(path, std::ios::binary);
    if (!output) throw IoError("cannot open report file: " + path.string());
    output << std::setprecision(12);
    const Eigen::VectorXd reconstruction = result.modes.colwise().sum().transpose();
    const Eigen::VectorXd residual = input.values - reconstruction;

    output << R"HTML(<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>VMD decomposition report</title>
<style>
:root{color-scheme:light;--ink:#18202a;--muted:#5d6875;--line:#d8dee6;--paper:#fff;--wash:#f4f7f9;--ok:#1c7253;--warn:#9b4b19;font-family:Inter,ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif}*{box-sizing:border-box}body{margin:0;background:var(--wash);color:var(--ink)}header{background:#10212b;color:#fff;padding:28px max(20px,calc((100vw - 1180px)/2)) 24px}h1{font-size:clamp(24px,4vw,38px);margin:0 0 8px;letter-spacing:0}header p{margin:0;color:#c7d7df;overflow-wrap:anywhere}.topline{display:flex;align-items:center;gap:12px;flex-wrap:wrap}.status{font-size:13px;font-weight:700;padding:4px 8px;border-radius:4px;background:#dff5ea;color:var(--ok)}.status.warn{background:#ffe9dc;color:var(--warn)}main{width:min(1180px,calc(100% - 32px));margin:22px auto 48px}.toolbar{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:18px}.button{display:inline-flex;align-items:center;gap:6px;padding:8px 11px;border:1px solid #b9c4ce;border-radius:5px;background:#fff;color:#24313a;text-decoration:none;font-size:14px;cursor:pointer}.button:hover{border-color:#657785;background:#f8fafb}.summary{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:1px;background:var(--line);border:1px solid var(--line);margin-bottom:20px}.metric{background:var(--paper);padding:14px}.metric b{display:block;font-size:19px;line-height:1.2}.metric span{font-size:12px;color:var(--muted)}section{background:var(--paper);border:1px solid var(--line);padding:18px;margin-top:14px}h2{font-size:18px;margin:0 0 4px}h3{font-size:14px;margin:0}.sub{margin:0 0 14px;color:var(--muted);font-size:13px}.chart{width:100%;height:auto;aspect-ratio:3.25/1;min-height:190px;display:block;background:#fff}.chart text{fill:#63707c;font-size:11px}.chart .grid{stroke:#e7ebef}.chart .axis{stroke:#9aa6b1}.chart path.series{fill:none;stroke-width:1.6;vector-effect:non-scaling-stroke}.chart .cursor{stroke:#5f6c78;stroke-dasharray:4 4;display:none}.readout{min-height:20px;font:12px ui-monospace,SFMono-Regular,Consolas,monospace;color:#34424c;margin-top:5px}.chart-actions{display:flex;justify-content:flex-end;margin-top:5px}.modes{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}.mode-panel{border-top:3px solid var(--mode-color);padding-top:9px;min-width:0}.mode-panel.hidden{display:none}.mode-panel .chart{aspect-ratio:3/1;min-height:155px}.toggles{display:flex;gap:8px 14px;flex-wrap:wrap;margin:10px 0 14px}.toggle{display:inline-flex;align-items:center;gap:6px;font-size:13px}.toggle input{accent-color:var(--mode-color)}footer{color:var(--muted);font-size:12px;margin-top:18px;text-align:center}@media(max-width:760px){header{padding:22px 16px}main{width:min(100% - 20px,1180px);margin-top:12px}.summary{grid-template-columns:repeat(2,minmax(0,1fr))}section{padding:13px}.modes{grid-template-columns:1fr}.chart{aspect-ratio:1.75/1;min-height:220px}.mode-panel .chart{aspect-ratio:1.8/1;min-height:200px}.toolbar{gap:7px}.button{padding:7px 9px}}
</style></head><body><header><div class="topline"><h1>VMD decomposition report</h1>)HTML";
    output << "<span class=\"status" << (result.diagnostics.converged ? "" : " warn") << "\">"
           << (result.diagnostics.converged ? "Converged" : "Not converged") << "</span></div><p>"
           << html_escape(input.source) << " · column " << html_escape(input.column_name) << " · "
           << input.values.size() << " samples at " << options.sample_rate_hz << " Hz</p></header>";
    output << R"HTML(<main><nav class="toolbar" aria-label="Report downloads"><a class="button" href="modes.csv" download>↓ Modes CSV</a><a class="button" href="spectra.csv" download>↓ Spectra CSV</a><a class="button" href="omega.csv" download>↓ Omega CSV</a><a class="button" href="summary.json" download>↓ Summary JSON</a></nav><div class="summary">)HTML";
    output << "<div class=\"metric\"><b>" << result.modes.rows() << "</b><span>modes</span></div>"
           << "<div class=\"metric\"><b>" << result.diagnostics.iterations << "</b><span>iterations</span></div>"
           << "<div class=\"metric\"><b>" << result.diagnostics.reconstruction_error << "</b><span>relative reconstruction error</span></div>"
           << "<div class=\"metric\"><b>" << result.timings.iteration_ms << " ms</b><span>main iteration</span></div></div>";
    output << R"HTML(<section><h2>Input and reconstruction</h2><p class="sub">Time in seconds. Hover for sample values.</p><div class="chart-shell" data-chart="overview"></div></section><section><h2>Reconstruction residual</h2><p class="sub">Input minus the sum of all modes.</p><div class="chart-shell" data-chart="residual"></div></section><section><h2>Intrinsic mode functions</h2><p class="sub">Mode controls also update spectrum and convergence views.</p><div class="toggles" id="mode-toggles"></div><div class="modes" id="mode-grid"></div></section><section><h2>Mode spectra</h2><p class="sub">One-sided magnitude spectrum in hertz.</p><div class="chart-shell" data-chart="spectra"></div></section><section><h2>Center-frequency convergence</h2><p class="sub">Estimated center frequency by iteration, in hertz.</p><div class="chart-shell" data-chart="omega"></div></section><footer>Generated by VMD_cpp. Plot data uses a min/max envelope; CSV files retain every sample.</footer></main><script>const DATA={overview:[)HTML";
    output << "{name:'Input',color:'#18202a',points:"; write_points(output, vector_series(input.values, options.sample_rate_hz)); output << "},";
    output << "{name:'Reconstruction',color:'#1769aa',points:"; write_points(output, vector_series(reconstruction, options.sample_rate_hz)); output << "}],residual:[{name:'Residual',color:'#c43c39',points:"; write_points(output, vector_series(residual, options.sample_rate_hz)); output << "}],modes:[";
    for (Eigen::Index k = 0; k < result.modes.rows(); ++k) { if (k != 0) output << ','; output << "{name:'IMF " << (k + 1) << "',color:'" << color(k) << "',points:"; write_points(output, mode_series(result, k, options.sample_rate_hz)); output << '}'; }
    output << "],spectra:[";
    for (Eigen::Index k = 0; k < result.modes.rows(); ++k) { if (k != 0) output << ','; output << "{name:'IMF " << (k + 1) << "',color:'" << color(k) << "',points:"; write_points(output, spectrum_series(result, k, options.sample_rate_hz)); output << '}'; }
    output << "],omega:[";
    for (Eigen::Index k = 0; k < result.modes.rows(); ++k) { if (k != 0) output << ','; output << "{name:'IMF " << (k + 1) << "',color:'" << color(k) << "',points:"; write_points(output, omega_series(result, k)); output << '}'; }
    output << R"HTML(]};
const visible=DATA.modes.map(function(){return true}),NS='http://www.w3.org/2000/svg';function el(name,attrs){const node=document.createElementNS(NS,name);Object.entries(attrs||{}).forEach(function(entry){node.setAttribute(entry[0],entry[1])});return node}function fmt(value){return Number.isFinite(value)?((Math.abs(value)>=1000||(Math.abs(value)<.001&&value!==0))?value.toExponential(3):value.toFixed(4)):'—'}
function renderChart(shell,series,title){shell.textContent='';const active=series.filter(function(_,index){return visible[index]===undefined||visible[index]});const svg=el('svg',{class:'chart',viewBox:'0 0 900 280',role:'img','aria-label':title}),margin={l:66,r:20,t:18,b:40},w=814,h=222,points=active.flatMap(function(item){return item.points});let xmin=points.length?Math.min.apply(null,points.map(function(p){return p[0]})):0,xmax=points.length?Math.max.apply(null,points.map(function(p){return p[0]})):1,ymin=points.length?Math.min.apply(null,points.map(function(p){return p[1]})):-1,ymax=points.length?Math.max.apply(null,points.map(function(p){return p[1]})):1;if(xmin===xmax)xmax=xmin+1;if(ymin===ymax){ymin-=1;ymax+=1}const pad=(ymax-ymin)*.06;ymin-=pad;ymax+=pad;const sx=function(x){return margin.l+(x-xmin)/(xmax-xmin)*w},sy=function(y){return margin.t+h-(y-ymin)/(ymax-ymin)*h};for(let i=0;i<=4;i++){const y=margin.t+h*i/4;svg.append(el('line',{class:'grid',x1:margin.l,x2:margin.l+w,y1:y,y2:y}));const label=el('text',{x:margin.l-9,y:y+4,'text-anchor':'end'});label.textContent=fmt(ymax-(ymax-ymin)*i/4);svg.append(label)}for(let i=0;i<=5;i++){const x=margin.l+w*i/5;svg.append(el('line',{class:'grid',x1:x,x2:x,y1:margin.t,y2:margin.t+h}));const label=el('text',{x:x,y:margin.t+h+24,'text-anchor':i===0?'start':(i===5?'end':'middle')});label.textContent=fmt(xmin+(xmax-xmin)*i/5);svg.append(label)}svg.append(el('line',{class:'axis',x1:margin.l,x2:margin.l+w,y1:margin.t+h,y2:margin.t+h}));svg.append(el('line',{class:'axis',x1:margin.l,x2:margin.l,y1:margin.t,y2:margin.t+h}));active.forEach(function(item){const d=item.points.map(function(p,i){return(i?'L':'M')+sx(p[0]).toFixed(2)+','+sy(p[1]).toFixed(2)}).join(' ');svg.append(el('path',{class:'series',d:d,stroke:item.color,'data-name':item.name}))});const cursor=el('line',{class:'cursor',x1:0,x2:0,y1:margin.t,y2:margin.t+h});svg.append(cursor);const readout=document.createElement('div');readout.className='readout';readout.textContent='Move over the plot to inspect values.';svg.addEventListener('pointermove',function(event){const box=svg.getBoundingClientRect(),vx=(event.clientX-box.left)/box.width*900,x=xmin+Math.max(0,Math.min(1,(vx-margin.l)/w))*(xmax-xmin);cursor.style.display='block';cursor.setAttribute('x1',sx(x));cursor.setAttribute('x2',sx(x));readout.textContent='x '+fmt(x)+' · '+active.map(function(item){let best=item.points[0];item.points.forEach(function(p){if(Math.abs(p[0]-x)<Math.abs(best[0]-x))best=p});return item.name+' '+fmt(best[1])}).join(' · ')});svg.addEventListener('pointerleave',function(){cursor.style.display='none'});const actions=document.createElement('div');actions.className='chart-actions';const button=document.createElement('button');button.className='button';button.type='button';button.textContent='↓ SVG';button.title='Download this plot as SVG';button.onclick=function(){const source=new XMLSerializer().serializeToString(svg),blob=new Blob([source],{type:'image/svg+xml'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=title.toLowerCase().replace(/[^a-z0-9]+/g,'-')+'.svg';a.click();URL.revokeObjectURL(url)};actions.append(button);shell.append(svg,readout,actions)}
function renderAll(){renderChart(document.querySelector('[data-chart=overview]'),DATA.overview,'Input and reconstruction');renderChart(document.querySelector('[data-chart=residual]'),DATA.residual,'Reconstruction residual');renderChart(document.querySelector('[data-chart=spectra]'),DATA.spectra,'Mode spectra');renderChart(document.querySelector('[data-chart=omega]'),DATA.omega,'Center frequency convergence');document.querySelectorAll('.mode-panel').forEach(function(panel,index){panel.classList.toggle('hidden',!visible[index])})}const toggles=document.getElementById('mode-toggles'),grid=document.getElementById('mode-grid');DATA.modes.forEach(function(mode,index){const label=document.createElement('label');label.className='toggle';label.style.setProperty('--mode-color',mode.color);const checkbox=document.createElement('input');checkbox.type='checkbox';checkbox.checked=true;checkbox.onchange=function(){visible[index]=checkbox.checked;renderAll()};label.append(checkbox,document.createTextNode(mode.name));toggles.append(label);const panel=document.createElement('article');panel.className='mode-panel';panel.style.setProperty('--mode-color',mode.color);const heading=document.createElement('h3');heading.textContent=mode.name;const shell=document.createElement('div');shell.className='chart-shell';panel.append(heading,shell);grid.append(panel);renderChart(shell,[mode],mode.name)});renderAll();
</script></body></html>)HTML";
}

}  // namespace vmd::cli
