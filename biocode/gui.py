"""Структурированный браузер результатов BIOCODE (блок B15), в духе Unipro UGENE.

Одностраничное десктоп-подобное приложение (http.server + inline HTML/CSS/JS, без
сборки и внешних зависимостей — как `pipeline_gui.py`). Открывает каталог прогона
(`<out>/run`) и даёт:

  • слева — НАВИГАТОР ГРУПП (аналог Project-дерева UGENE): поиск, сортировка,
    список клональных групп V+J с бейджами (размер / мутации / клады);
  • сверху — обзор прогона (метрики, топ-кандидаты по всем группам);
  • по центру — вкладки для выбранной группы (как объектные виды UGENE):
      Обзор · Дерево · Выравнивание(FR/CDR) · Мутации · Кандидаты · Клады.

Запуск:  python -m biocode gui --run EDU/results/biocode_demo/run --port 8766
"""
from __future__ import annotations

import csv
import http.server
import json
import socketserver
from pathlib import Path
from urllib.parse import parse_qs, urlparse

from .logging_ import get_logger

log = get_logger("gui")
_RUN: Path = Path(".")


def _read_tsv(path: Path, cap: int | None = None) -> list[dict]:
    if not path.is_file():
        return []
    rows = []
    with open(path, newline="") as f:
        for i, row in enumerate(csv.DictReader(f, delimiter="\t")):
            if cap and i >= cap:
                break
            rows.append(row)
    return rows


def _read_fasta(path: Path, cap: int | None = None) -> tuple[list[str], list[str]]:
    ids, seqs, cur = [], [], None
    if not path.is_file():
        return ids, seqs
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if cap and len(ids) >= cap:
                cur = None
                continue
            cur = line[1:].strip().split()[0]
            ids.append(cur)
            seqs.append("")
        elif cur is not None:
            seqs[-1] += line.strip()
    return ids, seqs


def _manifest() -> dict:
    mf = _RUN / "manifest.json"
    data = json.loads(mf.read_text()) if mf.is_file() else {"groups": [], "totals": {}}
    data["top_candidates"] = _read_tsv(_RUN / "candidates.tsv", cap=60)
    data["run_dir"] = str(_RUN)
    return data


def _group_payload(key: str, what: str) -> dict:
    gdir = _RUN / "groups" / key
    if what == "report":
        p = gdir / "report.json"
        return json.loads(p.read_text()) if p.is_file() else {}
    if what == "mutations":
        return {"rows": _read_tsv(gdir / "mutations.tsv", cap=8000)}
    if what == "candidates":
        return {"rows": _read_tsv(gdir / "candidates.tsv", cap=3000)}
    if what == "alignment":
        ids, seqs = _read_fasta(gdir / "align.fasta", cap=120)
        track = [r["region"] for r in _read_tsv(gdir / "align.regions.tsv")]
        return {"ids": ids, "seqs": seqs, "track": track,
                "capped": len(ids) >= 120}
    return {}


class Handler(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a):  # тихо
        pass

    def _send(self, code: int, ctype: str, body: bytes):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj):
        self._send(200, "application/json; charset=utf-8",
                   json.dumps(obj, ensure_ascii=False, default=str).encode("utf-8"))

    def do_GET(self):
        u = urlparse(self.path)
        q = parse_qs(u.query)
        try:
            if u.path == "/":
                self._send(200, "text/html; charset=utf-8", HTML.encode("utf-8"))
            elif u.path == "/api/manifest":
                self._json(_manifest())
            elif u.path == "/api/group":
                self._json(_group_payload(q.get("key", [""])[0], q.get("what", ["report"])[0]))
            elif u.path == "/api/image":
                key, name = q.get("key", [""])[0], q.get("name", ["tree.png"])[0]
                p = _RUN / "groups" / key / name
                if p.is_file() and name.endswith(".png"):
                    self._send(200, "image/png", p.read_bytes())
                else:
                    self._send(404, "text/plain", b"not found")
            else:
                self._send(404, "text/plain", b"not found")
        except Exception as e:  # noqa
            self._send(500, "text/plain", str(e).encode("utf-8"))


def launch(run_dir: str | Path, port: int = 8766) -> None:
    global _RUN
    _RUN = Path(run_dir)
    if not (_RUN / "manifest.json").is_file():
        log.warning("В %s нет manifest.json — сначала выполните `biocode run`.", _RUN)
    with socketserver.ThreadingTCPServer(("127.0.0.1", port), Handler) as httpd:
        log.info("BIOCODE GUI: http://127.0.0.1:%d  (прогон: %s)", port, _RUN)
        print(f"BIOCODE GUI → http://127.0.0.1:{port}   (Ctrl+C для выхода)")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nОстановлено.")


HTML = r"""<!doctype html><html lang="ru"><head><meta charset="utf-8">
<title>BIOCODE — браузер результатов</title>
<style>
 :root{--bg:#0f1420;--panel:#161d2e;--panel2:#1c2438;--line:#26304a;--tx:#e8edf3;--mut:#9fb0c8;
   --brand:#00add9;--cdr:#d94f3d;--cdr1:#f4c542;--cdr2:#f08a3c;--ok:#42a463}
 *{box-sizing:border-box} body{margin:0;font-family:Manrope,system-ui,sans-serif;background:var(--bg);color:var(--tx);height:100vh;overflow:hidden}
 .app{display:grid;grid-template-columns:300px 1fr;grid-template-rows:48px 1fr;height:100vh}
 .top{grid-column:1/3;display:flex;align-items:center;gap:14px;padding:0 16px;border-bottom:1px solid var(--line);background:var(--panel)}
 .top b{color:var(--brand)} .chips{display:flex;gap:8px;margin-left:auto;flex-wrap:wrap}
 .chip{background:var(--panel2);border:1px solid var(--line);border-radius:20px;padding:3px 10px;font-size:12px}
 .chip s{color:var(--mut);text-decoration:none} .chip b{color:var(--tx)}
 .nav{border-right:1px solid var(--line);background:var(--panel);overflow:auto;padding:10px}
 .nav h3{font-size:11px;text-transform:uppercase;color:var(--mut);margin:8px 4px}
 .nav input,.nav select{width:100%;background:var(--panel2);border:1px solid var(--line);color:var(--tx);border-radius:8px;padding:6px 8px;font-size:13px;margin-bottom:6px}
 .root{padding:8px;border-radius:8px;cursor:pointer;font-weight:700;margin-bottom:6px}
 .root:hover,.root.sel{background:var(--panel2)}
 .gitem{padding:7px 8px;border-radius:8px;cursor:pointer;font-size:13px;border:1px solid transparent}
 .gitem:hover{background:var(--panel2)} .gitem.sel{background:var(--panel2);border-color:var(--brand)}
 .gitem .k{font-family:ui-monospace,monospace;font-size:12px}
 .badges{display:flex;gap:5px;margin-top:3px} .b{font-size:10.5px;color:#0f1420;border-radius:8px;padding:0 6px;font-weight:700}
 .b.sz{background:#7c93b3} .b.mu{background:var(--cdr)} .b.cl{background:var(--ok)}
 .main{overflow:auto;padding:16px}
 .tabs{display:flex;gap:4px;border-bottom:1px solid var(--line);margin-bottom:14px;flex-wrap:wrap}
 .tab{padding:8px 14px;cursor:pointer;color:var(--mut);border-bottom:2px solid transparent;font-size:13px}
 .tab.sel{color:var(--tx);border-bottom-color:var(--brand)}
 .cards{display:flex;gap:12px;flex-wrap:wrap;margin-bottom:16px}
 .card{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:12px 16px;min-width:120px}
 .card .k{color:var(--mut);font-size:12px} .card .v{font-size:24px;font-weight:800}
 table{border-collapse:collapse;width:100%;font-size:12.5px} th,td{text-align:left;padding:5px 9px;border-bottom:1px solid var(--line)}
 th{color:var(--mut);position:sticky;top:0;background:var(--bg);cursor:pointer} td.n{text-align:right;font-variant-numeric:tabular-nums}
 .reg{padding:1px 7px;border-radius:9px;color:#10151f;font-weight:700;font-size:11px}
 .reg.CDR1{background:var(--cdr1)}.reg.CDR2{background:var(--cdr2)}.reg.CDR3{background:var(--cdr)}
 .reg.FR1,.reg.FR2,.reg.FR3,.reg.FR4{background:#cfd9e6}
 .muted{color:var(--mut)} .aln-wrap{overflow:auto;border:1px solid var(--line);border-radius:8px;max-height:70vh;background:#0b0f18}
 img.tree{max-width:100%;border-radius:8px;background:#fff}
 .clade{background:var(--panel);border:1px solid var(--line);border-left:3px solid var(--ok);border-radius:8px;padding:8px 12px;margin-bottom:8px}
 .pill{font-size:11px;background:var(--panel2);border-radius:12px;padding:1px 8px;margin-right:5px}
 .empty{color:var(--mut);padding:20px}
</style></head><body>
<div class="app">
 <div class="top"><b>BIOCODE</b><span class="muted">браузер результатов · фило­генетика антител</span>
   <div class="chips" id="chips"></div></div>
 <div class="nav">
   <div class="root sel" id="rootNode" onclick="showDash()">▣ Обзор прогона</div>
   <h3>Группы (V+J)</h3>
   <input id="q" placeholder="поиск: IGHV3-7…" oninput="renderNav()">
   <select id="sort" onchange="renderNav()">
     <option value="size">сортировать: по размеру</option>
     <option value="mut">по числу мутаций</option>
     <option value="clades">по уверенным кладам</option>
     <option value="key">по имени</option>
   </select>
   <div id="glist"></div>
 </div>
 <div class="main" id="main"><div class="empty">Загрузка…</div></div>
</div>
<script>
let MAN=null, CUR=null, TAB='overview';
const NUC={A:'#3aa564',C:'#3d7bd9',G:'#e0a327',T:'#d1553e','-':'#2a3550','.':'#2a3550',N:'#556'};
const REGc={FR1:'#cfd9e6',FR2:'#cfd9e6',FR3:'#cfd9e6',FR4:'#cfd9e6',CDR1:'#f4c542',CDR2:'#f08a3c',CDR3:'#d94f3d',other:'#223'};
const j=(u)=>fetch(u).then(r=>r.json());
const esc=(s)=>String(s==null?'':s).replace(/[&<>]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;'}[c]));

async function boot(){ MAN=await j('/api/manifest'); renderChips(); renderNav();
  const h=new URLSearchParams(location.hash.slice(1));
  if(h.get('group')){ await selGroup(h.get('group')); if(h.get('tab')) setTab(h.get('tab')); }
  else showDash();
}
function renderChips(){ const t=MAN.totals||{}; document.getElementById('chips').innerHTML=
  [['группы ok',t.groups_ok],['мутации',t.mutations],['кандидаты',t.candidates],['уверен. клады',t.confident_clades]]
  .map(([k,v])=>`<span class="chip"><s>${k}</s> <b>${v??0}</b></span>`).join(''); }

function groups(){ return (MAN.groups||[]).filter(g=>g.status==='ok'); }
function renderNav(){
  const q=(document.getElementById('q').value||'').toLowerCase();
  const s=document.getElementById('sort').value;
  let gs=groups().filter(g=>g.key.toLowerCase().includes(q));
  const key={size:g=>-(g.size||0),mut:g=>-(g.n_mutations||0),clades:g=>-(g.n_confident_clades||0),key:g=>g.key};
  gs.sort((a,b)=>{const ka=key[s](a),kb=key[s](b);return ka<kb?-1:ka>kb?1:0;});
  document.getElementById('glist').innerHTML=gs.map(g=>`
    <div class="gitem ${CUR===g.key?'sel':''}" onclick="selGroup('${g.key}')">
      <div class="k">${esc(g.key)}</div>
      <div class="badges"><span class="b sz">${g.size}</span>
        <span class="b mu">${g.n_mutations||0} мут</span>
        <span class="b cl">${g.n_confident_clades||0} клад</span></div>
    </div>`).join('') || '<div class="empty">нет групп</div>';
}

function showDash(){ CUR=null; TAB='overview'; document.getElementById('rootNode').classList.add('sel');
  renderNav();
  const t=MAN.totals||{}, cfg=MAN.config||{};
  const top=(MAN.top_candidates||[]);
  document.getElementById('main').innerHTML=`
   <div class="cards">
    ${card('групп проанализировано',t.groups_ok)}${card('пропущено/ошибок',(t.groups_skipped||0)+(t.groups_failed||0))}
    ${card('мутаций всего',t.mutations)}${card('кандидатов',t.candidates)}${card('уверенных клад',t.confident_clades)}
   </div>
   <p class="muted">Вход: ${esc(cfg.input||'')} · working_seq=${esc(cfg.working_seq||'')} · группировка=${esc(cfg.group_by||'')} · ML=IQ-TREE${cfg.run_bayes?' + Bayes=MrBayes':''}</p>
   <h3>Топ кандидатов-мутаций по всем группам</h3>
   ${candTable(top)}`;
}

async function selGroup(key){ CUR=key; TAB='overview';
  document.getElementById('rootNode').classList.remove('sel'); renderNav();
  const rep=await j('/api/group?key='+encodeURIComponent(key)+'&what=report'); CUR_REP=rep;
  renderGroup();
}
let CUR_REP=null;
const TABS=[['overview','Обзор'],['tree','Дерево'],['align','Выравнивание'],['mut','Мутации'],['cand','Кандидаты'],['clades','Клады']];
function renderGroup(){
  document.getElementById('main').innerHTML=
    `<h2 style="margin:0 0 4px;font-family:ui-monospace,monospace">${esc(CUR)}</h2>
     <div class="tabs">${TABS.map(([id,l])=>`<div class="tab ${TAB===id?'sel':''}" onclick="setTab('${id}')">${l}</div>`).join('')}</div>
     <div id="tabc"><div class="empty">…</div></div>`;
  drawTab();
}
function setTab(t){ TAB=t; renderGroup(); }

async function drawTab(){
  const c=document.getElementById('tabc'); const key=CUR, rep=CUR_REP||{};
  if(TAB==='overview'){
    const tr=rep.tree||{};
    c.innerHTML=`<div class="cards">
       ${card('последовательностей',rep.n_records)}${card('мутаций',rep.n_mutations)}
       ${card('кандидатов',rep.n_candidates)}${card('уверенных клад',(rep.confident_clades||[]).length)}</div>
       <p class="muted">Модель ML: <b>${esc(tr.model)}</b> · метод ${esc(tr.method)} · outgroup ${esc(tr.outgroup)} · время ${esc(tr.runtime_s)}с</p>
       <img class="tree" src="/api/image?key=${encodeURIComponent(key)}&name=tree.png">`;
  } else if(TAB==='tree'){
    c.innerHTML=`<img class="tree" src="/api/image?key=${encodeURIComponent(key)}&name=tree.png">
                 <h3>FR/CDR-разметка выравнивания</h3>
                 <img class="tree" src="/api/image?key=${encodeURIComponent(key)}&name=regions.png">`;
  } else if(TAB==='align'){
    c.innerHTML='<div class="empty">рендер выравнивания…</div>';
    const a=await j('/api/group?key='+encodeURIComponent(key)+'&what=alignment'); drawAlign(c,a);
  } else if(TAB==='mut'){
    const d=await j('/api/group?key='+encodeURIComponent(key)+'&what=mutations'); c.innerHTML=mutTable(d.rows);
  } else if(TAB==='cand'){
    const d=await j('/api/group?key='+encodeURIComponent(key)+'&what=candidates'); c.innerHTML=candTable(d.rows);
  } else if(TAB==='clades'){
    c.innerHTML=(rep.confident_clades||[]).map(cl=>`<div class="clade">
      <b>${esc(cl.clade)}</b> · размер ${cl.size} · UFBoot ${cl.ufboot} · aLRT ${cl.alrt}
      ${cl.posterior!=null?`· posterior ${cl.posterior} ${cl.confident_both_models?'✔ обе модели':''}`:''}
      <div class="muted" style="margin-top:4px">определяющих мутаций: ${cl.defining_mutations} (в CDR: ${cl.defining_cdr}) · изотипы: ${esc(JSON.stringify(cl.isotypes))}</div>
    </div>`).join('') || '<div class="empty">нет уверенных клад (UFBoot≥95 & aLRT≥80)</div>';
  }
}

function card(k,v){return `<div class="card"><div class="k">${k}</div><div class="v">${v??0}</div></div>`;}
function regBadge(r){return `<span class="reg ${esc(r)}">${esc(r)}</span>`;}
function mutTable(rows){ if(!rows||!rows.length) return '<div class="empty">нет мутаций</div>';
  return `<div class="muted">строк: ${rows.length}</div><div class="aln-wrap"><table>
   <tr><th>ветвь</th><th>поз</th><th>регион</th><th>замена</th><th>тип</th><th>внутр.</th><th>support</th></tr>
   ${rows.slice(0,3000).map(r=>`<tr><td>${esc(r.branch)}</td><td class="n">${esc(r.position)}</td>
     <td>${regBadge(r.region)}</td><td>${esc(r.ref)}→${esc(r.alt)}</td><td>${esc(r.kind)}</td>
     <td>${r.branch_internal==='True'?'да':''}</td><td class="muted">${esc(r.support)}</td></tr>`).join('')}
   </table></div>`; }
function candTable(rows){ if(!rows||!rows.length) return '<div class="empty">нет кандидатов</div>';
  return `<div class="aln-wrap"><table>
   <tr><th>группа</th><th>поз</th><th>регион</th><th>замена</th><th>тип</th><th>ветвей</th><th>UFBoot</th><th>скор</th></tr>
   ${rows.slice(0,1500).map(r=>`<tr><td class="muted">${esc(r.group)}</td><td class="n">${esc(r.position)}</td>
     <td>${regBadge(r.region)}</td><td>${esc(r.ref)}→${esc(r.alt)}</td><td>${esc(r.kind)}</td>
     <td class="n">${esc(r.n_branches)}</td><td class="n">${esc(r.max_support)}</td>
     <td class="n"><b>${esc(r.score)}</b></td></tr>`).join('')}
   </table></div>`; }

function drawAlign(c,a){
  if(!a.ids||!a.ids.length){ c.innerHTML='<div class="empty">нет выравнивания</div>'; return; }
  const cw=8, ch=13, lblW=150, W=a.track.length*cw, rows=a.ids.length;
  const wrap=document.createElement('div'); wrap.className='aln-wrap';
  const cv=document.createElement('canvas'); cv.width=lblW+W; cv.height=(rows+2)*ch+6;
  const x=cv.getContext('2d'); x.font='11px ui-monospace,monospace'; x.textBaseline='top';
  // FR/CDR трек
  for(let ci=0;ci<a.track.length;ci++){ x.fillStyle=REGc[a.track[ci]]||'#223'; x.fillRect(lblW+ci*cw,0,cw,ch-1); }
  x.fillStyle='#9fb0c8'; x.fillText('FR/CDR',4,1);
  // последовательности
  for(let ri=0;ri<rows;ri++){ const y=(ri+1)*ch+3;
    x.fillStyle='#9fb0c8'; x.fillText(a.ids[ri].slice(0,20),4,y);
    const s=a.seqs[ri];
    for(let ci=0;ci<s.length;ci++){ const chn=s[ci].toUpperCase();
      x.fillStyle=NUC[chn]||'#556'; x.fillRect(lblW+ci*cw,y,cw-0.5,ch-2);
      x.fillStyle='#0b0f18'; x.fillText(chn,lblW+ci*cw+1,y); }
  }
  wrap.appendChild(cv);
  c.innerHTML=`<div class="muted">последовательностей: ${a.ids.length}${a.capped?' (показаны первые 120)':''} · колонок: ${a.track.length} · IMGT-заякорено, метки FR/CDR перенесены на MSA</div>`;
  c.appendChild(wrap);
}
boot();
</script></body></html>"""
