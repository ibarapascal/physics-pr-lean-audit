# Inventory pass 1, stage 1: serial, resumable GitHub search for AI-attributed merged PRs. State in state.json.
# Requires an authenticated GitHub CLI (`gh`). Results depend on the live GitHub search API at run time.
import json, subprocess, time, os, sys, urllib.parse
D = os.path.dirname(os.path.abspath(__file__))
STATE = os.path.join(D, 'state.json')
LOG = os.path.join(D, 'calls.log')
REPOS = """modelica/ModelicaStandardLibrary lbl-srg/modelica-buildings ibpsa/modelica-ibpsa RWTH-EBC/AixLib UdK-VPT/BuildingSystems TransiEnt-Library/TransiEnt casella/ThermoPower xrg-simulation/ClaRa-official
SciML/ModelingToolkitStandardLibrary.jl SciML/ModelingToolkit.jl ClapeyronThermo/Clapeyron.jl CliMA/Thermodynamics.jl CliMA/Oceananigans.jl CliMA/ClimaAtmos.jl lanl-ansi/PowerModels.jl
pybamm-team/PyBaMM e2nIEE/pandapower PyPSA/PyPSA IDAES/idaes-pse Cantera/cantera CoolProp/CoolProp CalebBell/thermo CalebBell/fluids CalebBell/ht CalebBell/chemicals pvlib/pvlib-python wind-python/windpowerlib OpenMDAO/OpenMDAO usnistgov/fipy sfepy/sfepy devitocodes/devito FEniCS/dolfinx firedrakeproject/firedrake NREL/BLAST-Lite NREL/ssc NatLabRockies/EnergyPlus SINTEF/thermopack SanPen/GridCal
idaholab/moose OpenFOAM/OpenFOAM-dev""".split()
W = 'merged:2025-01-01..2026-10-16'
AUTH = ['app/devin-ai-integration','app/copilot-swe-agent','Copilot','app/chatgpt-codex-connector','app/google-labs-jules','openhands-agent','app/claude','app/sweep-ai']
TEXT = ['"Generated with Claude Code"','"Co-Authored-By: Claude"','"claude.ai/code"','"Devin run"','"chatgpt.com/codex"','"Codex"']
COMM = ['"Co-Authored-By: Claude"','"copilot-swe-agent"','"devin-ai-integration"']
def jobs():
    for r in REPOS:
        for a in AUTH: yield (r, 'author:'+a, 'issues', f'repo:{r} is:pr is:merged {W} author:{a}')
        for t in TEXT: yield (r, 'text:'+t, 'issues', f'repo:{r} is:pr is:merged {W} {t}')
        for c in COMM: yield (r, 'commit:'+c, 'commits', f'repo:{r} {c} committer-date:2025-01-01..2026-10-16')
def load():
    return json.load(open(STATE)) if os.path.exists(STATE) else {}
def save(s):
    json.dump(s, open(STATE+'.tmp','w')); os.replace(STATE+'.tmp', STATE)
def call(path, params, retries_log):
    for attempt in range(8):
        args = ['gh','api',path,'-X','GET']
        for k,v in params.items(): args += ['-f', f'{k}={v}']
        p = subprocess.run(args, capture_output=True, text=True)
        open(LOG,'a').write(f'{time.time():.0f} {path} {params.get("q","")[:90]} rc={p.returncode}\n')
        if p.returncode == 0: return json.loads(p.stdout)
        err = p.stderr + p.stdout
        if any(x in err.lower() for x in ['rate limit','403','429','secondary','abuse']):
            retries_log.append(f'attempt{attempt}: rate-limited, sleep 65s')
            time.sleep(65); continue
        if any(x in err for x in ['502','503','504','timeout']):
            retries_log.append(f'attempt{attempt}: server err'); time.sleep(10); continue
        return {'error': err[:300]}
    return {'error': 'retries exhausted'}
def main():
    s = load()
    for repo, qname, kind, q in jobs():
        key = f'{repo}|{qname}'
        if key in s and s[key].get('done'): continue
        ent = {'q': q, 'items': [], 'retries': [], 'done': False}
        page = 1
        while True:
            time.sleep(2.2)
            d = call(f'search/{kind}', {'q': q, 'per_page': 100, 'page': page}, ent['retries'])
            if 'error' in d:
                ent['error'] = d['error']; break
            ent['total'] = d.get('total_count')
            for it in d.get('items', []):
                if kind == 'issues':
                    ent['items'].append({'number': it['number'], 'url': it['html_url'], 'title': it['title'], 'author': it['user']['login'],
                        'merged_at': (it.get('pull_request') or {}).get('merged_at'), 'body': (it.get('body') or '')})
                else:
                    ent['items'].append({'sha': it['sha'], 'msg': it['commit']['message'], 'author': (it.get('author') or {}).get('login'), 'url': it['html_url']})
            if len(d.get('items', [])) < 100 or page >= 10: break
            page += 1
        ent['done'] = 'error' not in ent
        s[key] = ent; save(s)
        print(key, ent.get('total'), len(ent['items']), ent.get('error',''), flush=True)
main()
