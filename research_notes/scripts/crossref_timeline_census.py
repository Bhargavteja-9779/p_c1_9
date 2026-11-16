"""Census of submission(received)->acceptance days from Crossref 'assertion' metadata.
Covers all journal-articles with print/online publication date in the window.
Only publishers that deposit received/accepted assertions are measurable (Springer Nature, IET/Wiley, T&F)."""
import json, sys, time, urllib.request, urllib.parse, statistics
from datetime import datetime
UA = {'User-Agent': 'journal-timeline-census/1.0 (mailto:anonymous@example.org)'}
def get(u):
    for k in range(5):
        try:
            return json.load(urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60))
        except Exception as e:
            time.sleep(2 * (k + 1))
    raise RuntimeError(u)
def parse(s):
    for fmt in ('%d %B %Y', '%Y-%m-%d', '%d %b %Y', '%B %d, %Y'):
        try: return datetime.strptime(s.strip(), fmt)
        except Exception: pass
    return None
def census(issn, start='2025-11-15', end='2026-11-15', max_items=1500):
    cursor, items = '*', []
    while len(items) < max_items:
        q = f'https://api.crossref.org/journals/{issn}/works?filter=from-pub-date:{start},until-pub-date:{end},type:journal-article&rows=500&cursor={urllib.parse.quote(cursor)}&select=DOI,title,assertion,subject'
        m = get(q)['message']
        batch = m['items']
        if not batch: break
        items += batch; cursor = m['next-cursor']
        if len(batch) < 500: break
    days = []
    for it in items:
        a = {x.get('name'): x.get('value') for x in it.get('assertion', []) or []}
        r, ac = parse(a.get('received', '') or ''), parse(a.get('accepted', '') or '')
        if r and ac and ac >= r:
            days.append((ac - r).days)
    return len(items), days
if __name__ == '__main__':
    journals = json.loads(sys.argv[1])
    out = {}
    for name, issn in journals.items():
        n, d = census(issn)
        if d:
            d.sort(); q = statistics.quantiles(d, n=4)
            out[name] = dict(issn=issn, articles=n, with_dates=len(d), median=statistics.median(d), q1=q[0], q3=q[2], mean=round(statistics.mean(d),1),
                             pct_le_60=round(100*sum(x<=60 for x in d)/len(d),1), pct_le_100=round(100*sum(x<=100 for x in d)/len(d),1))
        else:
            out[name] = dict(issn=issn, articles=n, with_dates=0)
        print(name, out[name], flush=True)
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
