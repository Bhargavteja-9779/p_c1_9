"""Census of received->accepted days from PubMed <History> (PubMedPubDate received/accepted) for journals indexed in PubMed.
Uses NCBI E-utilities (no key; polite rate)."""
import json, sys, time, urllib.request, urllib.parse, statistics, xml.etree.ElementTree as ET
from datetime import date
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(u):
    for k in range(5):
        try:
            return urllib.request.urlopen(urllib.request.Request(u, headers={'User-Agent': 'timeline-census/1.0'}), timeout=60).read()
        except Exception:
            time.sleep(3 * (k + 1))
    raise RuntimeError(u)
def census(ta, start='2025/11/15', end='2026/11/15', retmax=400):
    term = f'"{ta}"[ta] AND ("{start}"[dp] : "{end}"[dp]) AND journal article[pt] NOT review[pt] NOT editorial[pt] NOT comment[pt] NOT erratum[pt]'
    ids = json.loads(get(E + 'esearch.fcgi?' + urllib.parse.urlencode(dict(db='pubmed', term=term, retmax=retmax, retmode='json'))))['esearchresult']['idlist']
    days = []
    for i in range(0, len(ids), 100):
        x = ET.fromstring(get(E + 'efetch.fcgi?' + urllib.parse.urlencode(dict(db='pubmed', id=','.join(ids[i:i+100]), retmode='xml'))))
        for art in x.findall('.//PubmedArticle'):
            h = {}
            for p in art.findall('.//History/PubMedPubDate'):
                try: h[p.get('PubStatus')] = date(int(p.findtext('Year')), int(p.findtext('Month')), int(p.findtext('Day')))
                except Exception: pass
            if 'received' in h and 'accepted' in h and h['accepted'] >= h['received']:
                days.append((h['accepted'] - h['received']).days)
        time.sleep(0.4)
    return len(ids), days
if __name__ == '__main__':
    out = {}
    for name, ta in json.loads(sys.argv[1]).items():
        try:
            n, d = census(ta)
        except Exception as e:
            print(name, 'ERR', e); continue
        if len(d) >= 5:
            q = statistics.quantiles(sorted(d), n=4)
            out[name] = dict(ta=ta, pubmed_articles=n, with_dates=len(d), median=statistics.median(d), q1=q[0], q3=q[2], mean=round(statistics.mean(d), 1),
                             pct_le_60=round(100*sum(x<=60 for x in d)/len(d),1), pct_le_100=round(100*sum(x<=100 for x in d)/len(d),1))
        else:
            out[name] = dict(ta=ta, pubmed_articles=n, with_dates=len(d))
        print(name, out[name], flush=True)
        time.sleep(0.5)
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
