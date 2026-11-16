# Journal Verification Record

## Method
Acceptance timelines were measured directly, not taken from marketing pages:
- **PubMed `<History>` census** (`research_notes/scripts/pubmed_timeline_census.py`): for each MEDLINE-indexed journal, the most recent ≤400 research articles (published 2025-11-15…2026-11-15, excluding reviews/editorials/errata) were pulled via NCBI E-utilities; days = `accepted` − `received` from the deposited `PubMedPubDate` history.
- **Crossref assertion census** (`research_notes/scripts/crossref_timeline_census.py`): for publishers that deposit `received`/`accepted` assertions (Springer Nature, IET/Wiley, T&F), all journal-articles in the window were pulled and the same difference computed.
- Full results: `timeline_census_table.md`; raw JSON in `research_notes/*census*.json`.

## JCIM (selected) — key facts
- **Submission→acceptance: median 84 days** (n=347; Q1 63 d, Q3 112 d; 64% ≤100 d, 22.5% ≤60 d). PRIMARY, from PubMed history. Satisfies the ≤100-day requirement.
- **JIF ≈ 6.4 (2025 JIF year), Q1** — secondary (aggregators). To re-confirm on the ACS journal page at submission.
- **Publisher:** ACS (legitimate society publisher; not predatory/hijacked). Low, sporadic retraction rate.
- **Data & Software Availability policy:** mandatory per ACS editorial 10.1021/acs.jcim.0c01389; our study ships full open code + uses only public data to comply.
- **Scope fit:** verified by a 235-paper Crossref scan of 2022–2026 JCIM methodology papers (`research_notes/jcim_methodology_corpus.json`).

## Distinctions explicitly respected (per the user's rules)
- "First decision" vs "acceptance": OUP/ACS report median *first-decision* times of ~8–20 days — these are **not** acceptance times and are **not** used here. We report submission→**acceptance** only.
- Where a journal does not deposit received/accepted dates (Elsevier, IEEE), its census row is marked "n/a" rather than guessed; such journals were assessed from publisher-reported "submission to acceptance" where available (ESWA 147–187 d) and otherwise marked UNVERIFIED.

## Timeline classification of the shortlist
- ≤100-day VERIFIED (primary census): **JCIM (84 d)**, J Chem Theory Comput (78 d), Comput Struct Biotechnol J (92 d), several ACS chemistry titles.
- 100–160 d: Journal of Cheminformatics (132), Briefings in Bioinformatics (124), Bioinformatics (147), Scientific Reports (136).
- OVER 100 d / UNVERIFIED-but-reported-high: Elsevier KBS/ESWA/Information Sciences, Springer DMKD/Machine Learning/NCA/Applied Intelligence/AI Review.
