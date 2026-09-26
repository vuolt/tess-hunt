"""tess-hunt: a local app for exploring and following up single-transit planet candidates.

    streamlit run app/main.py

The app only reads the pipeline's results. It writes nothing except its own
workflow-status file (app/workflow_status.json) and, when you start a search,
the run's log in work/app_runs/. It never logs in to, or submits anything to,
an external website: you do that yourself with the text and files it prepares.
"""

from __future__ import annotations

import os
import sys
from datetime import date

import numpy as np
import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app import data, runner, texts, workflow as wf  # noqa: E402

T = texts.term
st.set_page_config(page_title="tess-hunt", page_icon="🔭", layout="wide")


def md(s: str):
    st.markdown(s, unsafe_allow_html=True)


# ------------------------------------------------------------------ helpers

@st.cache_data(ttl=60)
def all_candidates(include_drop: bool) -> pd.DataFrame:
    rows = []
    for s in data.sectors():
        fu = data.followup(s)
        if fu is None:
            continue
        fu = fu[fu.role == "candidate"].copy()
        fu["sector"] = s
        rows.append(fu)
    if not rows:
        return pd.DataFrame()
    df = pd.concat(rows, ignore_index=True)
    p5 = data.phase5()
    if p5 is not None:
        m = p5[p5.role == "candidate"][["tic", "sector", "phase5", "verdict"]]
        df = df.merge(m, on=["tic", "sector"], how="left")
        df["status"] = df.phase5.fillna(df.category)
    else:
        df["status"], df["verdict"] = df.category, None
    if not include_drop:
        df = df[df.status.isin(["submit", "maybe"])]
    order = {"submit": 0, "maybe": 1, "drop": 2}
    return df.sort_values(["status", "sector", "order"], key=lambda c: c.map(order) if c.name == "status" else c)


def card_for(row) -> dict:
    return data.card(int(row["sector"]), row, data.phase5(), data.phase2_candidates(int(row["sector"])))


STATUS_COLOUR = {"submit": "🟢", "maybe": "🟡", "drop": "🔴"}


def show_card(c: dict, compact: bool = False):
    md(f"### {STATUS_COLOUR.get(c['status'], '')} {T('TIC')} {c['tic']} — "
       f"**{c['status']}** <span style='color:gray'>({T('sector', 'Sector')} {c['sector']})</span>")
    dist = f", about {c['distance_pc']:.0f} parsecs ({c['distance_pc'] * 3.26:.0f} light years) away" \
        if c.get("distance_pc") else ""
    md(f"**Star:** {c['star']}{dist}" + (f", {c['radius_rsun']:.2f}× the Sun's radius" if c.get("radius_rsun") else "") + ".")
    md(f"**{T('dip', 'Dip')}:** {c['depth_ppm']:.0f} {T('ppm')} deep "
       f"({c['depth_ppm'] / 1e4:.2f} % fainter), lasting {c['duration_h']:.1f} hours.")
    md(f"**Estimated size:** {c['size']}.")
    md(f"**Possible orbit:** {c['orbit']}")
    if c.get("verdict") and isinstance(c["verdict"], str):
        md(f"**Expert checks:** {c['verdict']}.")
    if c["reasons"]:
        md("**Why this status:** " + "; ".join(c["reasons"]) + ".")
    elif c["status"] == "submit":
        md("**Why this status:** every check we ran supports a planet.")
    if not compact:
        extra = []
        if c.get("snr") and np.isfinite(c["snr"]):
            extra.append(f"{T('SNR')} {c['snr']:.1f}")
        if c.get("fpp") is not None and np.isfinite(c["fpp"]):
            extra.append(f"{T('FPP')} {c['fpp']:.3f}")
        if extra:
            md("<span style='color:gray'>" + " · ".join(extra) + "</span>")


# ------------------------------------------------------------------ pages

def page_run():
    st.header("Run")
    md("Two kinds of run, both in the background: you can close this page and come back. When a run "
       "finishes, this page says in plain words what it found (or that nothing was found) and gives "
       "you a ready-made message to save the results. The app itself never commits or uploads anything.")
    cur = runner.current()
    running = bool(cur and cur["running"])
    procs_max = os.cpu_count() or 4
    t1, t2 = st.tabs(["🔭 Search a new sector for planets", "🕵️ Check other people's candidates"])
    with t1:
        md(f"Search one {T('TESS')} {T('sector')} for single {T('transit', 'transits')}: select the "
           "stars, search every light curve, vet the dips, follow up the best and run the expert "
           "checks. A full sector takes 5–8 hours.")
        searched = data.searched()
        if len(searched):
            md("Already searched: " + ", ".join(f"Sector {int(x)}" for x in searched.sector) + ".")
        c1, c2, c3 = st.columns([1, 1, 2])
        sector = c1.number_input("Sector", 1, 200, value=int(cur["sector"]) if cur and cur.get("sector")
                                 else 49, step=1, disabled=running)
        procs = c2.number_input("CPU cores to use", 1, procs_max, value=min(3, procs_max), disabled=running)
        skip = c3.checkbox("Skip the final sensitivity test (saves ~1–2 h)", value=False, disabled=running)
        done = data.searched_sector(int(sector))
        finished = bool(done and done["stage"] == "phase4")
        force = False
        if done:
            from tesshunt import ledger
            st.warning(ledger.describe(done) + (" Its results are on the Results and Candidates pages."
                                                if finished else " Starting it again resumes the run."))
            if finished:
                force = st.checkbox("Run this sector again anyway", value=False, disabled=running,
                                    help="Stars already searched are still not downloaded or searched "
                                         "again; the later steps are redone.")
        if st.button("▶ Start / resume the search", disabled=running or (finished and not force),
                     type="primary", key="start_sector",
                     help="Starting a sector that was interrupted continues where it stopped: every step "
                          "keeps its finished work."):
            try:
                runner.start(int(sector), int(procs), skip, force=force)
                st.success(f"Started Sector {int(sector)}.")
                st.rerun()
            except RuntimeError as e:
                st.error(str(e))
    with t2:
        fs = data.fp_status()
        md(f"Other people have submitted thousands of {T('TOI', 'TOIs')} and {T('CTOI', 'CTOIs')} that "
           "nobody has confirmed or ruled out yet. This check looks for signs that one is really two "
           "stars eclipsing each other: Gaia sees a companion star on the same period, alternate dips "
           "have different depths, or the dimming is off-centre. Every candidate checked is recorded, "
           "so nothing is checked twice, by you or anyone else.")
        md(f"So far: **{fs['checked']:,}** candidates' light curves checked, **{fs['flagged']}** flagged "
           f"as likely false positives ({fs['high']} with high confidence)."
           + (f" Lists downloaded on {fs['tables_date']}." if fs["tables_date"] else ""))
        c1, c2, c3 = st.columns([1, 1, 2])
        limit = c1.number_input("Candidates to check this run", 0, 20000, value=500, step=100,
                                disabled=running, help="0 = every candidate not checked yet. "
                                "About 15 per minute with 2 cores.")
        fprocs = c2.number_input("CPU cores to use ", 1, procs_max, value=min(2, procs_max), disabled=running)
        refresh = c3.checkbox("Download today's TOI and CTOI lists first", value=False, disabled=running,
                              help="Includes candidates submitted since the last download. The lists "
                                   "are downloaded once per run.")
        est = "all remaining candidates: several hours" if not limit else f"about {max(1, limit // 15)} min"
        st.caption(f"Estimated time: {est}.")
        if st.button("▶ Start the check", disabled=running, type="primary", key="start_fp"):
            try:
                runner.start_fp(int(fprocs), refresh, int(limit))
                st.success("Started the false-positive check.")
                st.rerun()
            except RuntimeError as e:
                st.error(str(e))
    if st.button("■ Stop the current run", disabled=not running, key="stop"):
        runner.stop()
        st.warning("Stopped. Start it again later to continue from where it stopped.")
        st.rerun()
    live_status()


@st.fragment(run_every=5)
def live_status():
    cur = runner.current()
    if not cur:
        st.info("Nothing has been started from the app yet.")
        return
    prog = runner.progress(runner.log_text(cur["log"]), cur.get("steps"))
    state = ("running" if cur["running"] else
             "finished" if prog["finished"] else
             f"failed at '{prog['failed']}'" if prog["failed"] else "stopped")
    st.subheader(f"{runner.label(cur)}: {state}")
    st.progress(prog["overall"], text=f"{prog['overall']:.0%} of all steps")
    if prog["current"] and cur["running"]:
        frac = f" ({prog['fraction']:.0%})" if prog["fraction"] is not None else ""
        st.write(f"Now: **{runner.STEP_WORDS.get(prog['current'], prog['current'])}**{frac}")
    for s in cur.get("steps") or runner.STEP_NAMES:
        mark = "✅" if s in prog["done"] else ("⏳" if s == prog["current"] else "·")
        st.write(f"{mark} {runner.STEP_WORDS.get(s, s)}")
    with st.expander("Log (last lines)"):
        st.code(runner.tail(cur["log"], 20) or "(empty)")
    if not cur["running"] and prog["finished"]:
        show_findings(cur)


def show_findings(cur):
    from tesshunt import findings
    f = findings.fp_findings() if cur.get("kind") == "fp" else findings.sector_findings(int(cur["sector"]))
    st.subheader("What this run found")
    if not f["done"]:
        st.info(f["text"])
        return
    (st.success if f["found"] else st.info)(f["title"])
    st.markdown(f["text"])
    st.subheader("Save the results")
    if not f["files"]:
        st.write("Nothing new to save: every result file is already committed.")
        return
    md(f"{len(f['files'])} result files are new or changed. Save them to the project's history "
       "(the app does not do this itself). **In a terminal**, in the project folder, paste:")
    st.code(findings.commit_commands(f["title"], f["body"]), language="bash")
    with st.expander("Using GitHub Desktop instead"):
        md("Tick the changed files under `results/` and `plots/`, paste this as the **summary**:")
        st.code(f["title"], language=None)
        md("and this as the **description**, then *Commit* and *Push origin*:")
        st.code(f["body"], language=None)
    with st.expander(f"Files to save ({len(f['files'])})"):
        st.code("\n".join(f["files"][:200]) + ("\n..." if len(f["files"]) > 200 else ""), language=None)


def page_results():
    st.header("Results by sector")
    secs = data.sectors()
    if not secs:
        st.info("No results yet. Start a search on the Run page.")
        return
    s = st.selectbox("Sector", secs, index=len(secs) - 1)
    steps = data.funnel(s)
    if not steps:
        st.info("This sector has no summary yet.")
        return
    prev = None
    for x in steps:
        n = x["n"]
        c1, c2 = st.columns([1, 3])
        c1.metric(x["step"], f"{n:,}" if n is not None else "–")
        if prev is not None and n is not None and prev >= n:
            c1.caption(f"{prev - n:,} removed")
        c2.write(x["text"])
        prev = n if n is not None else prev
    md("<span style='color:gray'>Each step removes things that are not planets. The last numbers are "
       "candidates worth a closer look, not confirmed planets.</span>")


def page_candidates():
    st.header("Candidates")
    show_drop = st.toggle("Also show dropped candidates", value=False)
    df = all_candidates(show_drop)
    if df.empty:
        st.info("No follow-up results yet (Phase 4).")
        return
    md(f"{len(df)} candidates. Each card explains the {T('dip')} in plain words; open the plots for the "
       "full detail. Hover over dotted words for an explanation.")
    for _, row in df.iterrows():
        c = card_for(row)
        with st.container(border=True):
            left, right = st.columns([3, 2])
            with left:
                show_card(c)
                if st.button("Open workflow →", key=f"wf_{row['sector']}_{row['tic']}"):
                    st.session_state["wf_key"] = f"s{int(row['sector']):04d}_{int(row['tic'])}"
                    st.session_state["page"] = "Workflow"
                    st.rerun()
            with right:
                pl = data.sheets(int(row["sector"]), int(row["tic"]))
                if "pht" in pl:
                    st.image(pl["pht"], caption="The light curve as Planet Hunters TESS shows it, "
                                                "with a zoom on the dip")
                if "followup" in pl:
                    st.image(pl["followup"], caption="Follow-up sheet: all TESS data and allowed orbits")
                with st.expander("More plots"):
                    for name, p in pl.items():
                        if name not in ("pht", "followup"):
                            st.image(p, caption=name)


def page_workflow():
    st.header("Workflow")
    df = all_candidates(True)
    if df.empty:
        st.info("No candidates yet.")
        return
    keys = [f"s{int(r.sector):04d}_{int(r.tic)}" for r in df.itertuples()]
    labels = {k: f"TIC {k.split('_')[1]} (S{int(k[1:5])}, {s})" for k, s in zip(keys, df.status)}
    default = st.session_state.get("wf_key", keys[0])
    key = st.selectbox("Candidate", keys, index=keys.index(default) if default in keys else 0,
                       format_func=lambda k: labels[k])
    st.session_state["wf_key"] = key
    row = df.iloc[keys.index(key)]
    c = card_for(row)
    try:
        statuses = wf.load()
    except (ValueError, OSError) as e:
        st.error(f"The status file {wf.status_path()} could not be read ({e}). Fix or move it; the app "
                 "will not overwrite it.")
        return
    s = wf.get(statuses, key)

    def save():
        wf.save(statuses)

    with st.container(border=True):
        show_card(c, compact=True)
    st.progress(s.progress(), text=f"{sum(s.done(x) for x in wf.STEPS)} of {len(wf.STEPS)} steps done")

    # a. feedback
    with st.expander(("✅ " if s.done("feedback") else "") + wf.STEP_TITLES["feedback"],
                     expanded=s.current_step() == "feedback"):
        md(f"Ask the {T('Planet Hunters TESS')} community for a second opinion before anything else.")
        post = texts.forum_post(c, row)
        st.text_area("Forum post (copy and paste)", post, height=320)
        plots = data.sheets(c["sector"], c["tic"])
        if "pht" in plots:
            st.image(plots["pht"], caption="Attach this one first: it looks like the light curves on the forum.")
            with open(plots["pht"], "rb") as fh:
                st.download_button("⬇ Light-curve plot (PNG)", fh.read(), file_name=os.path.basename(plots["pht"]),
                                   mime="image/png", key="pht_png")
        st.download_button("⬇ Plots to attach (zip)", texts.plot_bundle(plots),
                           file_name=f"tic{c['tic']}_plots.zip", mime="application/zip",
                           on_click=lambda: (s.mark_post_prepared() if not s.post_prepared else None, save()))
        st.link_button("Open the Planet Hunters TESS forum", texts.PHT_FORUM)
        if not s.post_prepared and st.button("I have prepared the post", key="prep"):
            s.mark_post_prepared()
            save()
            st.rerun()
        if s.post_prepared:
            with st.form("feedback"):
                url = st.text_input("Link to your forum post", s.post_url)
                fb = st.radio("What did the community say?", wf.FEEDBACK_VALUES,
                              index=wf.FEEDBACK_VALUES.index(s.feedback) if s.feedback else 2,
                              horizontal=True)
                notes = st.text_area("Notes", s.feedback_notes)
                if st.form_submit_button("Save feedback"):
                    try:
                        s.record_feedback(url, fb, notes)
                        save()
                        st.success("Saved.")
                        st.rerun()
                    except wf.WorkflowError as e:
                        st.error(str(e))

    # b. submission
    with st.expander(("✅ " if s.done("submission") else "") + wf.STEP_TITLES["submission"],
                     expanded=s.current_step() == "submission"):
        if s.feedback_done and s.feedback != "positive" and not s.override_reason:
            st.warning(f"The feedback was **{s.feedback}**, so this step is locked.")
            reason = st.text_input("Override anyway: type why you still want to submit")
            if st.button("Override"):
                try:
                    s.override(reason)
                    save()
                    st.rerun()
                except wf.WorkflowError as e:
                    st.error(str(e))
        if not s.unlocked("submission"):
            st.info(s.lock_reason("submission"))
        else:
            if s.override_reason:
                st.caption(f"Unlocked by override: {s.override_reason}")
            st.warning("**ExoFOP only accepts community candidates that are published in a refereed "
                       "journal**, and uploading needs an approved request. Check the guidelines "
                       "before you submit.")
            f = texts.ctoi_fields(c, row)
            md(f"**{T('CTOI')} fields** (copy into {T('ExoFOP')}'s form):")
            st.dataframe(pd.DataFrame([
                ("TIC", f["TIC"]), ("Transit epoch (BJD_TDB)", f"{f['epoch_bjd']} ± {f['epoch_unc']}"),
                ("Period (days)", f["period"] or "unknown (single transit)"),
                ("Period constraints", f["period_note"]),
                ("Depth (ppm)", f"{f['depth_ppm']} ± {f['depth_unc']}"),
                ("Duration (hours)", f"{f['duration_h']} ± {f['duration_unc']}"),
                ("Planet radius (Earth radii)", f["radius_re"]), ("Impact parameter", f["impact"]),
                ("Vetting summary", f["vetting"])], columns=["field", "value"]).astype(str),
                hide_index=True, use_container_width=True)
            tag = st.text_input("Your ExoFOP data tag (YYYYMMDD_username_description)", "")
            paper = st.text_input("URL of the published paper (required by ExoFOP)", "")
            row_txt = texts.exofop_row(f, tag or "YYYYMMDD_username_description", paper or "PAPER_URL")
            st.download_button("⬇ ExoFOP bulk-upload file", texts.exofop_file([row_txt]),
                               file_name=f"params_planet_{date.today():%Y%m%d}_001.txt", mime="text/plain")
            st.caption(f"Format: ExoFOP's planet-parameter template ({texts.EXOFOP_TEMPLATE}). Check "
                       f"the next free TIC{c['tic']}.nn number on ExoFOP before uploading.")
            l1, l2, l3 = st.columns(3)
            l1.link_button("Submission guidelines", texts.EXOFOP_GUIDELINES)
            l2.link_button("Request upload access", texts.EXOFOP_REQUEST)
            l3.link_button("ExoFOP bulk upload page", texts.EXOFOP_BULK_PARAMS)
            if not s.submission_done and st.button("Submission prepared — continue"):
                s.mark_submission_prepared()
                save()
                st.rerun()

    # c. record
    with st.expander(("✅ " if s.done("recorded") else "") + wf.STEP_TITLES["recorded"],
                     expanded=s.current_step() == "recorded"):
        if not s.unlocked("recorded"):
            st.info(s.lock_reason("recorded"))
        else:
            with st.form("record"):
                d = st.date_input("Date you submitted it to ExoFOP",
                                  date.fromisoformat(s.submitted_on) if s.submitted_on else date.today())
                ctoi = st.text_input("CTOI number, once ExoFOP assigns one (you can add it later)", s.ctoi)
                if st.form_submit_button("Save"):
                    try:
                        if s.recorded_done:
                            s.set_ctoi(ctoi)
                        else:
                            s.record_submission(d, ctoi)
                        save()
                        st.success("Saved.")
                        st.rerun()
                    except wf.WorkflowError as e:
                        st.error(str(e))

    # d. follow-up
    with st.expander(("✅ " if s.done("followup") else "") + wf.STEP_TITLES["followup"],
                     expanded=s.current_step() == "followup"):
        if not s.unlocked("followup"):
            st.info(s.lock_reason("followup"))
        else:
            pred, pred_md = data.predictions(c["tic"])
            if pred is not None:
                md("**Predicted future transits** (times in UTC; ±1σ in hours):")
                cols = [x for x in ("P_d", "tc_utc", "sigma_h", "tess_sector", "north_site",
                                    "north_transit_frac", "south_site", "south_transit_frac") if x in pred]
                st.dataframe(pred[cols], hide_index=True, use_container_width=True)
                with st.expander("Full prediction report"):
                    md(pred_md or "")
            else:
                md(f"No transit predictions yet: with a single {T('transit')} the {T('period')} is unknown, so "
                   "future transits cannot be predicted until a second one is found.")
            p2 = data.phase2_candidates(c["sector"])
            m = p2[p2.tic == c["tic"]] if p2 is not None else None
            if m is not None and len(m) and {"ra", "dec"} <= set(m.columns):
                fut = texts.future_tess(c["tic"], m.ra.iloc[0], m.dec.iloc[0],
                                        float(c["epoch_btjd"] or 0) + 30)
                if fut:
                    md(f"**{T('TESS')} will observe this star again** in: " +
                       ", ".join(f"Sector {x['sector']} (around {x['approx_mid']})" for x in fut) + ".")
                else:
                    md(f"According to {T('TESS')}'s published pointing plan (to Sector 134), it will **not** "
                       "observe this star again.")
            if not s.followup_done and st.button("Mark follow-up reviewed"):
                s.mark_followup_done()
                save()
                st.rerun()
    with st.expander("History"):
        for h in reversed(s.history):
            st.write(f"{h['at']}: {h['what']}")
        if st.button("Reset this candidate's workflow"):
            s.reset()
            save()
            st.rerun()


def page_glossary():
    st.header("Glossary")
    for k in sorted(texts.GLOSSARY, key=str.lower):
        md(f"**{k}** — {texts.GLOSSARY[k]}")


PAGES = {"Candidates": page_candidates, "Workflow": page_workflow, "Results": page_results,
         "Run": page_run, "Glossary": page_glossary}


def main():
    st.sidebar.title("🔭 tess-hunt")
    st.sidebar.caption("Hunting long-period planets in TESS data")
    page = st.sidebar.radio("Page", list(PAGES), index=list(PAGES).index(st.session_state.get("page", "Candidates")))
    st.session_state["page"] = page
    st.sidebar.markdown("---")
    st.sidebar.caption("Read-only: the app never changes results or pipeline code, and never submits "
                       "anything anywhere. It saves only your workflow progress "
                       "(app/workflow_status.json).")
    PAGES[page]()


main()
