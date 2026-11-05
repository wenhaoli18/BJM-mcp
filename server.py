"""MCP server exposing BJM (Backward Joint Model) dynamic prediction as tools.

Each tool shells out to R/backend.R, which only calls BJM's public API.
Run:  python server.py            (stdio transport, for Claude Desktop / Claude Code)
"""

from __future__ import annotations

import asyncio
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

try:  # mcp >= 2
    from mcp.server.mcpserver import Image, MCPServer as _Server
    from mcp.server.mcpserver.exceptions import ToolError
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP as _Server, Image
    from mcp.server.fastmcp.exceptions import ToolError

HERE = Path(__file__).resolve().parent
BACKEND = HERE / "R" / "backend.R"
RSCRIPT = os.environ.get("RSCRIPT", shutil.which("Rscript") or "Rscript")
TIMEOUT = float(os.environ.get("BJM_TIMEOUT", "300"))
# R must run in a UTF-8 locale to draw non-ASCII labels (e.g. Chinese figures).
R_ENV = {**os.environ, "LC_ALL": os.environ.get("BJM_LOCALE", "C.UTF-8")} if os.name != "nt" else None

INSTRUCTIONS = """\
Tools for dynamic risk prediction with a fitted Backward Joint Model (BJM).
Typical questions: "What is patient 216's 2-year risk?" -> predict_risk.
"Why did it change?" -> explain_risk_change (Shapley split of the change into
time elapsed plus each biomarker's new measurements; parts sum exactly to the change).
"How did the new data shift the risk over time?" / "compare methods" ->
explain_risk_change_path (method B: evidence curves + path decomposition).
Call model_info first if unsure about the outcome, markers or time unit.
Report risks as percentages and contributions in percentage points, name the
time window, and pass on the model's caveats. Do not present results as
clinical recommendations.
"""

mcp = _Server("bjm", instructions=INSTRUCTIONS)


async def _call(tool: str, **kw: Any) -> dict:
    req = {"tool": tool, **{k: v for k, v in kw.items() if v is not None}}
    proc = await asyncio.create_subprocess_exec(
        RSCRIPT, "--vanilla", str(BACKEND), json.dumps(req),
        stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE, env=R_ENV,
    )
    try:
        out, err = await asyncio.wait_for(proc.communicate(), TIMEOUT)
    except asyncio.TimeoutError:
        proc.kill()
        raise ToolError(f"BJM backend timed out after {TIMEOUT:.0f}s")
    try:
        res = json.loads(out.decode())
    except json.JSONDecodeError:
        raise ToolError(f"BJM backend failed: {err.decode()[-2000:]}")
    if not res.get("ok"):
        raise ToolError(res.get("error", "unknown BJM error"))
    return res["result"]


async def _call_with_plot(tool: str, plot: bool, plot_lang: str | None, **kw: Any) -> list:
    """Run a tool that can draw a figure; return [JSON text, PNG image]."""
    with tempfile.TemporaryDirectory() as tmp:
        png = Path(tmp) / "figure.png"
        res = await _call(tool, plot_file=str(png) if plot else None, plot_lang=plot_lang, **kw)
        out: list = [json.dumps(res, ensure_ascii=False)]
        if plot and png.exists():
            out.append(Image(data=png.read_bytes(), format="png"))
    return out


@mcp.tool()
async def model_info() -> dict:
    """Describe the loaded BJM: dataset, outcome, biomarkers, time unit, caveats."""
    return await _call("model_info")


@mcp.tool()
async def patient_history(patient_id: str) -> dict:
    """Baseline covariates and all biomarker measurements (by visit time) for one patient."""
    return await _call("patient_history", patient_id=patient_id)


@mcp.tool()
async def predict_risk(
    patient_id: str,
    horizon: float | None = None,
    prediction_time: float | None = None,
) -> dict:
    """Dynamic risk of the event within `horizon` (default 2, in the model's time
    unit) after `prediction_time`, given event-free at that time and the patient's
    biomarker history up to it. `prediction_time` defaults to the latest visit."""
    return await _call("predict_risk", patient_id=patient_id,
                       horizon=horizon, prediction_time=prediction_time)


@mcp.tool()
async def explain_risk_change(
    patient_id: str,
    from_time: float | None = None,
    to_time: float | None = None,
    horizon: float | None = None,
    plot: bool = True,
    plot_lang: str | None = None,
) -> list:
    """Explain why a patient's risk changed between two time points.

    Computes risk at `from_time` and `to_time` (defaults: the previous and the
    latest visit) and splits the change exactly into: 'time' (moving the
    prediction time forward with no new data) plus one Shapley contribution per
    biomarker for its new measurements in (from_time, to_time]. Returns
    contributions sorted by size, with marker values at both times. Takes ~5 s.

    With plot=True (default) also returns a PNG waterfall chart for the
    clinician: earlier risk -> time elapsed -> each marker (red raises, blue
    lowers, in pp) -> later risk. Show it to the user. plot_lang: "zh"
    (default, falls back to "en" if no Chinese font) or "en"."""
    return await _call_with_plot("explain_risk_change", plot, plot_lang, patient_id=patient_id,
                                 from_time=from_time, to_time=to_time, horizon=horizon)


@mcp.tool()
async def explain_risk_change_path(
    patient_id: str,
    from_time: float | None = None,
    to_time: float | None = None,
    horizon: float | None = None,
    plot: bool = True,
    plot_lang: str | None = None,
) -> list:
    """Method B: likelihood-ratio path decomposition of a risk change.

    Same question and defaults as explain_risk_change (method A), but explains
    the change through the new data's *evidence about event time*: per-marker
    evidence curves phi_k(t) (log-likelihood-ratio of new data over event-time
    bins) and psi_k = contribution of marker k, integrated along evidence
    strength 0 -> 1. Returns psi (method B) next to Shapley (method A) per
    component, the evidence curves (bins of time after to_time; positive values
    = new data make an event in that bin more likely), and numerical checks.
    Use when the user asks *how* a marker's new data shifted the timing of
    risk, or wants to compare the two attribution methods. Slower (~10-60 s).

    With plot=True (default) also returns a PNG of the evidence curves for the
    clinician: x = future event time, shaded = risk window, each curve labelled
    with that marker's contribution in pp. Show it to the user. plot_lang: "zh"
    (default, falls back to "en" if no Chinese font) or "en"."""
    return await _call_with_plot("explain_risk_change_path", plot, plot_lang, patient_id=patient_id,
                                 from_time=from_time, to_time=to_time, horizon=horizon)


if __name__ == "__main__":
    mcp.run()
