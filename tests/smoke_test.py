"""End-to-end smoke test: launch server.py over stdio and call each tool."""
import asyncio, base64, json, sys
from pathlib import Path
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

SERVER = str(Path(__file__).resolve().parents[1] / "server.py")

def save_png(r, name):
    imgs = [c for c in r.content if getattr(c, "type", "") == "image"]
    assert imgs, "no image returned"
    mime = getattr(imgs[0], "mime_type", None) or getattr(imgs[0], "mimeType", None)
    assert mime == "image/png", mime
    out = Path(__file__).with_name(name)
    out.write_bytes(base64.b64decode(imgs[0].data)); print("figure ->", out)

def payload(r):
    texts = [c for c in r.content if getattr(c, "type", "") == "text"]
    if len(r.content) > 1:  # mixed content (text + image)
        return json.loads(texts[0].text)
    if getattr(r, "structured_content", None) or getattr(r, "structuredContent", None):
        sc = getattr(r, "structured_content", None) or r.structuredContent
        return sc.get("result", sc)
    return json.loads(r.content[0].text)

async def main():
    params = StdioServerParameters(command=sys.executable, args=[SERVER])
    async with stdio_client(params) as (rd, wr):
        async with ClientSession(rd, wr) as s:
            await s.initialize()
            tools = await s.list_tools()
            print("tools:", [t.name for t in tools.tools])
            r = payload(await s.call_tool("predict_risk", {"patient_id": "216", "prediction_time": 4}))
            print("216 risk at t=4:", r["risk_percent"], "%")
            raw = await s.call_tool("explain_risk_change", {"patient_id": "216", "from_time": 3, "to_time": 4})
            save_png(raw, "waterfall_216.png")
            r = payload(raw)
            print(f"216: {r['risk_from_percent']}% -> {r['risk_to_percent']}%")
            for c in r["contributions"]:
                print(f"  {c['component']:12s} {c['contribution_pp']:+.1f} pp")
            assert abs(r["risk_from"] - 0.0921) < 1e-3 and abs(r["risk_to"] - 0.3867) < 1e-3
            assert abs(r["closure_error"]) < 1e-10
            raw = await s.call_tool("explain_risk_change_path", {"patient_id": "216", "from_time": 3, "to_time": 4})
            save_png(raw, "evidence_216.png")
            r = payload(raw)
            print(f"method B ({r['compute_seconds']}s):")
            for c in r["contributions"]:
                print(f"  {c['component']:12s} B {c['psi_B_pp']:+.2f}  A {c['shapley_A_pp']:+.2f} pp")
            b = {c["component"]: c["psi_B_pp"] for c in r["contributions"]}
            assert abs(b["serBilir"] - 25.1) < 0.2 and abs(b["platelets"] + 9.42) < 0.1
            assert abs(r["checks"]["sum_psi_minus_binned_dR_data"]) < 1e-6
            bad = await s.call_tool("predict_risk", {"patient_id": "99999"})
            print("unknown id -> isError:", getattr(bad, "is_error", getattr(bad, "isError", None)), bad.content[0].text[:80])
            print("OK")

asyncio.run(main())
