"""End-to-end smoke test: launch server.py over stdio and call each tool."""
import asyncio, json, sys
from pathlib import Path
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

SERVER = str(Path(__file__).resolve().parents[1] / "server.py")

def payload(r):
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
            r = payload(await s.call_tool("explain_risk_change", {"patient_id": "216", "from_time": 3, "to_time": 4}))
            print(f"216: {r['risk_from_percent']}% -> {r['risk_to_percent']}%")
            for c in r["contributions"]:
                print(f"  {c['component']:12s} {c['contribution_pp']:+.1f} pp")
            assert abs(r["risk_from"] - 0.0921) < 1e-3 and abs(r["risk_to"] - 0.3867) < 1e-3
            assert abs(r["closure_error"]) < 1e-10
            bad = await s.call_tool("predict_risk", {"patient_id": "99999"})
            print("unknown id -> isError:", getattr(bad, "is_error", getattr(bad, "isError", None)), bad.content[0].text[:80])
            print("OK")

asyncio.run(main())
