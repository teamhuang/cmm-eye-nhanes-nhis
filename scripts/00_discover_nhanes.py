"""Inspect NHANES .xpt files: rows, cols, and grep the exact variable codes we will use.
This is a DISCOVERY pass — no modeling. It exists so downstream scripts cite verified
column names, not memory."""
import pyreadstat, os, glob, sys, json

DATA = r"H:\眼科论文\眼科2\data\raw\nhanes"
files = sorted(glob.glob(os.path.join(DATA, "*.xpt")))

targets = {
    "VIQ": ["viq070","viq071","viq080","viq090","viq060","viq310","viq300","viq050","viq020","viq031","viq041","viq051"],
    "OPXRET": ["opddret","opddarma","opddarm","opddlarm","opddge","opddexu","opdga","opddaga","opddma","opddhe","opddhe","opdume","opdumec","opdsret","opddnvo"],
    "CRP": ["crpp"],
    "DEMO": ["wtmec2yr","sdmvpsu","sdmvstra","ridageyr","riagendr","ridreth1","dmdeduc","dmdhhin3","dmdmzx"],
    "DIQ": ["diq010","diq020","diq050","diq060","diq295"],
    "BPX": ["bpxsy1","bpxsy2","bpxsy3","bpxsy4","bpxdi1","bpxdi2","bpxdi3","bpxdi4"],
    "TCHOL": ["tchol"],
    "TRIGLY": ["trigly"],
    "HDL": ["hdl"],
    "LDL": ["ldl"],
    "BMX": ["bmi","bmxwt","bmxht"],
    "SMQ": ["smq040","smq020","smok30d"],
    "GLU": ["glu"],
}

summary = []
for path in files:
    name = os.path.basename(path)
    prefix = name.split("_")[0].replace(".xpt","")
    try:
        df, meta = pyreadstat.read_xport(path, metadataonly=False)
    except Exception as e:
        print(f"!! {name}: LOAD FAIL {e}")
        continue
    cols = [c.upper() for c in df.columns]
    n = len(df)
    m = meta.number_columns
    print(f"\n=== {name}: rows={n} cols={m} ===")
    # print first 30 columns for context
    print("first cols:", cols[:30])
    tgt = targets.get(prefix, [])
    for t in tgt:
        hit = [c for c in cols if t.upper() in c]
        if hit:
            print(f"  target {t} -> {hit[:6]}")
    # CRP special
    if prefix.upper() == "CRP":
        crp = [c for c in cols if c.startswith("CRP")]
        print("  CRP cols:", crp[:5])
    summary.append({"file": name, "rows": n, "cols": m})

with open(os.path.join(DATA, "_inventory.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\nDONE. inventory ->", os.path.join(DATA, "_inventory.json"))
