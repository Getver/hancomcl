from pathlib import Path
from typing import Optional
import sqlite3
import io, csv
from fastapi.responses import StreamingResponse

from fastapi import FastAPI, Query, HTTPException
from fastapi.responses import RedirectResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware

# ─────────────────────────────────────────────
# 기본 경로/DB
# ─────────────────────────────────────────────
BASE = Path(__file__).resolve().parents[1]  # /disk0/yang93/eQTL
DB_PATH = str(BASE / "data" / "eqtl.db")

TABLE = "eQTL_result"

HEADERS = ["CHR","BP","SNP","REF","ALT","GENE","DISTANCE","SOURCE","SLOPE","p-value","FDR"]

SORT_MAP = {
    "CHR":"s.chr","BP":"s.bp","SNP":"s.snp_name","REF":"s.ref","ALT":"s.alt",
    "GENE":"g.gene_name","DISTANCE":"e.distance","p-value":"e.npval",
    "SLOPE":"e.slope","FDR":"e.fdr","SOURCE":"e.source_file"
}

# ─────────────────────────────────────────────
# FastAPI 기본 셋업
# ─────────────────────────────────────────────
app = FastAPI(title="eQTL", version="0.1")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], allow_methods=["*"], allow_headers=["*"], allow_credentials=True
)

app.mount("/css",    StaticFiles(directory=str(BASE / "css")),    name="css")
app.mount("/image",  StaticFiles(directory=str(BASE / "image")),  name="image")
app.mount("/script", StaticFiles(directory=str(BASE / "script")), name="script")
app.mount("/html",   StaticFiles(directory=str(BASE / "html"), html=True), name="html")
app.mount("/text",   StaticFiles(directory=str(BASE / "text")),   name="text")

@app.get("/")
def root():
    return RedirectResponse(url="/html/main.html")

# ─────────────────────────────────────────────
# DB helpers
# ─────────────────────────────────────────────
def get_db():
    con = sqlite3.connect(DB_PATH)
    con.row_factory = sqlite3.Row
    return con

def ensure_table_exists(con: sqlite3.Connection):
    row = con.execute(
        "SELECT name FROM sqlite_master WHERE (type='table' OR type='view') AND name=?",
        (TABLE,)
    ).fetchone()
    if not row:
        raise HTTPException(
            status_code=500,
            detail=f"Required table/view '{TABLE}' not found."
        )

# ─────────────────────────────────────────────
# API: 컬럼 설명(프론트 tooltip)
# ─────────────────────────────────────────────
@app.get("/api/column-info")
def column_info():
    return {
        "CHR": "The chromosome where the SNP or gene is located.",
        "BP": "The base pair position on the chromosome (hg38 reference genome).",
        "REF": "The allele present in the reference genome at this position.",
        "ALT": "The variant allele observed at this position.",
        "SNP": "Reference SNP ID (rsID) for the variant tested in the eQTL analysis.",
        "GENE": "The gene whose expression is tested for association with the SNP.",
        "DISTANCE": "Distance between the SNP and the gene’s TSS (0 if within gene).",
        "p-value": "Nominal p-value for SNP–gene association.",
        "SLOPE": "Effect size (regression coefficient).",
        "FDR": "Adjusted p-value (False Discovery Rate).",
        "SOURCE": "Dataset or file where this association was found."
    }

# ─────────────────────────────────────────────
# API: 헬스체크
# ─────────────────────────────────────────────
@app.get("/api/health")
def health():
    con = get_db()
    try:
        ensure_table_exists(con)
        cnt = con.execute(f"SELECT COUNT(1) AS c FROM {TABLE}").fetchone()["c"]
        return {"ok": True, "table": TABLE, "rows": cnt}
    finally:
        con.close()

# ─────────────────────────────────────────────
# API: eQTL 목록 (서버 페이징/정렬/필터)
# ─────────────────────────────────────────────
@app.get("/api/eqtl")
def list_eqtl(
    page: int = Query(1, ge=1),
    size: int = Query(20, ge=1, le=200),
    sort: str = "CHR",
    order: str = "asc",
    gene: Optional[str] = None,
    snp: Optional[str] = None,
    distance: Optional[int] = None,
    pValueMin: Optional[float] = None,
    pValueMax: Optional[float] = None,
    fdrMin: Optional[float] = None,
    fdrMax: Optional[float] = None,
    search: Optional[str] = None,
):
    con = get_db()
    try:
        where, params = [], {}
        if gene:
            where.append("g.gene_name LIKE :gene")
            params["gene"] = f"%{gene}%"
        if snp:
            where.append("s.snp_name LIKE :snp")
            params["snp"] = f"%{snp}%"
        if distance is not None:
            where.append("e.distance = :distance")
            params["distance"] = distance
        if pValueMin is not None:
            where.append("e.npval >= :pmin")
            params["pmin"] = pValueMin
        if pValueMax is not None:
            where.append("e.npval <= :pmax")
            params["pmax"] = pValueMax
        if fdrMin is not None:
            where.append("e.fdr >= :fmin")
            params["fmin"] = fdrMin
        if fdrMax is not None:
            where.append("e.fdr <= :fmax")
            params["fmax"] = fdrMax
        if search:
            params["q"] = f"%{search}%"
            where.append("(s.snp_name LIKE :q OR g.gene_name LIKE :q)")

        where_sql = (" WHERE " + " AND ".join(where)) if where else ""
        sort_col = SORT_MAP.get(sort, "s.chr")
        order_sql = "DESC" if order.lower() == "desc" else "ASC"

        base = f"""
          FROM eQTL_result e
          JOIN SNP s ON e.snp_id = s.snp_id
          JOIN Gene g ON e.gene_id = g.gene_id
          {where_sql}
        """

        total = con.execute(f"SELECT COUNT(*) {base}", params).fetchone()[0]
        offset = (page - 1) * size

        sql = f"""
        SELECT
            s.chr       AS CHR,
            s.bp        AS BP,
            s.snp_name  AS SNP,
            s.ref       AS REF,
            s.alt       AS ALT,
            g.gene_name AS GENE,
            e.distance  AS DISTANCE,
            e.source_file AS SOURCE,
            e.npval     AS "p-value",
            e.slope     AS SLOPE,
            e.fdr       AS FDR
        {base}
        ORDER BY {sort_col} {order_sql}
        LIMIT :size OFFSET :offset
        """
        params.update({"size": size, "offset": offset})
        rows = con.execute(sql, params).fetchall()

        items = []
        for r in rows:
            row_dict = dict(r)
            if row_dict.get("FDR") is not None:
                try: row_dict["FDR"] = round(float(row_dict["FDR"]), 4)
                except: pass
            if row_dict.get("p-value") is not None:
                try: row_dict["p-value"] = round(float(row_dict["p-value"]), 4)
                except: pass
            items.append(row_dict)

        return {"headers": HEADERS, "items": items, "total": total, "page": page, "size": size}
    finally:
        con.close()

# ─────────────────────────────────────────────
# API: Statistics card
# ─────────────────────────────────────────────
@app.get("/api/stats")
def stats():
    con = get_db()
    try:
        cnt_genes  = con.execute("SELECT COUNT(*) FROM Gene").fetchone()[0]
        cnt_snps   = con.execute("SELECT COUNT(*) FROM SNP").fetchone()[0]
        cnt_eqtls  = con.execute("SELECT COUNT(*) FROM eQTL_result").fetchone()[0]

        return {"genes": cnt_genes, "snps": cnt_snps, "eqtls": cnt_eqtls}
    finally:
        con.close()

# ─────────────────────────────────────────────
# API: Export CSV
# ─────────────────────────────────────────────
@app.get("/api/eqtl/export")
def export_eqtl(
    sort: str = "CHR",
    order: str = "asc",
    gene: Optional[str] = None,
    snp: Optional[str] = None,
    distance: Optional[int] = None,
    pValueMin: Optional[float] = None,
    pValueMax: Optional[float] = None,
    fdrMin: Optional[float] = None,
    fdrMax: Optional[float] = None,
    search: Optional[str] = None,
):
    con = get_db()
    try:
        where, params = [], {}
        if gene:
            where.append("g.gene_name LIKE :gene")
            params["gene"] = f"%{gene}%"
        if snp:
            where.append("s.snp_name LIKE :snp")
            params["snp"] = f"%{snp}%"
        if distance is not None:
            where.append("e.distance = :distance")
            params["distance"] = distance
        if pValueMin is not None:
            where.append("e.npval >= :pmin")
            params["pmin"] = pValueMin
        if pValueMax is not None:
            where.append("e.npval <= :pmax")
            params["pmax"] = pValueMax
        if fdrMin is not None:
            where.append("e.fdr >= :fmin")
            params["fmin"] = fdrMin
        if fdrMax is not None:
            where.append("e.fdr <= :fmax")
            params["fmax"] = fdrMax
        if search:
            params["q"] = f"%{search}%"
            where.append("(s.snp_name LIKE :q OR g.gene_name LIKE :q)")

        where_sql = (" WHERE " + " AND ".join(where)) if where else ""
        sort_col = SORT_MAP.get(sort, "s.chr")
        order_sql = "DESC" if order.lower() == "desc" else "ASC"

        sql = f"""
        SELECT 
            s.chr       AS CHR,
            s.bp        AS BP,
            s.snp_name  AS SNP,
            s.ref       AS REF,
            s.alt       AS ALT,
            g.gene_name AS GENE,
            e.distance  AS DISTANCE,
            e.source_file AS SOURCE,
            e.npval     AS "p-value",
            e.slope     AS SLOPE,
            e.fdr       AS FDR
        FROM eQTL_result e
        JOIN SNP s ON e.snp_id = s.snp_id
        JOIN Gene g ON e.gene_id = g.gene_id
        {where_sql}
        ORDER BY {sort_col} {order_sql}
        """
        rows = con.execute(sql, params).fetchall()

        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(HEADERS)
        for r in rows:
            writer.writerow([r[h] for h in HEADERS])
        output.seek(0)

        return StreamingResponse(output, media_type="text/csv",
            headers={"Content-Disposition": "attachment; filename=snp_data.csv"})
    finally:
        con.close()
