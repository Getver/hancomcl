# -*- coding: utf-8 -*-
import sqlite3, csv, os, time

DB_PATH = "/disk0/yang93/eQTL/data/eqtl_idmap.db"
DATA_DIR = "/disk0/yang93/eQTL/data"

FILES = [
    ("1_CHR_ALL_T2D.txt", "T2D"),
    ("2_CHR_ALL_Non_T2D.txt", "Non_T2D"),
    ("3_CHR_ALL_Only_DR.txt", "Only_DR"),
    ("4_CHR_ALL_R.txt", "R_All"),
    ("5_DR_0.txt", "DR_0"),
    ("6_DR_5.txt", "DR_5"),
    ("7_Normal_0.txt", "Normal_0"),
    ("8_Normal_5.txt", "Normal_5"),
]

os.makedirs(os.path.dirname(DB_PATH), exist_ok=True)
conn = sqlite3.connect(DB_PATH)
cur = conn.cursor()

# 성능 최적화 옵션
cur.execute("PRAGMA journal_mode=WAL;")
cur.execute("PRAGMA synchronous=OFF;")     # 안전성 조금 줄이고 속도 ↑
cur.execute("PRAGMA cache_size=-200000;")
cur.execute("PRAGMA temp_store=MEMORY;")

# 스키마 생성
cur.executescript("""
DROP TABLE IF EXISTS SNP;
DROP TABLE IF EXISTS Gene;
DROP TABLE IF EXISTS Source;
DROP TABLE IF EXISTS eQTL;

CREATE TABLE SNP (
  id   INTEGER PRIMARY KEY,
  name TEXT UNIQUE
);
CREATE TABLE Gene (
  id     INTEGER PRIMARY KEY,
  symbol TEXT UNIQUE
);
CREATE TABLE Source (
  id    INTEGER PRIMARY KEY,
  label TEXT UNIQUE
);
CREATE TABLE eQTL (
    id        INTEGER PRIMARY KEY AUTOINCREMENT,
    chr       INTEGER,
    bp        INTEGER,
    ref       TEXT,
    alt       TEXT,
    snp_id    INTEGER,
    gene_id   INTEGER,
    distance  INTEGER,
    npval     REAL,
    slope     REAL,
    fdr       REAL,
    source_id INTEGER
);
""")

snp_map, gene_map, source_map = {}, {}, {}
snp_next, gene_next, source_next = 1, 1, 1

def get_or_create(mapping, key, counter_name):
    global snp_next, gene_next, source_next
    if key in mapping:
        return mapping[key]
    if counter_name == "snp":
        mapping[key] = snp_next; snp_next += 1
    elif counter_name == "gene":
        mapping[key] = gene_next; gene_next += 1
    elif counter_name == "source":
        mapping[key] = source_next; source_next += 1
    return mapping[key]

t_all = time.time()
total_rows = 0

for filename, label in FILES:
    path = os.path.join(DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"[SKIP] 파일 없음: {path}")
        continue

    print(f"[LOAD] {path} (label={label})")
    t0 = time.time()
    source_id = get_or_create(source_map, label, "source")

    conn.execute("BEGIN;")
    rows, count = [], 0

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                fdr_value = row.get("FDR") or row.get("nPval_FDR")
                snp_id = get_or_create(snp_map, row["SNP"], "snp")
                gene_id = get_or_create(gene_map, row["GENE"], "gene")
                rows.append((
                    int(row["CHR"]),
                    int(row["BP"]),
                    row["REF"], row["ALT"],
                    snp_id, gene_id,
                    int(row["DISTANCE"]),
                    float(row["nPval"]),
                    float(row["SLOPE"]),
                    float(fdr_value) if fdr_value else None,
                    source_id
                ))
                count += 1
                # --- 청크 단위 insert ---
                if len(rows) >= 50000:
                    cur.executemany("""
                        INSERT INTO eQTL(chr,bp,ref,alt,snp_id,gene_id,distance,npval,slope,fdr,source_id)
                        VALUES (?,?,?,?,?,?,?,?,?,?,?)
                    """, rows)
                    rows.clear()
            except Exception as e:
                print(f"[WARN] {filename} line {count}: {e}")

    if rows:
        cur.executemany("""
            INSERT INTO eQTL(chr,bp,ref,alt,snp_id,gene_id,distance,npval,slope,fdr,source_id)
            VALUES (?,?,?,?,?,?,?,?,?,?,?)
        """, rows)

    conn.commit()
    total_rows += count
    elapsed = time.time() - t0
    print(f"[DONE] {filename}: {count:,} rows inserted, elapsed {elapsed/60:.1f} min")

# 매핑 테이블 insert
cur.executemany("INSERT OR IGNORE INTO SNP(id,name) VALUES (?,?)", [(v,k) for k,v in snp_map.items()])
cur.executemany("INSERT OR IGNORE INTO Gene(id,symbol) VALUES (?,?)", [(v,k) for k,v in gene_map.items()])
cur.executemany("INSERT OR IGNORE INTO Source(id,label) VALUES (?,?)", [(v,k) for k,v in source_map.items()])
conn.commit()
conn.close()

print(f"\nAll done. Total rows: {total_rows:,} | Elapsed {(time.time()-t_all)/60:.1f} min")
