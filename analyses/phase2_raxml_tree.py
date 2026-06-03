"""
L4.13 Analysis -- Phase 2: RAxML-NG tree with bootstrap
========================================================

Ingroup: 258 L4.13 strains (61 L4.13.1 + 197 L4.13.2)
Outgroup: 5 L4.4.2 (direct sister of L4.13 per topology
          (L4.1, (L4.2, (L4.13, (L4.4, ...))))
          + 5 L4.2.1 (more distal) for a stable root

Model: BIN+G (binary 0/1 SNP matrix + Gamma site-rate variation)
Bootstrap: 100 (Felsenstein)
Seed: 42

Output:
    résultats/phase2_tree/alignment.fasta
    résultats/phase2_tree/strains.tsv
    résultats/phase2_tree/L4.13.raxml.support  (final tree)

Usage:
    python phase2_raxml_tree.py
"""

from __future__ import annotations

import random
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BDD = Path(__file__).resolve().parent.parent / "external" / "spdi_db"  # per-strain SPDI profiles; see external/README.md
RES = ROOT / "résultats"
TREE_DIR = RES / "phase2_tree"
TREE_DIR.mkdir(parents=True, exist_ok=True)

RAXML = "raxml-ng"  # RAxML-NG on PATH; https://github.com/amkozlov/raxml-ng

GROUPS = {
    "L4.13.1": BDD / "L4.13.1",
    "L4.13.2": BDD / "L4.13.2",
}
OUTGROUPS = {
    "L4.4.2": BDD / "L4.4.2",
    "L4.2.1": BDD / "L4.2.1",
}
OUTGROUP_SAMPLE_PER_SUBLIN = 5
SEED = 42


def list_strains(lineage_dir: Path) -> list[tuple[str, Path]]:
    out = []
    if not lineage_dir.is_dir():
        return out
    for sra_dir in sorted(lineage_dir.iterdir()):
        spdi_file = sra_dir / "NC_000962.3" / "spdi.txt"
        if spdi_file.is_file():
            out.append((sra_dir.name, spdi_file))
    return out


def load_spdis(spdi_file: Path) -> set[str]:
    with spdi_file.open() as f:
        return {line.strip() for line in f if line.strip()}


def main():
    rng = random.Random(SEED)

    print("=" * 70)
    print("L4.13 Analysis -- Phase 2: RAxML-NG tree")
    print("=" * 70)

    # ---- 1. Load ingroup ----
    strains = []
    for g, p in GROUPS.items():
        for sra, sf in list_strains(p):
            strains.append((f"{g}_{sra}", sra, g, load_spdis(sf)))
    print(f"\n[load] ingroup: {len(strains)} strains")

    # ---- 2. Load outgroup ----
    outgroup_labels = []
    for g, p in OUTGROUPS.items():
        sub = list_strains(p)
        rng.shuffle(sub)
        kept = sub[:OUTGROUP_SAMPLE_PER_SUBLIN]
        print(f"[load] outgroup {g}: {len(kept)} strains (of {len(sub)})")
        for sra, sf in kept:
            lbl = f"{g}_{sra}"
            strains.append((lbl, sra, g, load_spdis(sf)))
            outgroup_labels.append(lbl)
    print(f"[total] {len(strains)} strains "
          f"({len(strains) - len(outgroup_labels)} ingroup + "
          f"{len(outgroup_labels)} outgroup)")

    # ---- 3. Pan-SPDI + filter parsimony-informative ----
    all_spdis = set()
    for _, _, _, s in strains:
        all_spdis |= s
    print(f"\n[pan-SPDI] {len(all_spdis)} distinct SPDIs")

    counts = {s: 0 for s in all_spdis}
    for _, _, _, sset in strains:
        for s in sset:
            counts[s] += 1

    n = len(strains)
    variable = [s for s, c in counts.items() if 0 < c < n]
    pan = sorted(s for s in variable if 2 <= counts[s] <= n - 2)
    print(f"[filter] variable sites: {len(variable)}")
    print(f"[filter] parsimony-informative sites: {len(pan)}")

    # ---- 4. Write binary FASTA ----
    fasta = TREE_DIR / "alignment.fasta"
    with fasta.open("w") as f:
        for label, _, _, sset in strains:
            seq = "".join("1" if s in sset else "0" for s in pan)
            f.write(f">{label}\n{seq}\n")
    print(f"\n[write] {fasta} ({fasta.stat().st_size / 1e6:.2f} MB)")

    # ---- 5. Outgroup and strain tables ----
    (TREE_DIR / "outgroup.txt").write_text(",".join(outgroup_labels) + "\n")
    (TREE_DIR / "strains.tsv").write_text(
        "label\tsra\tgroup\n"
        + "\n".join(f"{l}\t{s}\t{g}" for l, s, g, _ in strains)
        + "\n"
    )
    (TREE_DIR / "pan_spdi.txt").write_text("\n".join(pan) + "\n")

    # ---- 6. Launch RAxML-NG ----
    prefix = TREE_DIR / "L4.13"
    outgroup_str = ",".join(outgroup_labels)
    cmd = [
        str(RAXML),
        "--all",
        "--msa", str(fasta),
        "--model", "BIN+G",
        "--prefix", str(prefix),
        "--seed", str(SEED),
        "--bs-trees", "100",
        "--threads", "auto{16}",
        "--outgroup", outgroup_str,
        "--redo",
    ]
    (TREE_DIR / "raxml_cmd.sh").write_text(" ".join(cmd) + "\n")
    log_file = TREE_DIR / "raxml_stdout.log"
    print(f"\n[raxml] launching...")
    print(f"  command: {RAXML.name} --all --msa alignment.fasta --model BIN+G "
          f"--prefix {prefix.name} --seed {SEED} --bs-trees 100 "
          f"--threads auto{{16}} --outgroup <{len(outgroup_labels)} outgroup>")
    with log_file.open("w") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    print(f"[raxml] exit code: {proc.returncode}")
    print(f"[raxml] log: {log_file}")

    support_tree = Path(str(prefix) + ".raxml.support")
    if support_tree.is_file():
        print(f"\n[OK] final tree with bootstrap support: {support_tree}")
    else:
        print(f"\n[warn] {support_tree.name} not found; check log")

    print("\n" + "=" * 70)
    print("Phase 2 complete.")
    print("=" * 70)


if __name__ == "__main__":
    main()
