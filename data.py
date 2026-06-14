"""
Generate random Erdős-Rényi directed graphs used in DST experiments.

uv sync
.venv/bin/python data.py --weighted both          # all standard sizes
.venv/bin/python data.py --sizes 10000000 --weighted true   # 10M only
.venv/bin/python data.py --sizes 100000000 --weighted true  # 100M only (~13 min)
"""

from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import networkx as nx
import numpy as np

N_NEIGHBOR = 10
DATASETS_DIR = Path("datasets")


def iter_fast_gnp_directed_edges(n: int, p: float, seed: int = 42):
    """Match networkx.fast_gnp_random_graph(..., directed=True) edge generation."""
    rng = random.Random(seed)
    lp = math.log(1.0 - p)

    v = 1
    w = -1
    while v < n:
        lr = math.log(1.0 - rng.random())
        w = w + 1 + int(lr / lp)
        while w >= v and v < n:
            w = w - v
            v = v + 1
        if v < n:
            yield w, v

    v = 1
    w = -1
    while v < n:
        lr = math.log(1.0 - rng.random())
        w = w + 1 + int(lr / lp)
        while w >= v and v < n:
            w = w - v
            v = v + 1
        if v < n:
            yield v, w


def write_graph_networkx(
    n: int,
    n_neighbor: int,
    weighted: bool,
    out_path: Path,
    seed: int = 42,
) -> int:
    p = n_neighbor / n
    g = nx.fast_gnp_random_graph(n, p, directed=True, seed=seed)

    np.random.seed(seed)
    edge_count = 0
    with out_path.open("w") as fout:
        for v, nbrs in g.adj.items():
            if weighted:
                weights = np.random.rand(len(nbrs))
            else:
                weights = [1] * len(nbrs)
            for nbr, w in zip(nbrs, weights):
                fout.write(f"{v},{nbr},{w}\n")
                edge_count += 1
    return edge_count


def write_graph_streaming(
    n: int,
    n_neighbor: int,
    weighted: bool,
    out_path: Path,
    seed: int = 42,
) -> int:
    p = n_neighbor / n
    rng = np.random.default_rng(seed)

    edge_count = 0
    with out_path.open("w") as fout:
        for src, dst in iter_fast_gnp_directed_edges(n, p, seed=seed):
            weight = rng.random() if weighted else 1
            fout.write(f"{src},{dst},{weight}\n")
            edge_count += 1
    return edge_count


def generate_graph(
    n: int,
    weighted: bool,
    n_neighbor: int = N_NEIGHBOR,
    seed: int = 42,
    streaming_threshold: int = 1_000_000,
) -> Path:
    DATASETS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = DATASETS_DIR / f"random_graph_w{weighted}_nghb{n_neighbor}_n{n}.csv"

    print(f"Generating {out_path.name} (n={n}, weighted={weighted}) ...", flush=True)
    if n <= streaming_threshold:
        edge_count = write_graph_networkx(n, n_neighbor, weighted, out_path, seed=seed)
    else: # avoids OOM on 30 GB RAM
        edge_count = write_graph_streaming(n, n_neighbor, weighted, out_path, seed=seed)

    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(
        f"  done: {edge_count:,} edges, {size_mb:.1f} MiB -> {out_path}",
        flush=True,
    )
    return out_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sizes",
        type=int,
        nargs="+",
        default=[1_000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000],
    )
    parser.add_argument(
        "--weighted",
        choices=("true", "false", "both"),
        default="both",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    weighted_values: list[bool]
    if args.weighted == "both":
        weighted_values = [False, True]
    elif args.weighted == "true":
        weighted_values = [True]
    else:
        weighted_values = [False]

    for weighted in weighted_values:
        for n in args.sizes:
            if not weighted and n > 1_000_000:
                continue
            generate_graph(n, weighted, seed=args.seed)


if __name__ == "__main__":
    main()
