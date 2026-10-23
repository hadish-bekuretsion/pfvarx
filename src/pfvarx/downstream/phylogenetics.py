"""Simple phylogenetic utilities.

A Neighbor-Joining tree is built from a pairwise distance matrix in pure
Python.  IQ-TREE 2 is available as an optional wrapper for maximum-likelihood
inference; the NJ tree is intended for rapid QC and visualisation only.
"""

from __future__ import annotations

from dataclasses import dataclass

from pfvarx.downstream.relatedness import RelatednessResult


@dataclass
class TreeNode:
    name: str | None = None
    children: list[TreeNode] | None = None
    branch_length: float = 0.0

    def to_newick(self) -> str:
        if self.children is None:
            return f"{self.name}:{self.branch_length:.6f}"
        inner = ",".join(c.to_newick() for c in self.children)
        label = self.name or ""
        return f"({inner}){label}:{self.branch_length:.6f}"


def distance_matrix(results: list[RelatednessResult], samples: list[str]) -> dict[tuple[str, str], float]:
    """Distance = 1 - IBS sharing."""
    dist = {(s, s): 0.0 for s in samples}
    for r in results:
        d = 1.0 - r.relatedness_estimate
        dist[(r.sample_a, r.sample_b)] = d
        dist[(r.sample_b, r.sample_a)] = d
    return dist


def neighbor_joining(samples: list[str], dist: dict[tuple[str, str], float]) -> TreeNode:
    """Classic Neighbor-Joining (Saitou and Nei, 1987) with a simple tree."""
    nodes: dict[str, TreeNode] = {s: TreeNode(name=s) for s in samples}
    active = list(samples)
    d = dict(dist)

    while len(active) > 2:
        n = len(active)
        r = {i: sum(d[(i, j)] for j in active if j != i) for i in active}
        best = None
        for i in range(n):
            for j in range(i + 1, n):
                a, b = active[i], active[j]
                q = (n - 2) * d[(a, b)] - r[a] - r[b]
                if best is None or q < best[0]:
                    best = (q, a, b)
        _, a, b = best
        dab = d[(a, b)]
        la = 0.5 * dab + (r[a] - r[b]) / (2 * (n - 2))
        lb = dab - la
        parent = TreeNode(children=[nodes[a], nodes[b]])
        nodes[a].branch_length = max(0.0, la)
        nodes[b].branch_length = max(0.0, lb)
        new_name = f"({a},{b})"
        nodes[new_name] = parent
        for k in active:
            if k in (a, b):
                continue
            dk = 0.5 * (d[(a, k)] + d[(b, k)] - dab)
            d[(new_name, k)] = dk
            d[(k, new_name)] = dk
        active = [x for x in active if x not in (a, b)] + [new_name]

    if len(active) == 2:
        a, b = active
        root = TreeNode(children=[nodes[a], nodes[b]])
        nodes[a].branch_length = d[(a, b)] / 2
        nodes[b].branch_length = d[(a, b)] / 2
        return root
    return nodes[active[0]]


def newick(samples: list[str], results: list[RelatednessResult]) -> str:
    dist = distance_matrix(results, samples)
    tree = neighbor_joining(samples, dist)
    return tree.to_newick() + ";"
