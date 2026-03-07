"""
Figure 1 Panels A & B -- Network Visualisation v5
==================================================
Creative solution to the "empty canvas + collapsed blobs" problem:

1. Community detection at |rho| > 0.8 (biological structure)
2. Layout computed at |rho| > 0.9 (fewer edges = nodes spread out)
3. Post-layout expansion: each community cluster is expanded
   proportional to sqrt(size) so they fill the canvas
4. Normalise to fill [-1,1] range
5. No stats bar (avoids contradicting manuscript)
6. Big nodes, community-coloured edges at low alpha

Author: Plant Systems Biology Lab | 2025 Revision
"""

import os
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgba
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")

# =============================================================================
# PATHS — USER CONFIGURATION: update these paths for your local environment
# =============================================================================
leaf_data_path = r"C:\Users\ms\Desktop\data_chem_3_10\data\data\n_p_l.csv"
root_data_path = r"C:\Users\ms\Desktop\data_chem_3_10\data\data\n_p_r.csv"
vip_file_path  = r"C:\Users\ms\Desktop\data_chem_3_10\output\results\spearman\VIP_mann_whitney_bonferroni_fdr_combine_above_one.csv"
output_dir     = r"C:\Users\ms\Desktop\data_chem_3_10\output\revision_pack\figures"
os.makedirs(output_dir, exist_ok=True)

vip_df = pd.read_csv(vip_file_path)
vip_set = set(vip_df["Metabolite"].tolist())
print(f"[OK] {len(vip_set)} VIP metabolites")

# =============================================================================
# COLOURS
# =============================================================================
COMM_COLS = [
    "#17AB81",   # dark green anchor
    "#00CED1",   # dark turquoise (bright blue-green)
    "#B8DC6F",   # yellow-green
    "#0EADAD",   # teal
    "#00BFFF",   # deep sky blue
    "#57D590",   # medium sea green
    "#4682B4",   # steel blue
    "#94CB64",   # lime
    "#08937E",   # nature green
    "#48D1CC",   # med turquoise
    "#2E8B57",   # dark sea green
    "#24B34F",   # spring green
    "#6495ED",   # cornflower
    "#307458",   # dark forest
    "#98FB98",   # pale green
]
VIP_RING = "#1a1a1a"


# =============================================================================
# BUILD GRAPHS (two thresholds)
# =============================================================================
def load_corr_matrix(data_path):
    """Load data, filter Treatment==1, compute Spearman correlation."""
    data = pd.read_csv(data_path)
    data = data[data["Treatment"] == 1]
    metab = data[data.columns[10:]]
    corr = metab.corr(method="spearman")
    return corr


def build_graph(corr, rho_thresh):
    """Build graph from correlation matrix at given threshold."""
    G = nx.Graph()
    nodes = corr.columns.tolist()
    G.add_nodes_from(nodes)

    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            rho = corr.iloc[i, j]
            if pd.notna(rho) and abs(rho) > rho_thresh:
                G.add_edge(nodes[i], nodes[j], weight=rho)

    isolates = list(nx.isolates(G))
    G.remove_nodes_from(isolates)
    print(f"  |rho|>{rho_thresh}: {G.number_of_nodes()} nodes, "
          f"{G.number_of_edges()} edges ({len(isolates)} isolates removed)")
    return G


# =============================================================================
# COMMUNITIES
# =============================================================================
def detect_communities(G):
    """Louvain community detection."""
    try:
        partition = nx.community.louvain_communities(G, seed=42, resolution=1.0)
    except AttributeError:
        try:
            import community as community_louvain
            pd_dict = community_louvain.best_partition(G, random_state=42)
            comms = {}
            for n, c in pd_dict.items():
                comms.setdefault(c, set()).add(n)
            partition = list(comms.values())
        except ImportError:
            partition = list(nx.community.greedy_modularity_communities(G))

    partition = sorted(partition, key=len, reverse=True)

    node_to_comm = {}
    for cid, members in enumerate(partition):
        for node in members:
            node_to_comm[node] = cid

    sizes = [len(c) for c in partition]
    print(f"  Communities: {len(partition)}, sizes: {sizes[:8]}")
    return partition, node_to_comm


# =============================================================================
# LAYOUT: spring on sparse graph + post-expansion
# =============================================================================
def compute_layout(G_layout, G_viz, partition, node_to_comm, seed=42):
    """
    1. Spring layout on G_layout (high-threshold, sparse)
    2. Assign positions to G_viz nodes not in G_layout via community centroid + jitter
    3. Expand each community cluster proportional to sqrt(size)
    4. Normalise to fill canvas
    """
    n = G_layout.number_of_nodes()
    k = 6.0 / np.sqrt(max(n, 1))
    print(f"  Spring layout on sparse graph: k={k:.3f}, 600 iterations...")
    pos_layout = nx.spring_layout(G_layout, k=k, iterations=600, seed=seed,
                                   weight=None)

    rng = np.random.RandomState(seed)

    # All nodes we want to draw (from G_viz)
    pos = {}

    # Copy positions from layout graph
    for node in G_viz.nodes():
        if node in pos_layout:
            pos[node] = np.array(pos_layout[node], dtype=float)

    # Nodes in G_viz but NOT in G_layout: place near their community centroid
    missing = [n for n in G_viz.nodes() if n not in pos]
    if missing:
        comm_coords = {}
        for node, xy in pos.items():
            cid = node_to_comm.get(node, -1)
            if cid >= 0:
                comm_coords.setdefault(cid, []).append(xy)

        comm_centroids = {}
        for cid, coords in comm_coords.items():
            comm_centroids[cid] = np.mean(coords, axis=0)

        global_centroid = np.mean(list(pos.values()), axis=0)
        for node in missing:
            cid = node_to_comm.get(node, -1)
            if cid in comm_centroids:
                center = comm_centroids[cid]
            else:
                center = global_centroid
            jitter = rng.randn(2) * 0.03
            pos[node] = center + jitter

    print(f"  Positioned: {len(pos)} nodes ({len(missing)} via community centroids)")

    # --- Post-layout expansion ---
    max_comm_size = max(len(c) for c in partition) if partition else 1

    for cid, members in enumerate(partition):
        members_in_pos = [n for n in members if n in pos]
        if len(members_in_pos) < 2:
            continue

        coords = np.array([pos[n] for n in members_in_pos])
        centroid = coords.mean(axis=0)

        # Bigger communities expand more
        expansion = 0.5 + 1.5 * np.sqrt(len(members_in_pos) / max_comm_size)

        for n in members_in_pos:
            pos[n] = centroid + (pos[n] - centroid) * expansion

    # --- Normalise to fill [-1, 1] ---
    all_coords = np.array(list(pos.values()))
    for dim in range(2):
        lo, hi = all_coords[:, dim].min(), all_coords[:, dim].max()
        span = hi - lo
        if span > 0:
            for node in pos:
                pos[node][dim] = 2.0 * (pos[node][dim] - lo) / span - 1.0

    # --- Post-layout repulsion: push overlapping nodes apart ---
    # Vectorised: only affects nodes closer than min_dist
    nodes_list = list(pos.keys())
    coords = np.array([pos[n] for n in nodes_list])
    min_dist = 0.04

    for iteration in range(30):
        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=min_dist)
        if not pairs:
            break
        for i, j in pairs:
            diff = coords[j] - coords[i]
            dist = np.linalg.norm(diff)
            if dist > 0:
                push = (min_dist - dist) / 2.0 * 0.6
                direction = diff / dist
                coords[i] -= direction * push
                coords[j] += direction * push

    # Re-normalise after repulsion
    for dim in range(2):
        lo, hi = coords[:, dim].min(), coords[:, dim].max()
        span = hi - lo
        if span > 0:
            coords[:, dim] = 2.0 * (coords[:, dim] - lo) / span - 1.0

    for i, node in enumerate(nodes_list):
        pos[node] = coords[i]

    print(f"  Repulsion: {iteration + 1} iterations, min_dist={min_dist}")

    return pos


# =============================================================================
# SAMPLE HELPER
# =============================================================================
def sample_idx(n, k, seed=42):
    rng = np.random.RandomState(seed)
    return rng.choice(n, size=min(k, n), replace=False)


# =============================================================================
# PLOT
# =============================================================================
def plot_network(G, partition, node_to_comm, pos, vip_set,
                 title, output_path, max_vis_edges=15000):

    fig, ax = plt.subplots(1, 1, figsize=(7, 4), facecolor="white")
    ax.set_facecolor("white")

    degrees = dict(G.degree())
    max_deg = max(degrees.values()) if degrees else 1
    n_comms = len(partition)

    # --- Edges ---
    intra_segs, intra_cols = [], []
    inter_segs = []

    for u, v, d in G.edges(data=True):
        if u not in pos or v not in pos:
            continue
        cu = node_to_comm.get(u, -1)
        cv = node_to_comm.get(v, -1)
        seg = [pos[u], pos[v]]

        if cu == cv and cu >= 0:
            col = COMM_COLS[cu % len(COMM_COLS)]
            intra_segs.append(seg)
            intra_cols.append(to_rgba(col, 0.05))
        else:
            inter_segs.append(seg)

    # Sample down if too many
    if len(intra_segs) > max_vis_edges:
        idx = sample_idx(len(intra_segs), max_vis_edges)
        intra_segs = [intra_segs[i] for i in idx]
        intra_cols = [intra_cols[i] for i in idx]
        print(f"  Sampled intra-edges to {max_vis_edges}")

    if len(inter_segs) > max_vis_edges // 4:
        limit = max_vis_edges // 4
        idx = sample_idx(len(inter_segs), limit)
        inter_segs = [inter_segs[i] for i in idx]
        print(f"  Sampled inter-edges to {limit}")

    # Draw inter (ghost)
    if inter_segs:
        lc = LineCollection(inter_segs, colors=to_rgba("#aaaaaa", 0.012),
                           linewidths=0.08, zorder=1)
        ax.add_collection(lc)

    # Draw intra (coloured)
    if intra_segs:
        lc = LineCollection(intra_segs, colors=intra_cols,
                           linewidths=0.27, zorder=2)
        ax.add_collection(lc)

    print(f"  Edges drawn: {len(intra_segs)} intra, {len(inter_segs)} inter")

    # --- Nodes (smallest communities first, biggest on top) ---
    for cid in reversed(range(n_comms)):
        members = [n for n in partition[cid] if n in pos]
        if not members:
            continue

        col = COMM_COLS[cid % len(COMM_COLS)]
        xs = np.array([pos[n][0] for n in members])
        ys = np.array([pos[n][1] for n in members])
        degs = np.array([degrees.get(n, 0) for n in members])

        # Aggressive sizing scaled for 6x7 canvas
        sizes = 10 + 149 * np.sqrt(degs / max_deg)

        vip_mask = np.array([n in vip_set for n in members])

        # Non-VIP
        if (~vip_mask).any():
            ax.scatter(xs[~vip_mask], ys[~vip_mask], s=sizes[~vip_mask],
                      c=col, alpha=0.7, edgecolors=to_rgba(col, 0.3),
                      linewidths=0.3, zorder=3)

        # VIP (black ring)
        if vip_mask.any():
            ax.scatter(xs[vip_mask], ys[vip_mask], s=sizes[vip_mask] * 1.4,
                      c=col, alpha=0.85, edgecolors=VIP_RING,
                      linewidths=0.7, zorder=4)

    # --- Community labels (>= 10 members) ---
    for cid, members in enumerate(partition):
        members_in = [n for n in members if n in pos]
        if len(members_in) < 10:
            continue
        coords = np.array([pos[n] for n in members_in])
        cx, cy = np.median(coords, axis=0)

        n_vip = sum(1 for n in members_in if n in vip_set)
        label = f"M{cid+1} (n={len(members_in)})"
        if n_vip > 0:
            label += f"\n{n_vip} VIP"

        ax.annotate(label, (cx, cy), ha="center", va="center",
                    fontsize=6.5, fontweight="bold", color="#2c3e50",
                    alpha=0.8,
                    bbox=dict(boxstyle="round,pad=0.2",
                              facecolor="white", edgecolor="#cccccc",
                              alpha=0.65, linewidth=0.4),
                    zorder=6)

    # --- No title (other Fig 1 panels don't have one) ---

    # --- Legend ---
    handles = []
    for cid in range(min(8, n_comms)):
        col = COMM_COLS[cid % len(COMM_COLS)]
        sz = len(partition[cid])
        if sz < 5:
            continue
        handles.append(
            mpatches.Patch(facecolor=col, edgecolor="#808080", linewidth=0.5,
                          label=f"Module {cid+1} (n={sz})")
        )
    if n_comms > 8:
        remaining = sum(1 for c in partition[8:] if len(c) >= 5)
        if remaining > 0:
            handles.append(
                mpatches.Patch(facecolor="#b3b3b3", edgecolor="#808080",
                              label=f"+ {remaining} more")
            )
    handles.append(
        mlines.Line2D([], [], color=VIP_RING, marker="o", linestyle="None",
                      markersize=5, markerfacecolor="white",
                      markeredgewidth=1.0, label="VIP feature")
    )

    leg = ax.legend(handles=handles, loc="lower center",
                    bbox_to_anchor=(0.5, 1.02),
                    ncol=3, fontsize=9, framealpha=0.9,
                    edgecolor="#b3b3b3", fancybox=True,
                    borderpad=0.5, labelspacing=0.4,
                    columnspacing=0.8, handletextpad=0.4,
                    handlelength=1.2)
    leg.get_frame().set_linewidth(0.4)

    # --- Axis (NO stats bar) ---
    ax.axis("off")
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)

    plt.subplots_adjust(bottom=0.12)
    plt.savefig(output_path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    pdf_path = output_path.replace(".png", ".pdf")
    plt.savefig(pdf_path, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close()
    print(f"  [OK] {output_path}")
    print(f"  [OK] {pdf_path}")


# =============================================================================
# MAIN
# =============================================================================
def process_tissue(data_path, tissue, filename, title):
    print(f"\n{'='*50}")
    print(f"  {tissue}")
    print(f"{'='*50}")

    corr = load_corr_matrix(data_path)
    print(f"  Correlation matrix: {corr.shape[0]} features")

    # Viz graph at 0.8 (community detection + drawing)
    G_viz = build_graph(corr, rho_thresh=0.8)

    # Layout graph at 0.93 (very sparse = hubs spread out more)
    G_layout = build_graph(corr, rho_thresh=0.93)

    # Communities from viz graph
    partition, node_to_comm = detect_communities(G_viz)

    # Layout from sparse graph + expand + normalise
    pos = compute_layout(G_layout, G_viz, partition, node_to_comm, seed=42)

    out = os.path.join(output_dir, filename)
    plot_network(G_viz, partition, node_to_comm, pos, vip_set,
                 title=title, output_path=out)


print("\nFigure 1 A & B -- Network v5")
print("=" * 60)

process_tissue(leaf_data_path, "Leaf",
               "Fig1A_Leaf_Network_Community.png",
               "A.  Leaf Metabolic Network")

process_tissue(root_data_path, "Root",
               "Fig1B_Root_Network_Community.png",
               "B.  Root Metabolic Network")

print(f"\n[DONE] {output_dir}")