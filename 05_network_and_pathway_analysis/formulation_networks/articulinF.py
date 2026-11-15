import matplotlib.pyplot as plt
import networkx as nx

# Initialize the graph
G = nx.Graph()

# Compound nodes
compounds = {
    "Withaferin A": "gold",
    "Guggulsterone-E": "gold",
    "Guggulsterone-Z": "gold",
    "Curcumin": "gold",
    "Demethoxycurcumin": "gold"
}

# First-level targets for compounds
compound_targets = {
    "Withaferin A": ["PTGS2", "IL6", "IL1B", "CCL2", "CXCL8", "JUN", "MAPK1", "MAPK14", "TNF", "NFE2L2", "RELA", "VEGFA"],
    "Guggulsterone-E": ["IL1B", "MMP9", "PTGS2", "RELA", "TNF", "NFKB1", "CASP3", "TGFB1", "VEGFA", "ALOX5", "CTSL", "CTSK", "NOS2"],
    "Guggulsterone-Z": ["ALOX5", "CTSK", "IKBKB", "NOS2"],
    "Curcumin": ["ALOX5", "IL6", "IL1B", "IL18", "TLR2", "CCL2", "CCL3", "CXCL1", "CXCL12", "MMP3", "MMP9", "CTSL", "PTGS2", "JUN", "STAT3", "MAPK1", "MAPK14", "RELA", "TNF", "NFKB1", "AKT1", "CASP3", "CAT", "CXCL8", "HIF1A", "IFNG", "IKBKB", "NFE2L2", "NOS2", "SOD1", "TGFB1", "VEGFA", "ALOX5AP"],
    "Demethoxycurcumin": ["ALOX5", "MMP3", "CTSL", "ALOX5AP", "STAT3", "TNF", "AKT1", "NFE2L2", "NOS2", "MMP9", "MAPK1", "BAX", "CASP3", "VEGFA"]
}

# Secondary interactors
secondary_interactions = {
    "IL6": ["IL1B", "PTGS2", "VEGFA", "JAK2", "JAK1", "IL2", "IL3"],
    "IL1B": ["VEGFA", "PTGS2", "IL6", "IL5", "IL2", "IL4", "IFNG", "MYD88"],
    "VEGFA": ["IL1B", "IL6", "CDKN2A"],
    "PTGS2": ["IL1B", "IL6", "CXCL8", "TP53", "ALOX5", "ALOX15", "TNF", "NOS2"],
    "CCL2": ["CXCL8", "CCL5", "CCL8", "CXCL13", "CCR1", "CCL4L2", "CCL15", "CCR2", "CXCL9", "CCR5"],
    "CXCL8": ["IL4", "CCL8", "CXCR3", "CXCR1", "CCR1", "CCR2", "PF4"],
    "RELA": ["CREBBP", "EP300", "NFKB1A", "NFKB2"],
    "TNF": ["CHUK", "TRADD", "NFKB1", "TRAF2", "TAB2"],
    "JUN": ["ATF2", "FOS", "EP300", "JUND", "JUNB"],
    "MAPK1": ["STAT3", "TP53", "MAP2K2", "MAP2K1"],
    "MAPK14": ["JUN", "ATF2", "TP53", "TAB1"],
    "ALOX5": ["PTGS2", "ALOX5AP", "LTA4H"],
    "MMP9": ["TNF", "TLR4", "IL1B", "MMP3", "TGFB1"],
    "NFKB1": ["IL1A", "CXCL1", "CXCL2", "IL18", "TLR4", "IL1B", "TNF", "RELA"],
    "TGFB1": ["IL10", "IL1B", "TNF", "MMP9"],
    "CASP3": ["TNF", "CYCS", "AKT1", "MAPK1"],
    "CYCS": ["BCL2"],
    "CTSK": ["CTSL", "TLR4", "ACP5", "TNFSF11"],
    "IKBKB": ["CHUK", "RELA", "NFKB1", "TNF", "TRAF6", "TNFRSF1A"],
    "IL18": ["IL1A", "IFNG", "JUN", "NFKB1", "RELA", "IL6", "CXCL8"],
    "SOD1": ["CAT"],
    "HIF1A": ["JUN", "MAPK1", "RELA", "STAT3"],
    "NFE2L2": ["JUN", "AKT1"]
}

# Add compound nodes
for compound, color in compounds.items():
    G.add_node(compound, color=color, size=800)

# Add first-level targets and edges
all_targets = set()
for compound, targets in compound_targets.items():
    for target in targets:
        G.add_node(target, color="skyblue", size=600)
        G.add_edge(compound, target)
        all_targets.add(target)

# Add secondary interactions
for target, secondaries in secondary_interactions.items():
    for interactor in secondaries:
        G.add_node(interactor, color="lightgreen", size=400)
        G.add_edge(target, interactor)

# Identify shared first-level targets
from collections import Counter
target_counts = Counter([t for targets in compound_targets.values() for t in targets])
shared_targets = {t for t, count in target_counts.items() if count > 1}

# Update shared targets to red and increase size
for target in shared_targets:
    if target in G.nodes:
        G.nodes[target]["color"] = "red"
        G.nodes[target]["size"] = 900

# Draw the graph
pos = nx.spring_layout(G, k=0.4, iterations=200)
node_colors = [G.nodes[n]["color"] for n in G.nodes]
node_sizes = [G.nodes[n]["size"] for n in G.nodes]

plt.figure(figsize=(24, 24))
nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=node_sizes, edge_color="gray", font_size=8)
plt.title("Articulin F: Compound-Target and Secondary Interactor Network", fontsize=20)
plt.show()

