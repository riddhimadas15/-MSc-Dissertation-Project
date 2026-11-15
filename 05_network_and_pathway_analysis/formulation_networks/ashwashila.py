import networkx as nx
import matplotlib.pyplot as plt

# Define the graph
G = nx.Graph()

# Define the compounds and their direct targets
compounds = {
    'Withaferin A': ['IL6', 'IL1B', 'PTGS2', 'CCL2', 'CXCL8', 'JUN', 'MAPK1', 'MAPK14', 'TNF', 'NFE2L2', 'RELA', 'IKBKB', 'VEGFA'],
    'Fulvic Acid': ['PTGS2', 'MAPK1', 'RELA', 'CXCL8', 'VEGFA'],
    'Humic Acid': ['MMP1', 'NOS2', 'IL1B', 'VEGFA', 'MAPK1', 'RORC', 'TLR4']
}

# Add compound nodes
for compound in compounds:
    G.add_node(compound)

# Map to track target appearances
from collections import defaultdict
target_map = defaultdict(set)
for compound, targets in compounds.items():
    for target in targets:
        G.add_edge(compound, target, connection='direct')
        target_map[target].add(compound)

# Secondary interaction partners
interactions = {
    'IL6': ['IL1B', 'PTGS2', 'VEGFA', 'JAK2', 'JAK1', 'IL2', 'IL3', 'IL6R', 'IL6ST'],
    'IL1B': ['VEGFA', 'PTGS2', 'IL6', 'IL5', 'IL2', 'IL4', 'IFNG', 'MYD88', 'IL1A', 'IL1R1', 'IL1R2'],
    'VEGFA': ['IL1B', 'IL6', 'CDKN2A'],
    'PTGS2': ['IL1B', 'IL6', 'CXCL8', 'TP53', 'ALOX5', 'ALOX15'],
    'CCL2': ['CXCL8', 'CCL5', 'CCL8', 'CXCL13', 'CCR1', 'CCL4L2', 'CCL15', 'CCR2', 'CXCL9', 'CCR5'],
    'CXCL8': ['IL4', 'CCL8', 'CXCR3', 'CXCR1', 'CCR1', 'CCR2', 'PF4', 'CCL2'],
    'RELA': ['CREBBP', 'EP300', 'NFKB1', 'NFKB2', 'NFKB1E', 'NFKB1B', 'RELB'],
    'TNF': ['CHUK', 'TRADD', 'NFKB1', 'TRAF2', 'TAB2','TNFRSF1B', 'TNFRSF1A', 'MAP3K7'],
    'JUN': ['ATF2', 'FOS', 'EP300', 'JUND', 'JUNB'],
    'MAPK1': ['JUN', 'STAT3', 'TP53', 'MAP2K2', 'MAP2K1'],
    'MAPK14': ['JUN', 'ATF2', 'TP53', 'TAB1'],
    'JAK1': ['STAT3', 'STAT1', 'IL2', 'IL4'],
    'JAK2': ['STAT3', 'STAT1', 'IL3'],
    'IL3': ['IL5', 'IL4', 'IL2'],
    'STAT3': ['EP300'],
    'EP300': ['TP53'],
    'TP53': ['CREBBP'],
    'MMP1': ['MMP3', 'MMP2', 'MMP9'],
    'NOS2': ['AKT1'],
    'RORC': ['IL17A', 'FOXP3', 'TBX21'],
    'TLR4': ['TLR2', 'TIRAP', 'LY96', 'TRAF6', 'IRAK4', 'TICAM1', 'TICAM2']
}

# Add secondary edges
for node, partners in interactions.items():
    for partner in partners:
        if not G.has_edge(node, partner):
            G.add_edge(node, partner, connection='indirect')

# Positioning
plt.figure(figsize=(18, 15))
pos = nx.spring_layout(G, k=0.45, seed=42)

# Node styling
node_colors = []
node_sizes = []
for node in G.nodes():
    if node in compounds:
        node_colors.append('gold')
        node_sizes.append(2000)
    elif len(target_map[node]) == 3:
        node_colors.append('red')
        node_sizes.append(1600)
    elif 'Humic Acid' in target_map[node]:
        node_colors.append('lightskyblue')
        node_sizes.append(1300)
    elif 'Fulvic Acid' in target_map[node]:
        node_colors.append('lightpink')
        node_sizes.append(1300)
    elif 'Withaferin A' in target_map[node]:
        node_colors.append('khaki')
        node_sizes.append(1300)
    else:
        node_colors.append('lightgray')
        node_sizes.append(700)

# Edge categories
direct_edges = [(u, v) for u, v, d in G.edges(data=True) if d['connection'] == 'direct']
indirect_edges = [(u, v) for u, v, d in G.edges(data=True) if d['connection'] == 'indirect']

# Draw nodes
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, edgecolors='black', linewidths=1.5)

# Draw edges
nx.draw_networkx_edges(G, pos, edgelist=direct_edges, width=3, edge_color='crimson')
nx.draw_networkx_edges(G, pos, edgelist=indirect_edges, width=1, style='dashed', edge_color='gray')

# Draw labels
nx.draw_networkx_labels(G, pos, font_size=9, font_family='sans-serif')

# Title
plt.title("Target Interaction Network of Withaferin A, Fulvic Acid & Humic Acid (RA Related)", fontsize=20, fontweight='bold', pad=20)

# Axis and background
plt.axis('off')
plt.gca().set_facecolor('whitesmoke')

plt.show()

