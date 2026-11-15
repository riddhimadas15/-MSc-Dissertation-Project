import matplotlib.pyplot as plt
import networkx as nx
import textwrap
from collections import defaultdict

# Initialize graph
G = nx.Graph()

# Formulation nodes
formulations = {
    "Ashwashila": "gold",
    "Articulin F": "gold",
    "Triphala Guggul": "gold"
}

# Pathway mappings
formulation_pathways = {
    "Ashwashila": [
        "Cytokine–cytokine receptor interaction",
        "NF-κB signaling",
        "MAPK signaling",
        "JAK–STAT signaling",
        "Chemokine signaling",
        "Toll-like receptor signaling",
        "Apoptosis",
        "VEGF signaling",
        "Arachidonic acid metabolism",
        "Nrf2 antioxidant response",
        "Th17 differentiation",
        "Nitric oxide synthase",
        "Extracellular matrix remodeling"
    ],
    "Articulin F": [
        "Cytokine–cytokine receptor interaction",
        "NF-κB signaling",
        "MAPK signaling",
        "JAK–STAT signaling",
        "Chemokine signaling",
        "Toll-like receptor signaling",
        "Apoptosis",
        "VEGF signaling",
        "Arachidonic acid metabolism",
        "Nrf2 antioxidant response",
        "PI3K–Akt signaling",
        "TGF-β signaling",
        "Extracellular matrix remodeling",
        "Nitric oxide synthase"
    ],
    "Triphala Guggul": [
        "Cytokine–cytokine receptor interaction",
        "NF-κB signaling",
        "MAPK signaling",
        "JAK–STAT signaling",
        "Chemokine signaling",
        "Toll-like receptor signaling",
        "Apoptosis",
        "VEGF signaling",
        "Arachidonic acid metabolism",
        "Nrf2 antioxidant response",
        "PI3K–Akt signaling",
        "TGF-β signaling",
        "Extracellular matrix remodeling",
        "Nitric oxide synthase",
        "Th17 differentiation",
        "Angiogenesis"
    ]
}

# Count how many formulations target each pathway
pathway_counts = defaultdict(int)
for paths in formulation_pathways.values():
    for path in paths:
        pathway_counts[path] += 1

# Add nodes with color based on type
for form, color in formulations.items():
    G.add_node(form, type='formulation', color=color)

# Define colors for pathway nodes based on how many formulations target them
def get_pathway_color(count):
    if count == 1:
        return "lightgray"
    elif count == 2:
        return "lightgreen"
    else:
        return "tomato"

# Add pathway nodes with color based on count
for path, count in pathway_counts.items():
    G.add_node(path, type='pathway', color=get_pathway_color(count))

# Add edges
for form in formulations:
    for pathway in formulation_pathways[form]:
        G.add_edge(form, pathway)

# Layout
pos = nx.spring_layout(G, k=0.5, iterations=200)

# Separate nodes
form_nodes = [n for n, attr in G.nodes(data=True) if attr['type'] == 'formulation']
path_nodes = [n for n, attr in G.nodes(data=True) if attr['type'] == 'pathway']

# Draw nodes
nx.draw_networkx_nodes(G, pos,
                       nodelist=form_nodes,
                       node_color=[G.nodes[n]['color'] for n in form_nodes],
                       node_shape='s',
                       node_size=1200)

nx.draw_networkx_nodes(G, pos,
                       nodelist=path_nodes,
                       node_color=[G.nodes[n]['color'] for n in path_nodes],
                       node_shape='s',
                       node_size=1000)

# Wrap labels
wrapped_labels = {n: textwrap.fill(n, width=18) for n in G.nodes()}

# Draw labels
nx.draw_networkx_labels(G, pos, labels=wrapped_labels,
                        font_size=8,
                        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", lw=0.5))

# Draw edges (same color for all now, since node color shows targeting)
nx.draw_networkx_edges(G, pos, edge_color='gray')

# Final plot
plt.title("Formulations → Targeted Pathways Network", fontsize=16)
plt.axis('off')
plt.tight_layout()
plt.show()

