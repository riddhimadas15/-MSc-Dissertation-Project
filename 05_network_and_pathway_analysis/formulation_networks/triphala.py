import networkx as nx
import matplotlib.pyplot as plt

# Initialize graph
G = nx.Graph()

# Define compounds and their primary targets
compounds = {
    "Ascorbic acid": ["IL6", "IL1B", "IL11", "CCL2", "MMP1", "MMP9", "TNFRSF11A", "PTGS2", "STAT3", "TYK2", "MAPK1", "MAPK14", "RELA", "TNF", "NFKB1", "BAX", "CASP3", "CAT", "CXCL8", "CXCL10", "HIF1A", "NFE2L2", "NOS2", "SOD1", "TGFB1", "VEGFA", "JUN", "RORC"],
    "Gallic acid": ["IL6", "IL1B", "MMP9", "PTGS2", "JUN", "MAPK1", "RELA", "TNF", "AKT1", "BAX", "CASP3", "CAT", "CXCL8", "HIF1A", "IFNG", "IKBKB", "NFE2L2", "VEGFA"],
    "Chebulagic acid": ["PTGS2", "MAPK14", "HIF1A", "IKBKB"],
    "Ellagic acid": ["IL1B", "CCL2", "MMP1", "MMP9", "JUN", "STAT3", "RELA", "TNF", "AKT1", "BAX", "CASP3", "CXCL8", "HIF1A", "NFE2L2", "NOS2", "TGFB1", "ALOX5", "TEK", "PTGS2", "NFKB1", "IKBKB"],
    "Guggulsterone-E": ["IL1B", "MMP9", "PTGS2", "RELA", "TNF", "NFKB1", "CASP3", "TGFB1", "VEGFA", "ALOX5", "CTSL", "CTSK", "NOS2"],
    "Guggulsterone-Z": ["ALOX5", "CTSK", "IKBKB", "NOS2"]
}

# Add compound-target edges
for compound, targets in compounds.items():
    for target in targets:
        G.add_edge(compound, target)

# Define secondary interactions (partial list to manage graph complexity)
interactions = [
    ("IL6", "IL1B"), ("IL6", "PTGS2"), ("IL6", "VEGFA"), ("IL6", "JAK2"), ("IL6", "JAK1"), ("IL6", "IL2"), ("IL6", "IL3"), ("IL6", "IL11"),
    ("IL1B", "VEGFA"), ("IL1B", "PTGS2"), ("IL1B", "IL5"), ("IL1B", "IL2"), ("IL1B", "IL4"), ("IL1B", "IFNG"), ("IL1B", "MYD88"),
    ("CCL2", "CXCL8"), ("CCL2", "CCL5"), ("CCL2", "CCL8"), ("CCL2", "CXCL13"), ("CCL2", "CCR1"), ("CCL2", "CCL4L2"), ("CCL2", "CCL15"), ("CCL2", "CCR2"), ("CCL2", "CXCL9"), ("CCL2", "CCR5"),
    ("TNFRSF11A", "TNFSF11"), ("TNFRSF11A", "TRAF6"), ("TNFRSF11A", "TNF"),
    ("TYK2", "JAK1"), ("TYK2", "STAT2"), ("TYK2", "SOCS1"),
    ("MAPK1", "STAT3"), ("MAPK1", "TP53"), ("MAPK1", "MAP2K2"), ("MAPK1", "MAP2K1"),
    ("MAPK14", "JUN"), ("MAPK14", "ATF2"), ("MAPK14", "TP53"), ("MAPK14", "TAB1"),
    ("MMP9", "TNF"), ("MMP9", "TLR4"), ("MMP9", "IL1B"), ("MMP9", "MMP3"), ("MMP9", "TGFB1"),
    ("SOD1", "CAT"),
    ("RELA", "CREBBP"), ("RELA", "EP300"), ("RELA", "NFKB1A"), ("RELA", "NFKB2"),
    ("TNF", "CHUK"), ("TNF", "TRADD"), ("TNF", "NFKB1"), ("TNF", "TRAF2"), ("TNF", "TAB2"),
    ("JUN", "ATF2"), ("JUN", "FOS"), ("JUN", "EP300"), ("JUN", "JUND"), ("JUN", "JUNB"),
    ("HIF1A", "STAT3"), ("HIF1A", "MAPK1"), ("HIF1A", "RELA"),
    ("RORC", "IL17A"), ("RORC", "FOXP3"),
    ("STAT3", "STAT1"), ("STAT3", "JAK1"), ("STAT3", "EP300"),
    ("NFE2L2", "JUN"), ("NFE2L2", "AKT1"),
    ("BAX", "CYCS"), ("BAX", "TP53"),
    ("PTGS2", "CXCL8"), ("PTGS2", "TP53"), ("PTGS2", "ALOX5"), ("PTGS2", "ALOX15"), ("PTGS2", "TNF"), ("PTGS2", "NOS2"),
    ("TGFB1", "IL10"), ("TGFB1", "IL1B"), ("TGFB1", "TNF"), ("TGFB1", "MMP9"),
    ("CASP3", "AKT1"), ("CASP3", "MAPK1"),
    ("IFNG", "IL10"), ("IFNG", "IL4"), ("IFNG", "TNF"), ("IFNG", "CD4"), ("IFNG", "JAK1"), ("IFNG", "JAK2")
]

# Add secondary interaction edges
for a, b in interactions:
    G.add_edge(a, b)

# Draw graph
plt.figure(figsize=(28, 22))
pos = nx.spring_layout(G, k=0.35)
nx.draw(G, pos, with_labels=True, node_size=900, node_color="lightblue", edge_color="gray", font_size=10)
plt.title("Triphala Target Interaction Network", fontsize=16)
plt.tight_layout()
plt.show()

