from textwrap import wrap
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from distinctipy import distinctipy

df = pd.read_csv("S_20_genes_for_figure.csv")
# Load gene-to-color mapping from a CSV file
gene_color_mapping = []
try:
    with open("gene_colors.csv", "r") as file:
        for line in file:
            gene_color_mapping.append(line.strip().split(","))
except Exception as e:
    raise RuntimeError(f"Error reading gene_colors.csv: {e}")

predefined_colors = [
(1,0.9333333333333333,0.5098039215686274),
(231/255, 76/255, 60/255,0.50),
(0.6352941176470588,0.3215686274509804,0.5019607843137255,0.50),
(0,0.4627450980392157,0.3058823529411765,0.50),
(0.09803921568627451,0.5490196078431373,0.7490196078431373,0.50),
(230/255, 126/255, 34/255,0.50),
(0.82421875,0.82421875,0.82421875),
]


# Assign predefined colors to gene groups
color_dict = {}
for group, color in zip(gene_color_mapping, predefined_colors):
    for gene in group:
        gene = gene.strip()
        color_dict[gene] = color

# Ensure the dataframe is sorted by the order of appearance in the file
df["Operon"] = pd.Categorical(df["Operon"], categories=df["Operon"].unique(), ordered=True)
df = df.sort_values("Operon")

# Assign colors from the mapping file to each gene name
def get_gene_color(gene_name):
    return color_dict.get(gene_name, "#FFFFFF")  # Default to grey if gene not in mapping

# Group data by operon
operons = []
operon_names = []
for operon_name, group in df.groupby("Operon", sort=False):  # Disable sorting in groupby
    operon_names.append(operon_name.split('_')[0])
    genes = []
    for _, row in group.iterrows():
        genes.append({
            "gene": row["Gene"],
            "start": row["Start"],
            "end": row["End"],
            "strand": row["Strand"],
            "color": get_gene_color(row["Gene"])
        })
    operons.append(genes)

# Create figure and axis
fig, ax = plt.subplots(figsize=(len(operons)*2, len(operons)))

# Define a fixed triangle length
triangle_length = 300

# Draw each operon
for idx, genes in enumerate(operons):
    y = len(operons) - idx  # Vertical position for each operon
    for gene in genes:
        height = 0.5  # Height of the rectangles

        # Adjust triangle length to ensure it is contained within the gene boundaries
        effective_triangle_length = min(triangle_length, (gene['end'] - gene['start']) / 2)

        # Create a tapered rectangle for the gene
        if gene['strand'] == "+":
            polygon = patches.Polygon(
                [
                    [gene['start'], y - height / 2],  # Bottom left
                    [gene['end'] - effective_triangle_length, y - height / 2],  # Bottom right (taper start)
                    [gene['end'], y],  # Tip of the taper
                    [gene['end'] - effective_triangle_length, y + height / 2],  # Top right (taper start)
                    [gene['start'], y + height / 2]  # Top left
                ],
                closed=True, edgecolor='grey', facecolor=gene['color']
            )
        else:
            polygon = patches.Polygon(
                [
                    [gene['end'], y - height / 2],  # Bottom right
                    [gene['start'] + effective_triangle_length, y - height / 2],  # Bottom left (taper start)
                    [gene['start'], y],  # Tip of the taper
                    [gene['start'] + effective_triangle_length, y + height / 2],  # Top left (taper start)
                    [gene['end'], y + height / 2]  # Top right
                ],
                closed=True, edgecolor='grey', facecolor=gene['color']
            )
        ax.add_patch(polygon)

        # Add gene label inside the rectangle
        if gene['gene'] == 'hpnS' or gene['gene'] == 'hpnT' or gene['gene'] == 'hpnL' or gene['gene'] == 'hpnM' or gene['gene'] == 'hpnM*':
            ax.text((gene['start'] + gene['end']) / 2, y, gene['gene'],ha='center', va='center', fontsize=18, color='black', fontweight='bold', style='italic')
        elif gene['gene'] != '-':
            ax.text((gene['start'] + gene['end']) / 2, y, gene['gene'],ha='center', va='center', fontsize=18, color='black', fontweight='bold')



# Adjust plot limits and labels
ax.set_xlim(df['Start'].min() - 500, df['End'].max() + 500)  # Adjust to accommodate negative positions
ax.set_ylim(0, len(operons) + 1)
ax.set_xlabel("Genomic Position", fontsize=24)
ax.set_yticks(range(1, len(operons) + 1))
ax.set_yticklabels(operon_names[::-1], fontsize=24, fontweight='bold')
ax.tick_params(axis='x', labelsize=24)
# Create legend at the bottom
fig.subplots_adjust(bottom=0.2)
legend_ax = fig.add_axes([0.125, 0.05, 0.8, 0.1])  # Custom axis for the legend
legend_ax.set_xlim(0, len(predefined_colors))
legend_ax.set_ylim(0, 1)
legend_ax.axis("off")

legend_labels = [
#"SHC, SHC*, hpnA, hpnA*, hpnB, hpnB*, hpnBC, hpnBD, hpnC, hpnD, hpnCD, hpnE, hpnE*, hpnG, hpnH, hpnH*, hpnI, hpnJ, hpnN, hpnN*, hpnO, hpnP*",
"Hopanoid",
#"IPPI, IPPI*, ispH, ispG, DXPS, DXPS*, yceI*",
"Isoprenoid",
"hpnM",
"hpnS",
"hpnT",
"hpnL",
"Other"
#"99, 23, 172, 134, 144, 219, 176, 69, 22, UNK3*, 25, 268, 15, 94, 110, 137, arsC, 88, glcB, 123, 81, 33, 82, 170, 161, 173, 177, UNK3, 14, 10, 147, 29, 248, 21, 41, 24, 149, 103, 130, 124, 119, 139, 159, 146, rpoD, 141, 60, 84, 213, 255, 114, 194, 42, 102, 13, 19, 58, 61, 192, 1, 72, 179, 135"
]
for i, color in enumerate(predefined_colors):
	rect = patches.Rectangle((i*0.6+1.5, 0.5), 0.220/2,0.25, facecolor=color, edgecolor='black')
	rect.set_clip_on(False)
	legend_ax.add_patch(rect)
	legend_ax.text(i*0.6+1.5 + 0.15, 0.625,  '\n'.join(wrap(legend_labels[i], 25)), va='center', fontsize=24, fontweight='bold')

# Save the plot as a PDF
plt.savefig("gene_neighborhood_figure.pdf")
plt.close(fig)
