#!/usr/bin/env python3
"""
Create PowerPoint presentation explaining the KAT8 Sex x Genotype
interaction analysis (part7.6) for presentation to PI / boss.

Run:  python create_presentation.py
Output:  Results/KAT8_sex_interaction_presentation.pptx
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
DARK_BLUE  = RGBColor(0x1B, 0x3A, 0x5C)
MED_BLUE   = RGBColor(0x2C, 0x6F, 0xAD)
LIGHT_BLUE = RGBColor(0xD6, 0xE8, 0xF7)
WHITE      = RGBColor(0xFF, 0xFF, 0xFF)
BLACK      = RGBColor(0x00, 0x00, 0x00)
GREY       = RGBColor(0x66, 0x66, 0x66)
LIGHT_GREY = RGBColor(0xDD, 0xDD, 0xDD)
RED_ACC    = RGBColor(0xC0, 0x39, 0x2B)
GREEN_ACC  = RGBColor(0x27, 0xAE, 0x60)

SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)


def set_slide_bg(slide, color):
    """Set solid background colour for a slide."""
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_textbox(slide, left, top, width, height, text,
                font_size=16, bold=False, color=BLACK, alignment=PP_ALIGN.LEFT,
                font_name="Calibri"):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    return tf


def add_bullet_list(slide, left, top, width, height, items,
                    font_size=14, color=BLACK, font_name="Calibri",
                    spacing_after=Pt(6)):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = item
        p.font.size = Pt(font_size)
        p.font.color.rgb = color
        p.font.name = font_name
        p.space_after = spacing_after
        p.level = 0
    return tf


def add_plot_placeholder(slide, left, top, width, height, filename):
    """Grey dashed box with instructions to insert the named plot."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = RGBColor(0xF2, 0xF2, 0xF2)
    shape.line.color.rgb = RGBColor(0x99, 0x99, 0x99)
    shape.line.dash_style = 2  # dash
    shape.line.width = Pt(1.5)
    tf = shape.text_frame
    tf.word_wrap = True
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    p = tf.paragraphs[0]
    p.text = "INSERT PLOT"
    p.font.size = Pt(14)
    p.font.bold = True
    p.font.color.rgb = GREY
    p2 = tf.add_paragraph()
    p2.text = filename
    p2.font.size = Pt(11)
    p2.font.color.rgb = GREY
    p2.font.italic = True
    p2.alignment = PP_ALIGN.CENTER


def add_section_header(slide, text, top=Inches(0.3)):
    add_textbox(slide, Inches(0.6), top, Inches(12), Inches(0.6),
                text, font_size=28, bold=True, color=DARK_BLUE)


def add_subtitle_text(slide, text, top=Inches(0.9)):
    add_textbox(slide, Inches(0.6), top, Inches(12), Inches(0.4),
                text, font_size=14, color=GREY)


# ---------------------------------------------------------------------------
# Build Presentation
# ---------------------------------------------------------------------------
prs = Presentation()
prs.slide_width = SLIDE_W
prs.slide_height = SLIDE_H
blank_layout = prs.slide_layouts[6]  # blank

# ====== SLIDE 1: Title ======
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, DARK_BLUE)
add_textbox(slide, Inches(1), Inches(2), Inches(11), Inches(1.5),
            "KAT8-KD RNA-seq:\nSex x Genotype Interaction Analysis",
            font_size=36, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(1), Inches(4), Inches(11), Inches(0.8),
            "Part 7.6  |  Formal test of sex-dependent transcriptional responses",
            font_size=18, color=RGBColor(0xAA, 0xCC, 0xEE), alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(1), Inches(5.5), Inches(11), Inches(0.5),
            "Adipose Depot Analysis: iWAT & gWAT",
            font_size=16, color=RGBColor(0x88, 0xAA, 0xCC), alignment=PP_ALIGN.CENTER)

# ====== SLIDE 2: Study Design ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Study Design")
add_bullet_list(slide, Inches(0.8), Inches(1.2), Inches(6), Inches(5), [
    "2 x 2 factorial design: Sex (F, M) x Genotype (CTL, KAT8-KD)",
    "Two adipose depots analyzed separately: iWAT and gWAT",
    "40 samples total, 2 outliers removed (JS_08, JS_28) -> 38 samples",
    "Per depot: 19 samples (after outlier removal)",
    "",
    "Group sizes per depot:",
    "   CTL Female:  5     |  KAT8-KD Female:  4",
    "   CTL Male:    5     |  KAT8-KD Male:    5",
], font_size=16)

# ====== SLIDE 3: The Question ======
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, RGBColor(0xF0, 0xF4, 0xF8))
add_textbox(slide, Inches(1), Inches(1.5), Inches(11), Inches(1),
            "The Central Question",
            font_size=32, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(1.5), Inches(3), Inches(10), Inches(1.2),
            "Does KAT8 knockdown affect gene expression\n"
            "differently in males versus females?",
            font_size=26, bold=True, color=MED_BLUE, alignment=PP_ALIGN.CENTER)
add_textbox(slide, Inches(1.5), Inches(4.8), Inches(10), Inches(1),
            "If the answer is NO (for most genes), we can justify pooling sexes\n"
            "in the primary DEG analysis, increasing statistical power.",
            font_size=16, color=GREY, alignment=PP_ALIGN.CENTER)

# ====== SLIDE 4: Method - LRT ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Method: Likelihood Ratio Test (LRT)")
add_bullet_list(slide, Inches(0.8), Inches(1.3), Inches(11.5), Inches(5.5), [
    "Tool: DESeq2 (Bioconductor), run separately per depot",
    "",
    "Full model:      ~ Sex + Genotype + Sex:Genotype",
    "Reduced model:  ~ Sex + Genotype",
    "",
    "The LRT asks: does adding the interaction term (Sex:Genotype) "
    "significantly improve the model fit for each gene?",
    "",
    "If YES -> that gene's KD response differs between M and F",
    "If NO  -> the KD effect is the same in both sexes",
    "",
    "We also run 4 complementary analyses to strengthen the conclusion:",
    "  (1) Pi0 estimation  (2) Variance attribution  "
    "(3) Sex-stratified concordance  (4) PCA"
], font_size=15)

# ====== SLIDE 5: iWAT LRT Results ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "iWAT: LRT Interaction Results")
add_bullet_list(slide, Inches(0.6), Inches(1.2), Inches(5.5), Inches(5), [
    "16,466 genes tested (after prefiltering)",
    "630 significant interactions (3.83%, FDR < 0.05)",
    "~96% of genes show NO sex-dependent KD response",
    "",
    "How to read the histogram:",
    "  - Flat distribution = genes with no interaction (null)",
    "  - Spike near p=0 = genes with real interaction signal",
    "  - Red dashed line = expected count if ALL genes were null",
    "  - The histogram is mostly flat -> most genes are null",
], font_size=14)
add_plot_placeholder(slide, Inches(6.5), Inches(1.2), Inches(6), Inches(5),
                     "pvalue_hist_interaction_iWAT_<run_tag>.png")

# ====== SLIDE 6: gWAT LRT Results ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "gWAT: LRT Interaction Results")
add_bullet_list(slide, Inches(0.6), Inches(1.2), Inches(5.5), Inches(5), [
    "16,827 genes tested (after prefiltering)",
    "374 significant interactions (2.22%, FDR < 0.05)",
    "~98% of genes show NO sex-dependent KD response",
    "",
    "Even fewer interaction genes than iWAT,",
    "supporting sex-independent KD effects in gWAT as well.",
], font_size=14)
add_plot_placeholder(slide, Inches(6.5), Inches(1.2), Inches(6), Inches(5),
                     "pvalue_hist_interaction_gWAT_<run_tag>.png")

# ====== SLIDE 7: Pi0 Estimation ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Pi0: Estimating the True Null Proportion")
add_bullet_list(slide, Inches(0.6), Inches(1.3), Inches(12), Inches(5.5), [
    "What is pi0?",
    "  Pi0 = the fraction of genes where the null hypothesis (no interaction) is truly correct",
    "  Not just 'failed to reject' -- this estimates the actual proportion of true nulls",
    "  Method: Storey & Tibshirani (qvalue package), uses shape of p-value distribution",
    "",
    "Results:",
    "  iWAT:  pi0 = 0.87  ->  87% of genes genuinely have NO sex-dependent KD response",
    "  gWAT:  pi0 = 0.85  ->  85% of genes genuinely have NO sex-dependent KD response",
    "",
    "Interpretation:",
    "  The vast majority of genes respond identically to KAT8-KD in both sexes",
    "  The remaining ~13-15% includes both true interactions AND residual noise",
    "  This strongly supports pooling sexes for the primary analysis",
], font_size=14)

# ====== SLIDE 8: Variance Attribution ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Variance Attribution: What Drives Gene Expression?")
add_subtitle_text(slide,
    "ANOVA decomposition on VST-transformed counts: "
    "lm(expression ~ Sex + Genotype + Sex:Genotype)")
add_bullet_list(slide, Inches(0.6), Inches(1.6), Inches(5.5), Inches(2.5), [
    "For each gene, how much variance is explained by each factor?",
    "",
    "         Term          iWAT (median)   gWAT (median)",
    "  Sex                       6.5%            11.5%",
    "  Genotype (KD effect)    33.2%            16.7%",
    "  Sex x Genotype            2.1%              2.4%",
    "  Residual                  46.3%            50.1%",
], font_size=14)
add_bullet_list(slide, Inches(0.6), Inches(4.5), Inches(5.5), Inches(2), [
    "Key takeaway: The interaction term explains only ~2% of variance",
    "  vs. 17-33% for the Genotype main effect (8-16x more)",
    "The KD signal dominates; sex-dependent effects are minimal",
], font_size=14, color=DARK_BLUE)
add_plot_placeholder(slide, Inches(6.5), Inches(1.5), Inches(6), Inches(5),
                     "variance_explained_<depot>_<run_tag>.png")

# ====== SLIDE 9: Concordance iWAT ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Sex-Stratified Concordance: iWAT")
add_subtitle_text(slide,
    "Run DESeq2 separately for F-only and M-only, then correlate KD log2FCs")
add_bullet_list(slide, Inches(0.6), Inches(1.6), Inches(5.5), Inches(4.5), [
    "Pearson r = 0.727 (moderate-to-strong concordance)",
    "",
    "How to read this scatter plot:",
    "  X-axis = Female KD log2FC   Y-axis = Male KD log2FC",
    "  Diagonal line = perfect concordance (M = F)",
    "",
    "  Purple dots: significant in BOTH sexes (strongest evidence)",
    "  Red dots: significant in F only",
    "  Blue dots: significant in M only",
    "  Grey dots: not significant in either",
    "",
    "Points along the diagonal = same response in both sexes",
    "r = 0.73 indicates good concordance -- both sexes largely agree",
], font_size=13)
add_plot_placeholder(slide, Inches(6.5), Inches(1.5), Inches(6), Inches(5),
                     "sex_concordance_scatter_iWAT_<run_tag>.png")

# ====== SLIDE 10: Concordance gWAT ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Sex-Stratified Concordance: gWAT")
add_subtitle_text(slide,
    "Run DESeq2 separately for F-only and M-only, then correlate KD log2FCs")
add_bullet_list(slide, Inches(0.6), Inches(1.6), Inches(5.5), Inches(4.5), [
    "Pearson r = 0.605 (moderate concordance)",
    "",
    "Slightly lower than iWAT -- visceral fat (gWAT) shows",
    "more baseline sex dimorphism than subcutaneous fat (iWAT),",
    "which is consistent with known adipose biology.",
    "",
    "However, the direction of KD effects is still positively",
    "correlated -- genes up in F are generally also up in M.",
    "",
    "Note: With only n=4-5 per sex-genotype group, per-sex",
    "analyses have limited power, so the true concordance",
    "is likely higher than observed.",
], font_size=13)
add_plot_placeholder(slide, Inches(6.5), Inches(1.5), Inches(6), Inches(5),
                     "sex_concordance_scatter_gWAT_<run_tag>.png")

# ====== SLIDE 11: PCA ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "PCA by Sex and Genotype")
add_bullet_list(slide, Inches(0.6), Inches(1.3), Inches(5.5), Inches(5), [
    "Principal Component Analysis of top 500 variable genes",
    "",
    "How to read this plot:",
    "  Color = Genotype (CTL vs KAT8-KD)",
    "  Shape = Sex (circle = F, triangle = M)",
    "",
    "What to look for:",
    "  If M and F of the same genotype cluster together",
    "  -> sexes are transcriptionally similar",
    "",
    "  If M and F of the same genotype separate",
    "  -> sex differences dominate over KD effects",
    "",
    "Expected: genotype drives separation more than sex",
], font_size=14)
add_plot_placeholder(slide, Inches(6.5), Inches(1.2), Inches(6), Inches(5),
                     "PCA_by_sex_<depot>_<run_tag>.png")

# ====== SLIDE 12: Interaction vs Genotype Scatter ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Interaction Effect vs. Genotype Main Effect")
add_bullet_list(slide, Inches(0.6), Inches(1.3), Inches(5.5), Inches(5), [
    "X-axis = Genotype main effect log2FC (KD vs CTL)",
    "Y-axis = Interaction log2FC (sex-dependent component)",
    "",
    "Red dots: significant interaction genes",
    "Blue dots: significant genotype DEGs (no interaction)",
    "Grey dots: neither significant",
    "",
    "What to look for:",
    "  Most dots clustered near y = 0 means the interaction",
    "  effect is small relative to the genotype effect.",
    "",
    "  The horizontal spread (x-axis) >> vertical spread (y-axis)",
    "  confirms that KD drives expression, not sex x KD.",
], font_size=14)
add_plot_placeholder(slide, Inches(6.5), Inches(1.2), Inches(6), Inches(5),
                     "interaction_vs_genotype_<depot>_<run_tag>.png")

# ====== SLIDE 13: Pathway Analysis (if applicable) ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Interaction Gene Pathways (ORA)")
add_subtitle_text(slide,
    "Over-representation analysis on genes with significant Sex x Genotype interaction")
add_bullet_list(slide, Inches(0.6), Inches(1.6), Inches(5.5), Inches(5), [
    "For each depot, the significant interaction genes were tested",
    "for enrichment in GO:BP and KEGG pathways.",
    "",
    "This tells us WHICH biological processes have sex-dependent",
    "KAT8-KD responses (even though most genes do not).",
    "",
    "Method: clusterProfiler enrichGO / enrichKEGG",
    "  Universe = all tested genes in that depot",
    "  FDR < 0.05 significance threshold",
    "",
    "Insert the dotplots from the output directory:",
    "  interaction_ORA_GOBP_dotplot_<depot>_<run_tag>.png",
    "  interaction_ORA_KEGG_dotplot_<depot>_<run_tag>.png",
], font_size=13)
add_plot_placeholder(slide, Inches(6.5), Inches(1.5), Inches(6), Inches(5),
                     "interaction_ORA_GOBP_dotplot_<depot>_<run_tag>.png")

# ====== SLIDE 14: Summary Table ======
slide = prs.slides.add_slide(blank_layout)
add_section_header(slide, "Summary: All Evidence at a Glance")

# Build a table
rows = 8
cols = 3
table_shape = slide.shapes.add_table(rows, cols,
    Inches(1), Inches(1.3), Inches(11), Inches(4.5))
tbl = table_shape.table

# Header row
headers = ["Metric", "iWAT", "gWAT"]
for i, h in enumerate(headers):
    cell = tbl.cell(0, i)
    cell.text = h
    p = cell.text_frame.paragraphs[0]
    p.font.bold = True
    p.font.size = Pt(14)
    p.font.color.rgb = WHITE
    p.font.name = "Calibri"
    cell.fill.solid()
    cell.fill.fore_color.rgb = DARK_BLUE

# Data rows
data = [
    ["Genes tested",                 "16,466",        "16,827"],
    ["Significant interactions",     "630 (3.83%)",   "374 (2.22%)"],
    ["Pi0 (true null proportion)",   "0.87",          "0.85"],
    ["Genotype DEGs",                "2,691",         "1,054"],
    ["Interaction/Genotype ratio",   "0.234",         "0.355"],
    ["Median variance: interaction", "2.1%",          "2.4%"],
    ["Sex-stratified concordance r", "0.727",         "0.605"],
]

for r_idx, row_data in enumerate(data):
    for c_idx, val in enumerate(row_data):
        cell = tbl.cell(r_idx + 1, c_idx)
        cell.text = val
        p = cell.text_frame.paragraphs[0]
        p.font.size = Pt(13)
        p.font.name = "Calibri"
        if c_idx == 0:
            p.font.bold = True
        if r_idx % 2 == 0:
            cell.fill.solid()
            cell.fill.fore_color.rgb = LIGHT_BLUE

add_textbox(slide, Inches(1), Inches(6), Inches(11), Inches(0.8),
            "All metrics converge: the KAT8-KD transcriptional response is "
            "overwhelmingly sex-independent in both depots.",
            font_size=16, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)

# ====== SLIDE 15: Conclusion ======
slide = prs.slides.add_slide(blank_layout)
set_slide_bg(slide, DARK_BLUE)
add_textbox(slide, Inches(1), Inches(1), Inches(11), Inches(1),
            "Conclusion",
            font_size=34, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

conclusion_items = [
    "Only 2-4% of genes show significant Sex x Genotype interaction",
    "Pi0 estimates confirm 85-87% of genes are true nulls (no interaction)",
    "The interaction term explains only ~2% of per-gene variance vs. 17-33% for Genotype",
    "Sex-stratified KD effects are positively correlated (r = 0.61-0.73)",
    "",
    "JUSTIFICATION: Male and female adipose tissue respond similarly",
    "to KAT8 knockdown. Pooling sexes for the primary DEG analysis is appropriate.",
    "",
    "CAVEAT: With n = 4-5 per group, subtle sex-dependent effects",
    "may be undetected. The 630 (iWAT) / 374 (gWAT) interaction genes",
    "warrant further exploration for sex-specific biology."
]

tf = add_bullet_list(slide, Inches(1.5), Inches(2.2), Inches(10), Inches(4.5),
                     conclusion_items, font_size=16, color=WHITE)

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Results")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "KAT8_sex_interaction_presentation.pptx")
prs.save(out_path)
print(f"[OK] Presentation saved to: {out_path}")
print(f"     {len(prs.slides)} slides created")
print()
print("Next steps:")
print("  1. Open the PPTX on your Windows machine")
print("  2. Replace the grey placeholder boxes with your actual plots")
print("     (Right-click placeholder -> Change Shape -> Insert Picture)")
print("  3. The plot filenames are listed in each placeholder")
