#!/usr/bin/env python3
"""
Batch script: For each subdirectory, combine Fdm + Dsuite PDFs into a single editable vector PDF.
- Tree is drawn as vector graphics
- Original PDFs are embedded as Form XObjects (no rasterization, no compression)
- All text labels are editable
"""

import fitz  # PyMuPDF
import os
import re
import time
import sys

base = r"E:\Project\Rosacea\Rosaceae_cytonuclear\03-final_result_to_publish\12-Get_introgress_gene\orthofinder\KEGG_GO_rich"

# Subdirectories to skip
SKIP_DIRS = {'.claude', 'enrichment_results', 'org.DR.eg.db', 'org.Cf.eg.db'}

def log(msg):
    print(f"  [{time.strftime('%H:%M:%S')}] {msg}")
    sys.stdout.flush()

def get_outgroup(subdir_path):
    """Read the outgroup from mapping_sp.txt or map_sp.txt."""
    for fname in ('mapping_sp.txt', 'map_sp.txt'):
        fpath = os.path.join(subdir_path, fname)
        if os.path.exists(fpath):
            with open(fpath) as f:
                content = f.read().strip()
            # Format: "1\tSpecies_name" or just "Species_name"
            parts = content.split('\t')
            if len(parts) >= 2:
                return parts[1].strip()
            return parts[0].strip()
    return None

def parse_samples(dirname):
    """
    Parse P1, P2, P3 from directory name.
    Format: SRR{id}_SRR{id}_{rest}  or  SRR{id}_SRR{id}_SRR{id}
    """
    m = re.match(r'^(SRR\d+)_(SRR\d+)_(.+)$', dirname)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return None, None, None

def draw_tree(page, tree_x0, tree_y0, tree_w, tree_h, p1_name, p2_name, p3_name, outgroup_name):
    """Draw phylogenetic tree as vector graphics on the PDF page."""
    tm_left = 15
    tm_right = 85
    tm_top = 50
    tm_bottom = 50

    draw_w = tree_w - tm_left - tm_right
    draw_h = tree_h - tm_top - tm_bottom

    # Tip y positions
    y_p3 = tree_y0 + tm_bottom + draw_h * 0.10
    y_p2 = tree_y0 + tm_bottom + draw_h * 0.40
    y_p1 = tree_y0 + tm_bottom + draw_h * 0.70
    y_og = tree_y0 + tm_bottom + draw_h * 0.90

    # Node x positions
    x_tip = tree_x0 + tm_left + draw_w * 0.90
    x_node_p2p3 = tree_x0 + tm_left + draw_w * 0.52
    x_node_p1p2p3 = tree_x0 + tm_left + draw_w * 0.32
    x_root = tree_x0 + tm_left + draw_w * 0.12

    c_og = (0.2, 0.2, 0.2)
    c_p1 = (0.91, 0.30, 0.24)
    c_p2 = (0.20, 0.60, 0.86)
    c_p3 = (0.18, 0.80, 0.44)
    c_line = (0.2, 0.2, 0.2)
    c_arrow = (0.91, 0.30, 0.24)
    lw = 1.5

    def dl(x1, y1, x2, y2, color=c_line, width=lw):
        page.draw_line(fitz.Point(x1, y1), fitz.Point(x2, y2), color=color, width=width)

    # Horizontal branches
    dl(x_root, y_og, x_tip, y_og, c_og)
    dl(x_node_p1p2p3, y_p1, x_tip, y_p1, c_p1)
    dl(x_node_p2p3, y_p2, x_tip, y_p2, c_p2)
    dl(x_node_p2p3, y_p3, x_tip, y_p3, c_p3)

    # Vertical lines
    dl(x_root, y_og, x_root, y_p1)
    dl(x_node_p1p2p3, y_p1, x_node_p1p2p3, y_p3)
    dl(x_node_p2p3, y_p2, x_node_p2p3, y_p3)

    # Horizontal connectors
    y_mid = (y_p1 + y_p3) / 2
    dl(x_root, y_mid, x_node_p1p2p3, y_mid)
    dl(x_node_p1p2p3, y_p2, x_node_p2p3, y_p2)

    # Introgression arrow
    arrow_x = x_node_p2p3 + draw_w * 0.06
    arrow_yt = y_p2 + 8
    arrow_yb = y_p3 - 8
    dl(arrow_x, arrow_yb, arrow_x, arrow_yt, c_arrow, 1.5)
    asz = 6
    page.draw_polyline([
        fitz.Point(arrow_x - asz, arrow_yt - asz),
        fitz.Point(arrow_x, arrow_yt),
        fitz.Point(arrow_x + asz, arrow_yt - asz),
    ], color=c_arrow, width=1.5)

    page.insert_text(
        fitz.Point(arrow_x + 8, (arrow_yt + arrow_yb) / 2),
        "Introgression", fontname="helv", fontsize=6, color=c_arrow, rotate=90
    )

    # Tip labels
    lox = 8
    fs_tip = 8
    fs_sm = 7

    # Outgroup
    og_lines = _wrap_outgroup(outgroup_name)
    for i, line in enumerate(og_lines):
        page.insert_text(
            fitz.Point(x_tip + lox, y_og - 4 + i * 10),
            line, fontname="helv", fontsize=fs_sm, color=c_og
        )
    page.insert_text(
        fitz.Point(x_tip + lox, y_og - 4 + len(og_lines) * 10),
        "(Outgroup)", fontname="helv", fontsize=fs_sm - 1, color=c_og
    )

    # P1
    page.insert_text(fitz.Point(x_tip + lox, y_p1 - 2), p1_name, fontname="helv", fontsize=fs_tip, color=c_p1)
    page.insert_text(fitz.Point(x_tip + lox, y_p1 + 9), "(P1)", fontname="helv", fontsize=fs_sm, color=c_p1)

    # P2
    page.insert_text(fitz.Point(x_tip + lox, y_p2 - 2), p2_name, fontname="helv", fontsize=fs_tip, color=c_p2)
    page.insert_text(fitz.Point(x_tip + lox, y_p2 + 9), "(P2)", fontname="helv", fontsize=fs_sm, color=c_p2)

    # P3
    page.insert_text(fitz.Point(x_tip + lox, y_p3 - 2), p3_name, fontname="helv", fontsize=fs_tip, color=c_p3)
    page.insert_text(fitz.Point(x_tip + lox, y_p3 + 9), "(P3)", fontname="helv", fontsize=fs_sm, color=c_p3)


def _wrap_outgroup(name):
    """Wrap long outgroup names into multiple lines."""
    if len(name) <= 20:
        return [name]
    # Split by underscore
    parts = name.split('_')
    lines = []
    current = ''
    for p in parts:
        if current and len(current + '_' + p) > 22:
            lines.append(current)
            current = p
        else:
            current = p if not current else current + '_' + p
    if current:
        lines.append(current)
    return lines if lines else [name]


def process_subdir(subdir_name):
    """Process a single subdirectory."""
    subdir_path = os.path.join(base, subdir_name)
    
    # Find required files
    fdm_files = [f for f in os.listdir(subdir_path) 
                 if f.startswith('Fdm_') and f.endswith('_introgress.pdf')]
    dsuite_files = sorted([f for f in os.listdir(subdir_path) 
                          if f.startswith('Dsuite_gomf_') and f.endswith('.pdf')])
    
    if not fdm_files or not dsuite_files:
        log(f"Skipping {subdir_name}: missing required PDFs")
        return False
    
    # Check if already processed
    out_name = f"Fdm_Dsuite_combined_{subdir_name}.pdf"
    out_path = os.path.join(subdir_path, out_name)
    if os.path.exists(out_path):
        log(f"Skipping {subdir_name}: already has combined PDF")
        return True
    
    # Get outgroup
    outgroup = get_outgroup(subdir_path)
    if not outgroup:
        log(f"WARNING: {subdir_name}: no mapping file found, using 'Outgroup'")
        outgroup = "Outgroup"
    
    # Parse sample names
    p1, p2, p3 = parse_samples(subdir_name)
    if not p1:
        log(f"ERROR: {subdir_name}: cannot parse sample names")
        return False
    
    log(f"Processing: {subdir_name}")
    log(f"  P1={p1}, P2={p2}, P3={p3}, Outgroup={outgroup}")
    
    try:
        # Open source PDFs
        fdm_path = os.path.join(subdir_path, fdm_files[0])
        fdm_doc = fitz.open(fdm_path)
        dsuite_docs = []
        for df in dsuite_files:
            dsuite_docs.append((df, fitz.open(os.path.join(subdir_path, df))))
        
        fdm_rect = fdm_doc[0].rect
        fdm_w, fdm_h = fdm_rect.width, fdm_rect.height
        ds_rect = dsuite_docs[0][1][0].rect
        ds_w, ds_h = ds_rect.width, ds_rect.height
        
        # Layout
        tree_w = 280
        gap_x = 20
        gap_y = 35
        margin = 10
        
        row1_w = tree_w + gap_x + fdm_w
        row2_w = ds_w * len(dsuite_files) + gap_x * (len(dsuite_files) - 1)
        canvas_w = max(row1_w, row2_w) + margin * 2
        
        row1_y = margin
        row2_y = row1_y + fdm_h + gap_y
        canvas_h = row2_y + ds_h + margin
        
        # Create output PDF
        out_doc = fitz.open()
        page = out_doc.new_page(width=canvas_w, height=canvas_h)
        
        # Draw tree
        draw_tree(page, margin, row1_y, tree_w, fdm_h, p1, p2, p3, outgroup)
        
        # Embed Fdm
        fdm_x = margin + tree_w + gap_x
        page.show_pdf_page(
            fitz.Rect(fdm_x, row1_y, fdm_x + fdm_w, row1_y + fdm_h),
            fdm_doc, 0, keep_proportion=True
        )
        
        # Embed Dsuite plots
        ds_total_w = ds_w * len(dsuite_files) + gap_x * (len(dsuite_files) - 1)
        ds_start_x = margin + (canvas_w - margin * 2 - ds_total_w) / 2
        
        for i, (df_name, ds_doc) in enumerate(dsuite_docs):
            ds_x = ds_start_x + i * (ds_w + gap_x)
            page.show_pdf_page(
                fitz.Rect(ds_x, row2_y, ds_x + ds_w, row2_y + ds_h),
                ds_doc, 0, keep_proportion=True
            )
        
        # Row labels
        page.insert_text(
            fitz.Point(margin + 5, row1_y + 16),
            "(A) Fdm introgression analysis",
            fontname="hebo", fontsize=12, color=(0, 0, 0)
        )
        page.insert_text(
            fitz.Point(margin + 5, row2_y + 16),
            "(B) GO enrichment of introgressed genes",
            fontname="hebo", fontsize=12, color=(0, 0, 0)
        )
        
        # Dsuite titles
        fs_title = 11
        title_map = {
            'Dsuite_gomf_0.03.pdf': ('(a)', 'Top 3% regions'),
            'Dsuite_gomf_0.05.pdf': ('(b)', 'Top 5% regions'),
        }
        
        for i, (df_name, ds_doc) in enumerate(dsuite_docs):
            if df_name in title_map:
                prefix, label = title_map[df_name]
                title_text = f"{prefix} {label}"
            else:
                # Extract percentage from filename
                m = re.search(r'(\d+\.?\d*)', df_name)
                pct = m.group(1) if m else '?'
                title_text = f"({chr(97+i)}) Top {pct}% regions"
            
            tw = fitz.get_text_length(title_text, fontname="hebo", fontsize=fs_title)
            ds_x = ds_start_x + i * (ds_w + gap_x)
            page.insert_text(
                fitz.Point(ds_x + ds_w / 2 - tw / 2, row2_y - 14),
                title_text, fontname="hebo", fontsize=fs_title, color=(0, 0, 0)
            )
        
        # Separator lines
        sc = (0.7, 0.7, 0.7)
        sw = 0.5
        
        def dl(x1, y1, x2, y2, c=sc, w=sw):
            page.draw_line(fitz.Point(x1, y1), fitz.Point(x2, y2), color=c, width=w)
        
        dl(margin + tree_w + gap_x / 2, row1_y, margin + tree_w + gap_x / 2, row1_y + fdm_h)
        dl(margin, row1_y + fdm_h + gap_y / 2, canvas_w - margin, row1_y + fdm_h + gap_y / 2)
        
        for i in range(len(dsuite_files) - 1):
            sep_x = ds_start_x + (i + 1) * ds_w + (i + 0.5) * gap_x
            dl(sep_x, row2_y, sep_x, row2_y + ds_h)
        
        # Save
        out_doc.save(out_path, garbage=0, deflate=True)
        out_doc.close()
        fdm_doc.close()
        for _, ds_doc in dsuite_docs:
            ds_doc.close()
        
        log(f"  -> Saved: {out_name}")
        return True
        
    except Exception as e:
        log(f"  ERROR: {e}")
        return False


# ============================================
# Main
# ============================================
if __name__ == '__main__':
    print("=" * 60)
    print("Batch combining Fdm + Dsuite figures for all subdirectories")
    print("=" * 60)
    
    all_dirs = sorted([d for d in os.listdir(base) 
                       if os.path.isdir(os.path.join(base, d)) 
                       and d not in SKIP_DIRS])
    
    success = 0
    skipped = 0
    failed = 0
    
    for d in all_dirs:
        print(f"\n--- {d} ---")
        result = process_subdir(d)
        if result is True:
            success += 1
        elif result is False:
            failed += 1
        else:
            skipped += 1
    
    print(f"\n{'=' * 60}")
    print(f"Summary: {success} succeeded, {skipped} skipped, {failed} failed")
    print(f"{'=' * 60}")