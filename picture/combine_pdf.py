from pypdf import PdfWriter, PdfReader, Transformation
import glob, os

pdf_dir = r"E:\Project\Rosacea\Rosaceae_cytonuclear\03-final_result_to_publish\12-Get_introgress_gene\orthofinder\KEGG_GO_rich\SRR10377315_SRR15691173_SRR24154117\comb"
out_pdf = r"E:\Project\Rosacea\Rosaceae_cytonuclear\03-final_result_to_publish\12-Get_introgress_gene\orthofinder\KEGG_GO_rich\SRR10377315_SRR15691173_SRR24154117\comb\combined_output.pdf"
COLS = 3
ROWS = 1
GAP = 5

pdfs = sorted(glob.glob(os.path.join(pdf_dir, "*.pdf")))[:COLS * ROWS]

if not pdfs:
    raise RuntimeError("No PDF files found")

pages = [PdfReader(p).pages[0] for p in pdfs]

cell_w = float(pages[0].mediabox.width)
cell_h = float(pages[0].mediabox.height)

canvas_w = COLS * cell_w + (COLS - 1) * GAP
canvas_h = ROWS * cell_h + (ROWS - 1) * GAP

writer = PdfWriter()
canvas_page = writer.add_blank_page(width=canvas_w, height=canvas_h)

for idx, page in enumerate(pages):
    row = idx // COLS
    col = idx % COLS
    x = col * (cell_w + GAP)
    y = canvas_h - (row + 1) * cell_h - row * GAP

    # ✅ pypdf 6.x 的正确写法
    canvas_page.merge_transformed_page(
        page,
        Transformation().translate(x, y)
    )

writer.write(out_pdf)
print(f"Saved: {out_pdf}")