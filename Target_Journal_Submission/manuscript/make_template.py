"""Build an ACS-style reference.docx that pandoc uses for fonts/heading styles."""
from docx import Document
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
import os

d = Document()
# base body font
style = d.styles["Normal"]
style.font.name = "Times New Roman"; style.font.size = Pt(12)
style.paragraph_format.line_spacing = 2.0  # double-spaced for submission
style.paragraph_format.space_after = Pt(0)

def conf(name, size, bold=True, italic=False, before=12, after=6, color=None):
    s = d.styles[name]
    s.font.name = "Times New Roman"; s.font.size = Pt(size); s.font.bold = bold; s.font.italic = italic
    s.paragraph_format.space_before = Pt(before); s.paragraph_format.space_after = Pt(after)
    s.paragraph_format.line_spacing = 1.5
    if color: s.font.color.rgb = RGBColor(*color)

for h, sz in [("Heading 1", 13), ("Heading 2", 12), ("Heading 3", 12)]:
    conf(h, sz, color=(0, 0, 0))
conf("Title", 16, before=0, after=12)
# caption style
try:
    cap = d.styles["Caption"]; cap.font.name = "Times New Roman"; cap.font.size = Pt(10); cap.font.italic = True
except KeyError:
    pass
d.save(os.path.join(os.path.dirname(__file__), "assets", "acs_template.docx"))
print("template written")
