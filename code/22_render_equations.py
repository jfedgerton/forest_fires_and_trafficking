#!/usr/bin/env python3
"""
22_render_equations.py

Renders the deck's equations with real LaTeX and writes transparent white-on-
nothing PNGs for the black slides.

Why not native PowerPoint equations. omml.py (kept in this folder) injects real
OMML, and it validates, but LibreOffice silently drops a14:m when it renders,
so there is no way to check in this environment what the slide actually looks
like before the meeting. An image renders identically everywhere and can be
inspected. The cost is that the equations are no longer editable in PowerPoint;
if that matters more than verifiability, run omml.py instead.

Output: output/figures/eq_<key>.png, transparent, 600 dpi.
Run from code/.
"""

import os
import shutil
import subprocess
import tempfile

FIGDIR = "output/figures"
DPI = 600

EQUATIONS = {
    "es": r"e_s(T)\;=\;0.6108\,\exp\!\left(\frac{17.27\,T}{T+237.3}\right)"
          r"\quad[\mathrm{kPa}]",

    "vpd": r"\mathrm{VPD}\;=\;\max\!\left(0,\;e_s(T_{\mathrm{air}})"
           r"-e_s(T_{\mathrm{dew}})\right)",

    "exc3": r"\mathrm{VPD}^{\mathrm{exc}}\;=\;"
            r"\max\!\left(0,\;\mathrm{VPD}^{\max}-3.0\right)",

    "stage1": r"\mathrm{Fire}_{mt}\;=\;\pi\,\mathrm{VPD}_{mt}"
              r"\;+\;X_{mt}'\gamma\;+\;\alpha_m\;+\;\delta_t\;+\;u_{mt}",

    "stage2": r"Y_{mt}\;=\;\beta\,\widehat{\mathrm{Fire}}_{mt}"
              r"\;+\;X_{mt}'\gamma\;+\;\alpha_m\;+\;\delta_t\;+\;\varepsilon_{mt}",

    "rf": r"Y_{mt}\;=\;\rho\,\mathrm{VPD}_{mt}\;+\;\alpha_m\;+\;\delta_t"
          r"\;+\;\eta_{mt}",
}

PREAMBLE = r"""\documentclass[12pt]{article}
\usepackage[paperwidth=30cm,paperheight=6cm,margin=0pt]{geometry}
\usepackage{amsmath,amssymb}
\usepackage[T1]{fontenc}
\usepackage{xcolor}
\pagestyle{empty}
\begin{document}
\color{white}
\noindent$\displaystyle %s $
\end{document}
"""


def render(key, tex, outdir):
    with tempfile.TemporaryDirectory() as td:
        src = os.path.join(td, "eq.tex")
        with open(src, "w") as f:
            f.write(PREAMBLE % tex)
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                            "-output-directory", td, src],
                           capture_output=True, text=True)
        pdf = os.path.join(td, "eq.pdf")
        if r.returncode != 0 or not os.path.exists(pdf):
            tail = "\n".join(r.stdout.splitlines()[-25:])
            raise RuntimeError(f"pdflatex failed for {key}:\n{tail}")
        png = os.path.join(outdir, f"eq_{key}.png")
        # pdftocairo rather than ImageMagick: the default ImageMagick policy in
        # this image refuses to read PDF, and -transp gives a genuinely
        # transparent background so the white glyphs sit on the black slide.
        stem = os.path.join(td, "eq_out")
        c = subprocess.run(["pdftocairo", "-png", "-r", str(DPI), "-transp",
                            "-singlefile", pdf, stem],
                           capture_output=True, text=True)
        raw = stem + ".png"
        if c.returncode != 0 or not os.path.exists(raw):
            raise RuntimeError(f"pdftocairo failed for {key}:\n{c.stderr[-800:]}")
        # crop to the glyphs: the LaTeX page is deliberately oversized
        from PIL import Image
        im = Image.open(raw).convert("RGBA")
        bbox = im.getchannel("A").getbbox()
        if bbox is None:
            raise RuntimeError(f"{key} rendered empty")
        pad = int(DPI * 0.02)
        l, t, r_, b = bbox
        im = im.crop((max(0, l - pad), max(0, t - pad),
                      min(im.width, r_ + pad), min(im.height, b + pad)))
        im.save(png)
        return png


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    if shutil.which("pdflatex") is None:
        raise SystemExit("pdflatex not found")
    from PIL import Image
    # The printed ratios are pasted into EQ_RATIO in 21_build_deck_edgerton.js,
    # which places each equation by width and derives its height, so an equation
    # can never end up stretched on a slide.
    print("aspect ratios, for sizing the slide boxes:")
    for key, tex in EQUATIONS.items():
        p = render(key, tex, FIGDIR)
        w, h = Image.open(p).size
        print(f"  {key:8s} {w:5d} x {h:4d}  ratio {w/h:6.2f}  -> {p}")


if __name__ == "__main__":
    main()
