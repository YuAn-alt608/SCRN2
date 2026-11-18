"""
Shared plotting configuration: colour-blind-safe palette, Arial, 300 dpi vector export
Imported by every publication-grade figure in this deposit
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ============ Colour-blind-safe palette (modified Okabe-Ito, journal style) ============
# blue / orange / cyan / green / purple / yellow / grey / red, all separable for deuteranopes
C_BLUE   = "#0072B2"   # primary blue
C_ORANGE = "#E69F00"   # orange
C_CYAN   = "#56B4E9"   # cyan
C_GREEN  = "#009E73"   # green (up / positive direction)
C_PURPLE = "#CC79A7"   # purple
C_YELLOW = "#F0E442"   # yellow
C_GRAY   = "#999999"   # grey
C_RED    = "#D55E00"   # red-brown (down / caution; not pure red, colour-blind-safe)
C_DARK   = "#333333"

PALETTE = [C_BLUE, C_ORANGE, C_CYAN, C_GREEN, C_PURPLE, C_YELLOW, C_GRAY, C_RED]

# ============ Fonts and canvas ============
rcParams["font.family"] = "sans-serif"
# Arial first; CJK faces kept only as a fallback for stray non-ASCII labels
rcParams["font.sans-serif"] = ["Arial", "Microsoft YaHei", "SimHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 9
rcParams["axes.linewidth"] = 0.8
rcParams["axes.edgecolor"] = C_DARK
rcParams["xtick.direction"] = "out"
rcParams["ytick.direction"] = "out"
rcParams["axes.spines.top"] = False
rcParams["axes.spines.right"] = False

DPI = 300

def panel_label(ax, text, x=-0.12, y=1.08, size=14, weight="bold"):
    """Place a panel letter (A/B/C ...) at a fixed offset from the axes."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size,
            fontweight=weight, va="top", ha="right", color=C_DARK)

def save(fig, path):
    """Save as 300 dpi PNG and, from the same figure, as vector PDF."""
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    pdf = path.replace(".png", ".pdf")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"saved: {path}  +  {pdf}")
