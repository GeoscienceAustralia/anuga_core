"""Schematic of the Yeh--Petroff dam-break flume.

The report used to include a scan of the figure from Silvester and Cleary,
which is not redistributable and was never committed, so the figure came out
blank. This draws the same geometry from the numbers the case actually uses
(see numerical_Yeh_Petroff.py).
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

L, W = 1.6, 0.61               # flume, m
X_DAM = 0.4                    # gate
H_UP, H_DOWN = 0.3, 0.01       # stage either side of the gate, m
X_COL, Y_COL, S_COL = 0.9, 0.25, 0.12   # square column


def setup(filename='Yeh_Petroff.png'):
    fig, (ax_p, ax_s) = plt.subplots(
        2, 1, figsize=(9, 5.2), gridspec_kw={'height_ratios': [2.0, 1.0]})

    # ---- plan view -------------------------------------------------------
    ax_p.add_patch(Rectangle((0, 0), X_DAM, W, fc='#cfe3f5', ec='none'))
    ax_p.add_patch(Rectangle((X_DAM, 0), L - X_DAM, W, fc='#f2f2f2', ec='none'))
    ax_p.add_patch(Rectangle((X_COL, Y_COL), S_COL, S_COL, fc='0.35', ec='k'))
    ax_p.plot([X_DAM, X_DAM], [0, W], 'k--', lw=1.6)
    ax_p.add_patch(Rectangle((0, 0), L, W, fc='none', ec='k', lw=1.4))

    ax_p.annotate('gate at $x = %.1f$ m' % X_DAM, xy=(X_DAM, W),
                  xytext=(X_DAM - 0.02, W + 0.06), ha='center', fontsize=9)
    ax_p.annotate('reservoir\n$h = %.2f$ m' % H_UP, xy=(0.2, W / 2),
                  ha='center', va='center', fontsize=9)
    ax_p.annotate('wet bed, $h = %.2f$ m' % H_DOWN, xy=(1.3, W / 2 - 0.12),
                  ha='center', va='center', fontsize=9)
    ax_p.annotate('%.2f m column' % S_COL,
                  xy=(X_COL + S_COL, Y_COL + S_COL),
                  xytext=(X_COL + 0.34, Y_COL + S_COL + 0.10), fontsize=9,
                  arrowprops=dict(arrowstyle='->', lw=1))

    # flume dimensions
    ax_p.annotate('', xy=(0, -0.05), xytext=(L, -0.05),
                  arrowprops=dict(arrowstyle='<->', lw=1))
    ax_p.text(L / 2, -0.10, '%.2f m' % L, ha='center', va='top', fontsize=9)
    ax_p.annotate('', xy=(L + 0.04, 0), xytext=(L + 0.04, W),
                  arrowprops=dict(arrowstyle='<->', lw=1))
    ax_p.text(L + 0.07, W / 2, '%.2f m' % W, ha='left', va='center',
              rotation=90, fontsize=9)

    ax_p.set_xlim(-0.05, L + 0.20); ax_p.set_ylim(-0.20, W + 0.20)
    ax_p.set_aspect('equal'); ax_p.axis('off')
    ax_p.set_title('Plan view', fontsize=10, loc='left')

    # ---- side view -------------------------------------------------------
    ax_s.add_patch(Rectangle((0, 0), X_DAM, H_UP, fc='#cfe3f5', ec='none'))
    ax_s.add_patch(Rectangle((X_DAM, 0), L - X_DAM, H_DOWN, fc='#cfe3f5', ec='none'))
    ax_s.add_patch(Rectangle((X_COL, 0), S_COL, 0.28, fc='0.35', ec='k'))
    ax_s.plot([X_DAM, X_DAM], [0, H_UP + 0.04], 'k--', lw=1.6)
    ax_s.plot([0, L], [0, 0], 'k-', lw=1.4)
    ax_s.annotate('%.2f m' % H_UP, xy=(0.2, H_UP / 2), ha='center',
                  va='center', fontsize=9)

    ax_s.set_xlim(-0.05, L + 0.05); ax_s.set_ylim(-0.03, 0.40)
    ax_s.set_aspect('equal'); ax_s.axis('off')
    ax_s.set_title('Side view', fontsize=10, loc='left')

    fig.tight_layout()
    fig.savefig(filename, dpi=130)
    plt.close(fig)


if __name__ == '__main__':
    setup()
