import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import io

def make_vortex_antivortex_gif(
        filename="Vortex-Antivortex_movie.gif",
        n_frames=60,
        nx=35, ny=50):
    """
    Generate an animated GIF of a vortex–antivortex pair advected in a
    uniform background flow (qualitatively matching the GIF you showed).
    """

    # ----- spatial grid (tall rectangle, like your GIF) -----
    x = np.linspace(-2.0, 2.0, nx)
    y = np.linspace(-3.0, 3.0, ny)
    X, Y = np.meshgrid(x, y)

    # vortex strength and background flow
    Gamma = 1.0         # circulation magnitude
    U0 = 0.4            # uniform flow to the right
    eps2 = 0.05**2      # core regularization to avoid singularities

    frames = []

    for k in range(n_frames):
        t = k / (n_frames - 1.0)

        # ----- vortex and antivortex positions vs time -----
        # one comes from upper-right towards centre/bottom-left
        x1 = 1.5 - 2.2 * t
        y1 = 1.8 - 2.2 * t

        # the other comes from mid-left towards centre/top-right
        x2 = -1.2 + 2.2 * t
        y2 = 0.2 - 1.2 * t

        # ----- velocity field from two point vortices + uniform flow -----
        dx1 = X - x1
        dy1 = Y - y1
        dx2 = X - x2
        dy2 = Y - y2

        r1sq = dx1**2 + dy1**2 + eps2
        r2sq = dx2**2 + dy2**2 + eps2

        # point-vortex velocities: u = -Γ/(2π) * (y - y0)/r², v = Γ/(2π) * (x - x0)/r²
        u = (U0
             - Gamma / (2*np.pi) * dy1 / r1sq
             + Gamma / (2*np.pi) * dy2 / r2sq)
        v = (Gamma / (2*np.pi) * dx1 / r1sq
             - Gamma / (2*np.pi) * dx2 / r2sq)

        speed = np.sqrt(u**2 + v**2)

        # ----- make one frame -----
        fig = plt.figure(figsize=(3, 4), dpi=200)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_axis_off()
        ax.set_aspect("equal")
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())

        # colored quiver (similar to the original GIF)
        ax.quiver(X, Y, u, v, speed,
                  cmap="viridis",
                  angles="xy",
                  scale_units="xy",
                  scale=12)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0)
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf))

    # ----- save animated GIF -----
    frames[0].save(
        filename,
        save_all=True,
        append_images=frames[1:],
        duration=80,   # ms per frame
        loop=0
    )
    print(f"Saved: {filename}")

if __name__ == "__main__":
    make_vortex_antivortex_gif()
