import numpy as np
import csv
from pathlib import Path

# --------------------------
# Problem construction (half-domain with left-right symmetry)
# --------------------------
def build_problem_half(
    N=241, L=9.0, L_inner=3.0,
    T_top=100.0, T_bottom=32.0, T_inner=212.0, init_guess=90.0,
    bath_height=None,
    bath_fraction=4.0/9.0
):
    Ny = N
    Nx_full = N
    h_y = L / (Ny - 1)
    Nx_half = Nx_full // 2 + 1
    h_x = (L/2.0) / (Nx_half - 1) if Nx_half > 1 else 0.0

    T = np.full((Ny, Nx_half), init_guess, dtype=np.float64)
    fixed = np.zeros_like(T, dtype=bool)

    y = np.linspace(0.0, L, Ny)

    # Dirichlet: bottom/top
    T[0, :]  = T_bottom; fixed[0, :]  = True
    T[-1, :] = T_top;    fixed[-1, :] = True

    # Outer wall profile: 32°F up to bath height, then linear to 100°F at top
    y0 = bath_height if bath_height is not None else bath_fraction * L
    y0 = np.clip(y0, 0.0, L)
    side_profile = np.where(
        y <= y0,
        T_bottom,
        T_bottom + (T_top - T_bottom) * ((y - y0) / (L - y0))
    )
    T[:, -1] = side_profile
    fixed[:, -1] = True

    # Inner hot square (centered)
    half = L / 2.0
    inner_half = L_inner / 2.0
    j_low  = int(round((half - inner_half) / h_y))
    j_high = int(round((half + inner_half) / h_y))
    i_left, i_right = 0, int(round(inner_half / h_x))
    T[j_low:j_high+1, i_left:i_right+1] = T_inner
    fixed[j_low:j_high+1, i_left:i_right+1] = True

    return T, fixed, h_y, h_x, L

# --------------------------
# Mirror half -> full
# --------------------------
def mirror_full_from_half(T_half):
    # exclude centerline when mirroring to avoid duplicate column
    left = T_half[:, 1:][:, ::-1]
    return np.hstack([left, T_half])

# --------------------------
# Export full structured CSV with physical coords (0..9 m)
# --------------------------
def export_structured_csv_full(T_full, L, out_path):
    """
    Writes a structured table for ParaView:
      columns: i,j,k,x,y,Temperature
    X, Y are scaled from 0..L
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    Ny, Nx = T_full.shape
    h_x = L / (Nx - 1)
    h_y = L / (Ny - 1)

    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["i", "j", "k", "x", "y", "Temperature"])
        for j in range(Ny):
            y = j * h_y
            for i in range(Nx):
                x = i * h_x
                w.writerow([i, j, 0, x, y, float(T_full[j, i])])

# --------------------------
# Jacobi solver on half with periodic snapshots
# --------------------------
def jacobi_laplace_half_with_snapshots(
    T, fixed,
    h_y, h_x, L,
    delta=1e-3,
    max_iters=2_000_000,
    report_every=500,
    snap_every=500,
    out_dir="paradise_series",
    prefix="full_structured"
):
    T = T.copy()
    Ny, Nx = T.shape

    # which cells get updated
    var = ~fixed
    var[:, 0]  = False   # symmetry line
    var[0, :]  = False   # bottom
    var[-1, :] = False   # top
    var[:, -1] = False   # right wall

    epochs = 0
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    while True:
        T_new = T.copy()

        avg = 0.25 * (T[:-2, 1:-1] + T[2:, 1:-1] + T[1:-1, :-2] + T[1:-1, 2:])
        interior = (slice(1, -1), slice(1, -1))
        mask = var[1:-1, 1:-1]
        T_new[interior][mask] = avg[mask]

        # Neumann (∂T/∂x = 0) at centerline
        T_new[:, 0] = T_new[:, 1]

        # keep Dirichlet fixed
        T_new[fixed] = T[fixed]

        # convergence metric
        max_change = np.max(np.abs(T_new - T)[var]) if np.any(var) else 0.0

        T = T_new
        epochs += 1

        if epochs % report_every == 0:
            print(f"epoch {epochs} | max_change={max_change:.6f}")

        if epochs % snap_every == 0:
            T_full = mirror_full_from_half(T)
            export_structured_csv_full(T_full, L, out_dir / f"{prefix}_{epochs:06d}.csv")

        if max_change < delta or epochs >= max_iters:
            T_full = mirror_full_from_half(T)
            export_structured_csv_full(T_full, L, out_dir / f"{prefix}_final.csv")
            return T, epochs, max_change

# --------------------------
# Main
# --------------------------
def main():
    T0_half, fixed_half, h_y, h_x, L = build_problem_half(
        N=241, L=9.0, L_inner=3.0,
        T_top=100.0, T_bottom=32.0, T_inner=212.0, init_guess=90.0
    )

    T_half, epochs, last_change = jacobi_laplace_half_with_snapshots(
        T0_half, fixed_half,
        h_y=h_y, h_x=h_x, L=L,
        delta=1e-3,
        report_every=500,
        snap_every=500,
        out_dir="paradise_series",
        prefix="full_structured"
    )

    print(f"Converged in {epochs} iterations, max change = {last_change:.6f} °F")
    print("CSV time series written to ./paradise_series")

if __name__ == "__main__":
    main()
