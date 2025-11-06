# CVVR

Assignment Week 2 for **Scientific Visualization and Virtual Reality** at the **University of Amsterdam (UvA)**.

The task was to create a **2D Steady-State Heat Conduction** domain and compute the equilibrium temperature distribution.  
The simulation models a 9×9 m square plate discretized into 241×241 blocks — balancing computational efficiency with sufficient spatial resolution — with an inner 3×3 m square fixed at 212 °F.  
The bottom 4 m of the outer wall is held at 32 °F, increasing linearly up to 100 °F at the top boundary.  
The initial guess for all unknown cells was 90 °F, chosen as a balanced starting point, though any value would converge to equilibrium.

The solver uses the **Jacobi relaxation method** to solve the **Laplace equation** on a half-domain with left–right symmetry to reduce computation time.  
Snapshots are saved every 500 epochs and mirrored to reconstruct the full temperature field.  
Visualization of the temperature evolution and steady-state result was created in **ParaView** and exported as an **MP4/AVI** animation.

### Final Converged Temperature Field
The steady-state solver converged to a tolerance of **Δ = 1×10⁻³** within **8664 epochs**.  
The final temperature distribution is shown below:
![final](media/final.png)

## Temperature Progression
The animation below shows the steady-state solver approaching convergence.  
Each frame corresponds to one snapshot every **500 iterations** (18 frames total).

![Temperature evolution](media/temperature_progression.gif)

[ Download full animation (AVI)](media/temperature_progression.avi)
