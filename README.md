# AlignTraj

`AlignTraj.m` loads a trajectory dataset, aligns all trajectories to a common solute-based reference, puts everything into one consistent 3D frame (and optionally makes a movie of the result to check).
The main point is simple: different trajectories can have arbitrary overall translation and rotation, which makes direct comparison and anisotropy `S2` simulations messy. This script removes that rigid-body motion so the remaining differences are the actual structural dynamics.
The files are commented. The test data set (`test_traj_dataset.mat`) is just some truncated trajectories of a solute (38 atoms) in solvent, overall 500 atoms. The original system had more atoms, but I trimmed it down to keep the example compact and GitHub-friendly.
We use the first time frame of the first trajectory as a reference, and align the solute part of every trajectory to that common structure with a mass-weighted rigid transform. That gives each trajectory one rotation and one translation, which are then applied to all time frames of that trajectory, so the internal motion is preserved.
After that, the script defines one final global frame from the aligned solute reference by centering it at its mass-weighted center of mass and orienting it along its principal axes. This just makes the final coordinates look clean and consistent, so all trajectories are centered and pointing the same way.
The main output is `xyzt_final`, which has the same size as the input `xyzt` but now lives in one shared frame. 
If you comment out the `return` line near the end, the script also writes a small MP4 movie (`ensemble_movie.mp4`) showing how the aligned solute ensemble evolves in time.
