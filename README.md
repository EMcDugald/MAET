This project contains methods to solve the forward problem of Magneto Acousto Electric Tomography (MAET).

MAET is a medical imaging modality, which recovers the conductivity of tissue as a contrast agent. 
A MAET image is obtained by placing the tissue of interest in a conductive medium, and coupling acoustic waves with a constant magnetic field to generate electrical signals.
Let $M(t)$ denote our measurement at time $t$. If we are using a pair of electrodes at points $z_1, z_2$ on the boundary of the apparatus, our measurements are given by 
$$
\begin{align}
M(t) = u(t,z_1)-u(t,z_2)
\end{align}
$$
where $u(t,x)$ denotes the electric potential of the tissue.

