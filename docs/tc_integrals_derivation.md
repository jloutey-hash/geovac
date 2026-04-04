# Transcorrelated Integrals in Level 3 Hyperspherical Coordinates

**Track BX-2, Sub-agent 1 | April 3, 2026**

## 1. Conventions and Setup

### 1.1 Jastrow Convention

We use **Convention A**: the exact wavefunction is related to the smooth TC wavefunction by

$$\psi_{\text{exact}} = e^{J} \, \psi_{\text{TC}}$$

and the TC Hamiltonian is obtained by the similarity transformation

$$H_{\text{TC}} = e^{-J} \, H \, e^{J}.$$

The Kato cusp condition requires $\psi \sim 1 + \tfrac{1}{2} r_{12}$ near $r_{12} = 0$, so $e^J \sim 1 + \tfrac{1}{2} r_{12}$, giving

$$J = \tfrac{1}{2} r_{12}$$

(positive sign). This is the standard choice in the TC literature (Ten-no 2004, Luo 2010, Dobrautz et al. 2019).

**Sign check:** With $J = +\tfrac{1}{2} r_{12}$, the cusp factor $e^{J}$ has positive cusp slope $\partial e^J / \partial r_{12}|_0 = +1/2$, matching Kato.

### 1.2 Level 3 Coordinates

From Paper 13, the two-electron system in hyperspherical coordinates uses:

- $R$ = hyperradius, $R^2 = r_1^2 + r_2^2$
- $\alpha \in [0, \pi/4]$ = hyperangle, $\tan\alpha = r_2/r_1$ (singlet symmetry restricts to $[0, \pi/4]$)
- $\theta_{12} \in [0, \pi]$ = interelectron angle
- Three Euler angles for total angular momentum (decouple for $L=0$)

The interelectron distance is:
$$r_{12} = R\sqrt{u}, \qquad u \equiv 1 - \sin 2\alpha \cos\theta_{12}$$

At the coalescence point $r_{12} = 0$: $\alpha = \pi/4$, $\theta_{12} = 0$ (where $u = 0$).

### 1.3 Kinetic Energy

The full 6D Laplacian (Paper 13, Eq. 4):
$$\nabla_1^2 + \nabla_2^2 = \frac{\partial^2}{\partial R^2} + \frac{5}{R}\frac{\partial}{\partial R} + \frac{1}{R^2}\Lambda^2$$

where $\Lambda^2$ is the grand angular momentum (Laplace-Beltrami on $S^5$, Eq. 5):
$$\Lambda^2 = -\frac{1}{\sin^2\alpha\cos^2\alpha}\frac{\partial}{\partial\alpha}\left(\sin^2\alpha\cos^2\alpha\frac{\partial}{\partial\alpha}\right) + \frac{\hat{l}_1^2}{\cos^2\alpha} + \frac{\hat{l}_2^2}{\sin^2\alpha}$$

The kinetic energy is $T = -\frac{1}{2}(\nabla_1^2 + \nabla_2^2)$.

### 1.4 Angular Basis

The Level 3 angular basis for $L=0$ singlet is (Paper 13 Sec III, `algebraic_angular.py`):
$$\phi_{\nu,l}(\alpha, \theta_{12}) \propto (\sin\alpha\cos\alpha)^{l+1} C_k^{l+1}(\cos 2\alpha) \cdot P_l(\cos\theta_{12})$$

where $\nu = 2(l+k+1)$ is the hyperspherical quantum number, $k = 0, 2, 4, \ldots$ (singlet), and $C_k^\lambda$ is the Gegenbauer polynomial. The free eigenvalues are $\mu_{\text{free}}(\nu) = \nu(\nu+4)/2$.

---

## 2. BCH Expansion

Since $J = \tfrac{1}{2} r_{12}$ is a function of coordinates only, $[V, J] = 0$. The BCH expansion is:

$$H_{\text{TC}} = e^{-J} H e^{J} = H + [H, J] + \tfrac{1}{2}[[H, J], J] + \cdots$$

Since $[V, J] = 0$, only the kinetic part contributes:
$$H_{\text{TC}} = T + V + [T, J] + \tfrac{1}{2}[[T, J], J]$$

For two electrons with $T = T_1 + T_2$, $T_i = -\frac{1}{2}\nabla_i^2$:

$$[T, J] = [T_1 + T_2, J] = -\sum_{i=1}^{2} \left[(\nabla_i J) \cdot \nabla_i + \tfrac{1}{2}\nabla_i^2 J\right]$$

$$\tfrac{1}{2}[[T, J], J] = \tfrac{1}{2}\sum_{i=1}^{2} |\nabla_i J|^2$$

The BCH series **terminates exactly at second order** for a linear Jastrow $J = \tfrac{1}{2}r_{12}$, because $[[T, J], J]$ involves $|\nabla_i J|^2$ which is a multiplicative operator, and $[[[T,J],J],J] = 0$ since $\nabla_i(|\nabla_j J|^2)$ contracted with $\nabla_j J$ gives zero when $|\nabla_i r_{12}|^2 = 1$ (constant).

---

## 3. Term-by-Term Derivation

### 3.1 Gradient of J

$$\nabla_1 J = \tfrac{1}{2}\nabla_1 r_{12} = \tfrac{1}{2}\frac{\mathbf{r}_1 - \mathbf{r}_2}{r_{12}} = \tfrac{1}{2}\hat{r}_{12}$$

$$\nabla_2 J = \tfrac{1}{2}\nabla_2 r_{12} = -\tfrac{1}{2}\hat{r}_{12}$$

These are unit vectors (up to the factor 1/2), well-defined everywhere except $r_{12} = 0$.

### 3.2 Term 3 (Easiest): The Double Commutator $\frac{1}{2}[[T,J],J]$

$$\tfrac{1}{2}[[T,J],J] = \tfrac{1}{2}\sum_{i=1}^{2}|\nabla_i J|^2 = \tfrac{1}{2}\left(\tfrac{1}{4} + \tfrac{1}{4}\right) = \frac{1}{4}$$

**This is a CONSTANT.** It shifts all eigenvalues by $+1/4$ uniformly.

**Coordinate invariance check:** $|\nabla_i r_{12}|^2$ is a scalar (coordinate-independent). In Cartesian coordinates, $|\nabla_1 r_{12}|^2 = |\hat{r}_{12}|^2 = 1$. This is true in any coordinate system. In hyperspherical coordinates, the metric factors in $|\nabla r_{12}|^2$ produce the same scalar 1.

**Result:** $\frac{1}{2}[[T,J],J] = \frac{1}{4}$ (trivial energy shift, no matrix elements needed).

### 3.3 Term 2: The Laplacian $-\frac{1}{2}\sum_i \nabla_i^2 J$

For a single electron in 3D, the Laplacian of the distance function is:
$$\nabla_1^2 r_{12} = \frac{2}{r_{12}}$$

This is the standard result: for $f(\mathbf{r}) = |\mathbf{r} - \mathbf{r}_0|$ in 3D, $\nabla^2 f = 2/f$.

Similarly $\nabla_2^2 r_{12} = 2/r_{12}$.

Therefore:
$$-\frac{1}{2}\sum_{i=1}^{2}\nabla_i^2 J = -\frac{1}{2}\cdot\frac{1}{2}\cdot\frac{4}{r_{12}} = -\frac{1}{r_{12}}$$

**This exactly cancels the electron-electron repulsion $V_{ee} = +1/r_{12}$!**

Combined:
$$V_{ee} + \left(-\frac{1}{2}\sum_i \nabla_i^2 J\right) = \frac{1}{r_{12}} - \frac{1}{r_{12}} = 0$$

**The $1/r_{12}$ singularity is completely removed from $H_{\text{TC}}$.**

### 3.4 Term 1: The Gradient Term $-\sum_i (\nabla_i J)\cdot\nabla_i$

This is the only non-trivial TC correction term. It is a first-order differential operator:

$$G \equiv -\sum_{i=1}^{2}(\nabla_i J)\cdot\nabla_i = -\frac{1}{2}\hat{r}_{12}\cdot\nabla_1 + \frac{1}{2}\hat{r}_{12}\cdot\nabla_2$$

$$= -\frac{1}{2}(\hat{r}_{12}\cdot\nabla_1 - \hat{r}_{12}\cdot\nabla_2) = -\frac{1}{2}\frac{\partial}{\partial r_{12}}$$

where $\partial/\partial r_{12}$ is the derivative along the interelectron separation direction. This is the **relative momentum operator** projected onto the interelectron axis.

More precisely, defining the relative coordinate $\mathbf{s} = \mathbf{r}_1 - \mathbf{r}_2$ and center-of-mass $\mathbf{S} = (\mathbf{r}_1 + \mathbf{r}_2)/2$:

$$G = -\frac{1}{2}\hat{s}\cdot\nabla_s = -\frac{1}{2}\frac{\partial}{\partial s}$$

where $s = |\mathbf{s}| = r_{12}$ and $\partial/\partial s$ is the radial derivative in relative coordinates.

---

## 4. Full TC Hamiltonian

Collecting all terms:

$$\boxed{H_{\text{TC}} = T + V_{\text{nuc}} + G + \frac{1}{4}}$$

where:
- $T = -\frac{1}{2}(\nabla_1^2 + \nabla_2^2)$ — original kinetic energy
- $V_{\text{nuc}} = -Z/r_1 - Z/r_2$ — nuclear attraction (unchanged)
- $G = -\frac{1}{2}\frac{\partial}{\partial r_{12}}$ — the TC gradient operator (first-order, non-Hermitian)
- $+1/4$ — constant energy shift

**The $V_{ee} = 1/r_{12}$ term has been completely eliminated.**

This is the central result: the TC transformation with the exact Kato Jastrow removes the electron-electron Coulomb singularity entirely. The price is a first-order differential operator $G$ that makes $H_{\text{TC}}$ non-Hermitian.

---

## 5. The Gradient Operator in Level 3 Coordinates

### 5.1 Derivatives of $r_{12}$ in Hyperspherical Coordinates

With $r_{12} = R\sqrt{u}$, $u = 1 - \sin 2\alpha\cos\theta_{12}$:

$$\frac{\partial r_{12}}{\partial R} = \sqrt{u} = \frac{r_{12}}{R}$$

$$\frac{\partial r_{12}}{\partial\alpha} = R\frac{\partial\sqrt{u}}{\partial\alpha} = R\frac{-2\cos 2\alpha\cos\theta_{12}}{2\sqrt{u}} = -\frac{R\cos 2\alpha\cos\theta_{12}}{\sqrt{u}}$$

$$\frac{\partial r_{12}}{\partial\theta_{12}} = R\frac{\partial\sqrt{u}}{\partial\theta_{12}} = R\frac{\sin 2\alpha\sin\theta_{12}}{2\sqrt{u}}$$

### 5.2 The Operator $\hat{r}_{12}\cdot\nabla_1 - \hat{r}_{12}\cdot\nabla_2$

Rather than expressing the 6D gradient in hyperspherical coordinates (which involves the complicated $S^5$ metric), we use the identity from Section 3.4:

$$\hat{r}_{12}\cdot\nabla_1 - \hat{r}_{12}\cdot\nabla_2 = 2\frac{\partial}{\partial r_{12}}$$

where $\partial/\partial r_{12}$ is the derivative with respect to the interelectron distance, holding all other coordinates fixed. In the $(R, \alpha, \theta_{12})$ system, this can be expressed using the chain rule. Since $r_{12} = R\sqrt{u(\alpha, \theta_{12})}$, the coordinate $r_{12}$ is a function of all three coordinates $R, \alpha, \theta_{12}$. We need the inverse: how does $\partial/\partial r_{12}$ act on functions of $(R, \alpha, \theta_{12})$?

**Change of variable approach.** Consider holding $R$ and one other angular variable fixed. At fixed $R$ and fixed $\alpha$, we have $r_{12}^2 = R^2(1 - \sin 2\alpha\cos\theta_{12})$, so $\theta_{12}$ and $r_{12}$ are in 1-to-1 correspondence. Then:

$$\left.\frac{\partial}{\partial r_{12}}\right|_{R,\alpha} = \frac{1}{\partial r_{12}/\partial\theta_{12}}\frac{\partial}{\partial\theta_{12}} = \frac{2\sqrt{u}}{R\sin 2\alpha\sin\theta_{12}}\frac{\partial}{\partial\theta_{12}}$$

However, $\partial/\partial r_{12}$ at fixed $R$ and $\alpha$ is NOT the same as the physical $\partial/\partial r_{12}$ appearing in the TC operator. The physical derivative is at fixed center-of-mass and fixed angular orientation of the interelectron vector.

### 5.3 Direct 6D Gradient Approach

It is cleaner to work directly in the 6D formalism. The TC gradient operator is:

$$G = -\frac{1}{2}(\hat{r}_{12}\cdot\nabla_1 - \hat{r}_{12}\cdot\nabla_2) = -\frac{\partial}{\partial r_{12}}$$

In the hyperspherical $(R, \alpha, \theta_{12}, \text{Euler angles})$ system, we express $G$ using the 6D metric. The 6D line element in hyperspherical coordinates is:

$$ds^2 = dR^2 + R^2\left(d\alpha^2 + \cos^2\alpha\, d\Omega_1^2 + \sin^2\alpha\, d\Omega_2^2\right)$$

where $d\Omega_i^2$ are the solid angle elements for each electron. For $L = 0$, after integrating over the three Euler angles, we work in the $(R, \alpha, \theta_{12})$ subspace.

The operator $G$ in this subspace can be decomposed using the chain rule. Since $r_{12} = R\sqrt{u(\alpha, \theta_{12})}$:

$$G = -\frac{\partial}{\partial r_{12}} = -\frac{1}{\nabla_s r_{12}\cdot\nabla_s r_{12}}\nabla_s r_{12}\cdot\nabla_s$$

where $\nabla_s$ is the gradient in the relative coordinate $\mathbf{s} = \mathbf{r}_1 - \mathbf{r}_2$ (3D). Since $|\nabla_s r_{12}|^2 = 1$ (the relative-coordinate gradient of the relative distance is a unit vector):

$$G = -\nabla_s r_{12}\cdot\nabla_s = -\hat{s}\cdot\nabla_s = -\frac{\partial}{\partial s}$$

This is the **radial derivative in relative coordinates**. To express it in $(R, \alpha, \theta_{12})$, we use the Jacobian. The relative and center-of-mass coordinates are:

$$\mathbf{s} = \mathbf{r}_1 - \mathbf{r}_2, \qquad \mathbf{S} = \tfrac{1}{2}(\mathbf{r}_1 + \mathbf{r}_2)$$

In hyperspherical coordinates, the interelectron distance depends on $(R, \alpha, \theta_{12})$ as $r_{12} = R\sqrt{u}$. However, the physical derivative $\partial/\partial r_{12}$ holding the center-of-mass and angular orientation of $\hat{r}_{12}$ fixed does NOT correspond simply to a derivative in $(R, \alpha, \theta_{12})$ holding two of them fixed.

### 5.4 Resolution: Matrix Elements via Projection

The cleanest approach for implementation is to compute the matrix elements of $G$ directly, using the angular basis. For the angular eigenvalue problem at fixed $R$, the TC correction modifies the angular Hamiltonian. We need:

$$\langle\phi_{\nu,l}|G|\phi_{\nu',l'}\rangle$$

where the inner product is over the $S^5$ angular coordinates. Since $G = -\partial/\partial r_{12}$ is a first-order differential operator, we can compute its matrix elements by expanding $G$ in the angular variables.

**Key insight:** At fixed $R$, the operator $G$ acts only on the angular part. We can write the action of $G$ on a function $f(R, \alpha, \theta_{12})$ at fixed $R$ as a differential operator in $(\alpha, \theta_{12})$.

Using the chain rule at fixed $R$, $\hat{\Omega}_1$, and $\hat{\Omega}_2$ orientations (which is what the L=0 partial wave reduction does):

$$G\Big|_R = -\frac{\partial r_{12}}{\partial\alpha}\Big|_{R,\theta_{12}}\cdot g^{\alpha\alpha}\frac{\partial}{\partial\alpha} - \frac{\partial r_{12}}{\partial\theta_{12}}\Big|_{R,\alpha}\cdot g^{\theta_{12}\theta_{12}}\frac{\partial}{\partial\theta_{12}} + \ldots$$

Wait — this is the wrong way to think about it. $G$ is defined as $-\hat{s}\cdot\nabla_s$ in the 3D relative space. After the $L=0$ partial wave expansion, we integrate out the angular part of $\hat{s}$ and are left with the **radial** part in relative coordinates, which is $-\partial/\partial s$ where $s = r_{12}$.

For the angular eigenvalue problem at fixed $R$, what matters is the angular projection. Let us re-derive from scratch.

### 5.5 Derivation via the L=0 Partial Wave Expansion

The full 6D wavefunction for $L=0$ singlet is:
$$\Psi(\mathbf{r}_1, \mathbf{r}_2) = R^{-5/2}\sum_l f_l(R, \alpha)\frac{P_l(\cos\theta_{12})}{\sqrt{4\pi/(2l+1)}}$$

The TC gradient operator is (using the original Cartesian form):
$$G = -\frac{1}{2}\hat{r}_{12}\cdot(\nabla_1 - \nabla_2) = -\hat{s}\cdot\nabla_s$$

In spherical coordinates $(s, \hat{s})$ for the relative vector, where $s = r_{12}$:
$$\hat{s}\cdot\nabla_s = \frac{\partial}{\partial s} + \text{(angular terms in $\hat{s}$ that vanish for L=0 after averaging)}$$

Actually, for $L=0$, after integrating over the overall orientation, $\hat{s}\cdot\nabla_s$ acting on a function of $(R, \alpha, \theta_{12})$ becomes:

$$G = -\frac{\partial}{\partial s}\bigg|_{\hat{s}, \mathbf{S}} = -\frac{\partial}{\partial r_{12}}\bigg|_{\text{fixed CM and orientation}}$$

To convert this to $(R, \alpha)$ coordinates, we use the relations (at fixed $\theta_{12}$ and Euler angles):

$$R^2 = r_1^2 + r_2^2, \quad \tan\alpha = r_2/r_1$$
$$r_1 = R\cos\alpha, \quad r_2 = R\sin\alpha$$
$$r_{12}^2 = r_1^2 + r_2^2 - 2r_1r_2\cos\theta_{12} = R^2(1 - \sin 2\alpha\cos\theta_{12})$$

The Jacobian for the transformation $(r_1, r_2) \to (R, \alpha)$ at fixed angular variables gives:

$$\frac{\partial}{\partial r_{12}} = \frac{\partial R}{\partial r_{12}}\frac{\partial}{\partial R} + \frac{\partial\alpha}{\partial r_{12}}\frac{\partial}{\partial\alpha} + \frac{\partial\theta_{12}}{\partial r_{12}}\frac{\partial}{\partial\theta_{12}}$$

However, these partial derivatives depend on what is held fixed, which differs between the relative/CM coordinates and the hyperspherical coordinates. This makes a direct analytical transformation unwieldy.

### 5.6 Practical Resolution: Matrix Elements via Integration by Parts

For numerical implementation, the most robust approach is to compute $\langle\phi_i|G|\phi_j\rangle$ using integration by parts in the natural angular coordinates. The key identity is:

$$\langle\phi_i|G|\phi_j\rangle = -\int_{S^5} \phi_i^*\left(\frac{\partial}{\partial r_{12}}\right)\phi_j\, d\Omega_5$$

where the integral is over $S^5$ with the appropriate Jacobian measure $\sin^2\alpha\cos^2\alpha\sin\theta_{12}\,d\alpha\,d\theta_{12}\,d(\text{Euler})$.

At fixed $R$, the derivative $\partial/\partial r_{12}$ acting on the angular function $\phi_j(\alpha, \theta_{12})$ can be computed as:

$$\frac{\partial\phi_j}{\partial r_{12}} = \frac{\partial\phi_j}{\partial\alpha}\frac{\partial\alpha}{\partial r_{12}} + \frac{\partial\phi_j}{\partial\theta_{12}}\frac{\partial\theta_{12}}{\partial r_{12}}$$

where at fixed $R$ (and fixed CM direction):

$$\frac{\partial\alpha}{\partial r_{12}}\bigg|_{R} \quad \text{and} \quad \frac{\partial\theta_{12}}{\partial r_{12}}\bigg|_{R}$$

**These are not well-defined** without specifying what else is held constant besides $R$. This is the fundamental difficulty: the TC operator $\partial/\partial r_{12}$ mixes the hyperspherical coordinates in a way that depends on what is held constant.

### 5.7 Resolution via Explicit 6D Gradient

The cleanest derivation uses the explicit 6D gradient in hyperspherical coordinates. For a scalar function $f(\mathbf{r}_1, \mathbf{r}_2)$ that depends only on $r_{12} = |\mathbf{r}_1 - \mathbf{r}_2|$, the 6D gradient is:

$$\nabla_{6D} f(r_{12}) = f'(r_{12})\nabla_{6D} r_{12}$$

The TC operator is:
$$G = -\frac{1}{2}(\nabla_1 r_{12}\cdot\nabla_1 + \nabla_2 r_{12}\cdot\nabla_2) + \frac{1}{2}(\nabla_2 r_{12}\cdot\nabla_2 - \nabla_1 r_{12}\cdot\nabla_1)$$

Wait, let me be more careful. We have:
$$G = -(\nabla_1 J)\cdot\nabla_1 - (\nabla_2 J)\cdot\nabla_2 = -\frac{1}{2}\hat{r}_{12}\cdot\nabla_1 + \frac{1}{2}\hat{r}_{12}\cdot\nabla_2$$

In the 6D space with coordinates $q = (R, \alpha, \theta_{12}, \phi_1, \theta_1', \phi_2')$ (the last four being Euler-type angles), the gradient of $r_{12}$ has components:

$$(\nabla_{6D} r_{12})^R = \frac{\partial r_{12}}{\partial R} = \sqrt{u}$$

$$(\nabla_{6D} r_{12})^\alpha = g^{\alpha\alpha}\frac{\partial r_{12}}{\partial\alpha} = \frac{1}{R^2}\cdot\frac{-R\cos 2\alpha\cos\theta_{12}}{\sqrt{u}}$$

$$(\nabla_{6D} r_{12})^{\theta_{12}} = g^{\theta_{12}\theta_{12}}\frac{\partial r_{12}}{\partial\theta_{12}}$$

The metric factor for $\theta_{12}$ in the $L=0$ reduced problem requires knowledge of the reduced metric after integrating out the Euler angles. For the $L=0$ sector, the effective metric on the $(R, \alpha, \theta_{12})$ subspace is:

$$ds^2_{L=0} = dR^2 + R^2\,d\alpha^2 + R^2\sin^2\alpha\cos^2\alpha\,d\theta_{12}^2$$

(The last term comes from the dot product $\hat{r}_1\cdot\hat{r}_2 = \cos\theta_{12}$, with the metric factor being $r_1 r_2 = R^2\sin\alpha\cos\alpha$ from the geometry.)

**Actually this metric is not quite right** — the $\theta_{12}$ metric factor should account for the solid angle integration properly. The correct reduced line element in the $(R, \alpha, \theta_{12})$ subspace after $L=0$ projection involves the Jacobian $J = R^5\sin^2\alpha\cos^2\alpha\sin\theta_{12}$.

Rather than pursuing this further (which requires careful differential geometry of the $S^5 \to (S^5/\text{Euler})$ quotient), let us take a **different approach** that is more directly useful for implementation.

---

## 6. The TC Operator as a Modified Charge Function

### 6.1 Key Reformulation

Instead of working with the abstract $\partial/\partial r_{12}$ operator, we note that the TC Hamiltonian can be written as:

$$H_{\text{TC}} = T + V_{\text{nuc}} + G + \frac{1}{4}$$

where $V_{ee} = 1/r_{12}$ has been cancelled and replaced by $G$. The operator $G$ is first-order in derivatives. In the adiabatic framework where $\Psi = R^{-5/2}F(R)\Phi(\Omega)$, the angular eigenvalue problem at fixed $R$ becomes:

$$\left[\frac{\Lambda^2}{2} + R\,C_{\text{nuc}}(\alpha) + R\,G_{\text{ang}}(\alpha, \theta_{12})\right]\Phi = \mu_{\text{TC}}(R)\,\Phi$$

where $C_{\text{nuc}} = -Z(1/\cos\alpha + 1/\sin\alpha)$ is the nuclear charge function (unchanged), and $G_{\text{ang}}$ is the angular part of $G$ (the $V_{ee}$ part of the charge function has been replaced).

**The V_ee Neumann expansion $1/\sqrt{1 - \sin 2\alpha\cos\theta_{12}}$ is completely removed.** This is the term that couples different $l$-channels and creates the convergence bottleneck.

### 6.2 Structure of the Gradient Operator G

The operator $G = -\frac{1}{2}(\hat{r}_{12}\cdot\nabla_1 - \hat{r}_{12}\cdot\nabla_2)$ in the $L=0$ sector can be decomposed into:

**(a) A hyperradial component** $G_R$ that acts on $F(R)$:
$$G_R = -\frac{r_{12}}{2R}\frac{\partial}{\partial R} + \ldots$$

which modifies the non-adiabatic coupling between adiabatic channels but does not affect the angular eigenvalue problem at fixed $R$.

**(b) An angular component** $G_{\Omega}$ that acts within the angular eigenvalue problem at fixed $R$. This is the piece we need.

### 6.3 Angular Component of G

To extract the angular part, consider: at fixed $R$, $r_{12} = R\sqrt{u(\alpha, \theta_{12})}$. The operator $G$ projected onto the angular subspace at fixed $R$ involves derivatives with respect to $\alpha$ and $\theta_{12}$.

Using the relative/CM decomposition and projecting onto $L=0$, the angular part of $G$ can be written as:

$$G_{\Omega} = -\frac{1}{2R}\left[A(\alpha, \theta_{12})\frac{\partial}{\partial\alpha} + B(\alpha, \theta_{12})\frac{1}{\sin\alpha\cos\alpha}\frac{\partial}{\partial\theta_{12}}\right]$$

where $A$ and $B$ encode the projection of $\hat{r}_{12}$ onto the $\alpha$ and $\theta_{12}$ directions with the appropriate metric weights.

To derive $A$ and $B$, we use the explicit form of $\nabla_i$ in hyperspherical coordinates. For electron $i$ with position $\mathbf{r}_i = r_i\hat{r}_i$:

$$\nabla_i = \hat{r}_i\frac{\partial}{\partial r_i} + \frac{1}{r_i}\nabla_{\hat{r}_i}$$

For $L=0$ (integrating out the angular momentum directions), $\hat{r}_{12}\cdot\nabla_{\hat{r}_i}$ contributes through the $\theta_{12}$ derivative.

The interelectron unit vector is:
$$\hat{r}_{12} = \frac{\mathbf{r}_1 - \mathbf{r}_2}{r_{12}} = \frac{r_1\hat{r}_1 - r_2\hat{r}_2}{r_{12}} = \frac{R\cos\alpha\,\hat{r}_1 - R\sin\alpha\,\hat{r}_2}{R\sqrt{u}}$$

Then:
$$\hat{r}_{12}\cdot\nabla_1 = \frac{1}{\sqrt{u}}\left[\cos\alpha\frac{\partial}{\partial r_1} + \frac{(-\sin\alpha\cos\theta_{12})}{r_1}\frac{\partial}{\partial\theta_{12}} + \ldots\right]$$

where we used $\hat{r}_{12}\cdot\hat{r}_1 = (r_1 - r_2\cos\theta_{12})/r_{12} = (\cos\alpha - \sin\alpha\cos\theta_{12})/\sqrt{u}$ and $\hat{r}_{12}\cdot\hat{\theta}_{12,1} = (r_2\sin\theta_{12})/(r_{12}\cdot r_1 \cdot ...)$ (the angular component).

This is getting complicated. Let me take the most direct route.

### 6.4 Direct Computation via $\partial r_{12}/\partial$ coordinates

We know the physical operator is $G = -\partial/\partial r_{12}$. We want its matrix elements in the angular basis at fixed $R$. Using the Jacobian:

$$\frac{\partial}{\partial r_{12}}\bigg|_{\text{other indep. coords}} = \sum_q \frac{\partial q}{\partial r_{12}}\frac{\partial}{\partial q}$$

The key question is: at fixed $R$ and fixed center-of-mass position $\mathbf{S}$, what is $\partial/\partial r_{12}$?

For two electrons with $r_1 = R\cos\alpha$, $r_2 = R\sin\alpha$, and the angle between them $\theta_{12}$, the center-of-mass is $\mathbf{S} = \frac{1}{2}(r_1\hat{r}_1 + r_2\hat{r}_2)$, which depends on all the variables. Fixing $\mathbf{S}$ and $r_{12}$ simultaneously with fixing $R$ is overconstrained.

**The correct approach for the angular eigenvalue problem** is to note that within the adiabatic separation $\Psi = R^{-5/2}F(R)\Phi(\Omega)$, the angular eigenproblem is:

$$\left[\frac{\Lambda^2}{2} + R\,C(\alpha, \theta_{12})\right]\Phi = \mu(R)\Phi$$

In the TC version, the charge function $C$ is modified. The V_ee contribution to $C$ is:

$$C_{ee} = \frac{1}{\sqrt{1 - \sin 2\alpha\cos\theta_{12}}} = \frac{1}{\sqrt{u}}$$

This is replaced by the angular part of $G$. At fixed $R$, scaling arguments show:

$$G\Big|_{\text{fixed }R} = \frac{1}{R}\,G_{\text{ang}}(\alpha, \theta_{12})$$

where $G_{\text{ang}}$ is a first-order differential operator in $(\alpha, \theta_{12})$ that scales as $1/R$ (like all terms in the charge function). This is because $G$ involves one derivative and $r_{12} \sim R$, so $\partial/\partial r_{12} \sim 1/R$.

**The TC angular Hamiltonian is:**

$$H_{\text{ang}}^{\text{TC}} = \frac{\Lambda^2}{2} + R\,C_{\text{nuc}}(\alpha) + G_{\text{ang}}(\alpha, \theta_{12};\partial_\alpha, \partial_{\theta_{12}})$$

where $G_{\text{ang}}$ is a first-order operator with **smooth** coefficients (no $1/\sqrt{u}$ singularity).

---

## 7. Explicit Form of $G_{\text{ang}}$

### 7.1 From Cartesian to Hyperspherical

We derive $G_{\text{ang}}$ by expressing $\hat{r}_{12}\cdot\nabla_1$ and $\hat{r}_{12}\cdot\nabla_2$ in terms of $\partial/\partial\alpha$ and $\partial/\partial\theta_{12}$ at fixed $R$.

**Step 1: Radial derivatives in terms of $(R, \alpha)$.**

From $r_1 = R\cos\alpha$, $r_2 = R\sin\alpha$:

$$\frac{\partial}{\partial r_1}\bigg|_{r_2} = \frac{\cos\alpha}{1}\frac{\partial}{\partial R}\bigg|_\alpha + \frac{-\sin\alpha}{R}\frac{\partial}{\partial\alpha}\bigg|_R + \ldots$$

More carefully, the Jacobian of $(r_1, r_2) \to (R, \alpha)$ is:
$$r_1 = R\cos\alpha \implies dr_1 = \cos\alpha\,dR - R\sin\alpha\,d\alpha$$
$$r_2 = R\sin\alpha \implies dr_2 = \sin\alpha\,dR + R\cos\alpha\,d\alpha$$

Inverting:
$$dR = \cos\alpha\,dr_1 + \sin\alpha\,dr_2$$
$$d\alpha = \frac{-\sin\alpha}{R}dr_1 + \frac{\cos\alpha}{R}dr_2$$

Therefore:
$$\frac{\partial}{\partial r_1}\bigg|_{r_2,\hat{r}_1,\hat{r}_2} = \cos\alpha\frac{\partial}{\partial R} - \frac{\sin\alpha}{R}\frac{\partial}{\partial\alpha}$$

$$\frac{\partial}{\partial r_2}\bigg|_{r_1,\hat{r}_1,\hat{r}_2} = \sin\alpha\frac{\partial}{\partial R} + \frac{\cos\alpha}{R}\frac{\partial}{\partial\alpha}$$

**Step 2: Components of $\hat{r}_{12}\cdot\nabla_i$.**

The gradient of electron $i$ is $\nabla_i = \hat{r}_i\frac{\partial}{\partial r_i} + \frac{1}{r_i}\nabla_{\Omega_i}$. For the $L=0$ sector, the angular gradient terms contribute through $\theta_{12}$ derivatives. Specifically:

$$\hat{r}_{12}\cdot\hat{r}_1 = \frac{r_1 - r_2\cos\theta_{12}}{r_{12}} = \frac{\cos\alpha - \sin\alpha\cos\theta_{12}}{\sqrt{u}}$$

$$\hat{r}_{12}\cdot\hat{r}_2 = \frac{r_2\cos\theta_{12} - r_1}{r_{12}} \cdot (-1) = \frac{-(sin\alpha - \cos\alpha\cos\theta_{12})}{\sqrt{u}}$$

Wait, let me be more careful:
$$\hat{r}_{12} = \frac{\mathbf{r}_1 - \mathbf{r}_2}{r_{12}}$$

$$\hat{r}_{12}\cdot\hat{r}_2 = \frac{(\mathbf{r}_1 - \mathbf{r}_2)\cdot\hat{r}_2}{r_{12}} = \frac{r_1\cos\theta_{12} - r_2}{r_{12}} = \frac{\cos\alpha\cos\theta_{12} - \sin\alpha}{\sqrt{u}}$$

For the angular derivative, in the $L=0$ sector, $\hat{r}_{12}$ has a component perpendicular to both $\hat{r}_1$ and $\hat{r}_2$ that vanishes upon $L=0$ averaging. The remaining angular contribution is through $\theta_{12}$:

$$(\hat{r}_{12})_{\perp,1} = \frac{-r_2\sin\theta_{12}}{r_{12}} \cdot \hat{\theta}_{12,1} = \frac{-\sin\alpha\sin\theta_{12}}{\sqrt{u}} \cdot \hat{\theta}_{12,1}$$

where $\hat{\theta}_{12,1}$ is the unit vector in the direction of increasing $\theta_{12}$ at electron 1's position. Then $\hat{r}_{12}\cdot\frac{1}{r_1}\nabla_{\Omega_1} = \frac{-\sin\alpha\sin\theta_{12}}{r_1\sqrt{u}}\frac{\partial}{\partial\theta_{12}} = \frac{-\sin\alpha\sin\theta_{12}}{R\cos\alpha\sqrt{u}}\frac{\partial}{\partial\theta_{12}}$.

Similarly for electron 2:
$$(\hat{r}_{12})_{\perp,2}\cdot\frac{1}{r_2}\nabla_{\Omega_2} = \frac{+\cos\alpha\sin\theta_{12}}{R\sin\alpha\sqrt{u}}\frac{\partial}{\partial\theta_{12}}$$

(The sign is positive because $\hat{r}_{12}\cdot(-\hat{\theta}_{12,2}) = +\cos\alpha\sin\theta_{12}/\sqrt{u}$, and $\nabla_2$ picks up the minus from $\nabla_2 r_{12} = -\hat{r}_{12}$... let me be very careful.)

**Step 3: Assembling $G$.**

$$G = -\frac{1}{2}\hat{r}_{12}\cdot\nabla_1 + \frac{1}{2}\hat{r}_{12}\cdot\nabla_2$$

Radial contributions:
$$G_{\text{rad}} = -\frac{1}{2}\left[(\hat{r}_{12}\cdot\hat{r}_1)\frac{\partial}{\partial r_1} - (\hat{r}_{12}\cdot\hat{r}_2)\frac{\partial}{\partial r_2}\right]$$

$$= -\frac{1}{2\sqrt{u}}\left[(\cos\alpha - \sin\alpha\cos\theta_{12})\left(\cos\alpha\frac{\partial}{\partial R} - \frac{\sin\alpha}{R}\frac{\partial}{\partial\alpha}\right) - (\cos\alpha\cos\theta_{12} - \sin\alpha)\left(\sin\alpha\frac{\partial}{\partial R} + \frac{\cos\alpha}{R}\frac{\partial}{\partial\alpha}\right)\right]$$

The $\partial/\partial R$ coefficient:
$$\cos\alpha(\cos\alpha - \sin\alpha\cos\theta_{12}) - \sin\alpha(\cos\alpha\cos\theta_{12} - \sin\alpha)$$
$$= \cos^2\alpha - \sin\alpha\cos\alpha\cos\theta_{12} - \sin\alpha\cos\alpha\cos\theta_{12} + \sin^2\alpha$$
$$= 1 - \sin 2\alpha\cos\theta_{12} = u$$

So the $\partial/\partial R$ part gives $-\frac{u}{2\sqrt{u}}\frac{\partial}{\partial R} = -\frac{\sqrt{u}}{2}\frac{\partial}{\partial R}$, which is the hyperradial component (does not enter the angular eigenproblem at fixed $R$).

The $\partial/\partial\alpha$ coefficient:
$$\frac{1}{R}\left[\sin\alpha(\cos\alpha - \sin\alpha\cos\theta_{12}) + \cos\alpha(\cos\alpha\cos\theta_{12} - \sin\alpha)\right]$$
$$= \frac{1}{R}\left[\sin\alpha\cos\alpha - \sin^2\alpha\cos\theta_{12} + \cos^2\alpha\cos\theta_{12} - \sin\alpha\cos\alpha\right]$$
$$= \frac{1}{R}\left[\cos\theta_{12}(\cos^2\alpha - \sin^2\alpha)\right] = \frac{\cos 2\alpha\cos\theta_{12}}{R}$$

So the $\partial/\partial\alpha$ part gives $-\frac{1}{2\sqrt{u}}\cdot\frac{\cos 2\alpha\cos\theta_{12}}{R}\frac{\partial}{\partial\alpha}$.

Angular ($\theta_{12}$) contributions:
$$G_{\theta} = -\frac{1}{2}\left[\frac{-\sin\alpha\sin\theta_{12}}{R\cos\alpha\sqrt{u}}\right]\frac{\partial}{\partial\theta_{12}} + \frac{1}{2}\left[\frac{\cos\alpha\cos\theta_{12} - \sin\alpha}{R\sin\alpha\sqrt{u}} \cdot \frac{(-\sin\theta_{12})}{?}\right]$$

Let me redo this more carefully. For electron 1:
$$\hat{r}_{12}\cdot\nabla_1 = (\hat{r}_{12}\cdot\hat{r}_1)\frac{\partial}{\partial r_1} + \frac{1}{r_1}(\hat{r}_{12}\cdot\hat{e}_{\theta_{12},1})\frac{\partial}{\partial\theta_{12}}$$

where $\hat{e}_{\theta_{12},1}$ is the direction at $\hat{r}_1$ toward $\hat{r}_2$ (the direction of increasing $\theta_{12}$). In the plane of $\hat{r}_1$ and $\hat{r}_2$:

$$\hat{r}_{12} = \frac{R\cos\alpha\hat{r}_1 - R\sin\alpha\hat{r}_2}{R\sqrt{u}}$$

Decomposing $\hat{r}_2$ from the perspective of electron 1: $\hat{r}_2 = \cos\theta_{12}\hat{r}_1 + \sin\theta_{12}\hat{e}_{\theta_{12},1}$

$$\hat{r}_{12} = \frac{\cos\alpha(\hat{r}_1) - \sin\alpha(\cos\theta_{12}\hat{r}_1 + \sin\theta_{12}\hat{e}_{\theta_{12},1})}{\sqrt{u}}$$

$$= \frac{(\cos\alpha - \sin\alpha\cos\theta_{12})\hat{r}_1 - \sin\alpha\sin\theta_{12}\hat{e}_{\theta_{12},1}}{\sqrt{u}}$$

So:
$$\hat{r}_{12}\cdot\hat{e}_{\theta_{12},1} = \frac{-\sin\alpha\sin\theta_{12}}{\sqrt{u}}$$

And for electron 2, decomposing $\hat{r}_1 = \cos\theta_{12}\hat{r}_2 + \sin\theta_{12}\hat{e}_{\theta_{12},2}$ where $\hat{e}_{\theta_{12},2}$ points from $\hat{r}_2$ toward $\hat{r}_1$:

$$\hat{r}_{12} = \frac{\cos\alpha(\cos\theta_{12}\hat{r}_2 + \sin\theta_{12}\hat{e}_{\theta_{12},2}) - \sin\alpha\hat{r}_2}{\sqrt{u}}$$

$$= \frac{(\cos\alpha\cos\theta_{12} - \sin\alpha)\hat{r}_2 + \cos\alpha\sin\theta_{12}\hat{e}_{\theta_{12},2}}{\sqrt{u}}$$

So:
$$\hat{r}_{12}\cdot\hat{e}_{\theta_{12},2} = \frac{\cos\alpha\sin\theta_{12}}{\sqrt{u}}$$

Note: The $\theta_{12}$ derivative from electron 2's angular gradient has the opposite sign convention: $\frac{1}{r_2}\hat{e}_{\theta_{12},2}\frac{\partial}{\partial\theta_{12}}$ should be $\frac{1}{r_2}\hat{e}_{\theta_{12},2}\cdot(-\frac{\partial}{\partial\theta_{12}})$ because increasing $\theta_{12}$ from the perspective of electron 2 means $\hat{r}_2$ moving TOWARD $\hat{r}_1$, but the shared angle $\theta_{12}$ increases in the opposite direction. Actually, $\theta_{12}$ is the angle between $\hat{r}_1$ and $\hat{r}_2$, and both electrons see the same $\theta_{12}$, but the angular derivative from each electron's perspective has opposite sign in the direction of decreasing $\theta_{12}$.

For $L=0$, the $\theta_{12}$ derivative from electron 1's angular momentum:
$$\hat{e}_{\theta_{12},1}\cdot\frac{1}{r_1}\nabla_{\Omega_1}\,f(\theta_{12}) = \frac{1}{r_1}\frac{\partial f}{\partial\theta_{12}}$$

and from electron 2 (where moving $\hat{r}_2$ toward $\hat{r}_1$ DECREASES $\theta_{12}$):
$$\hat{e}_{\theta_{12},2}\cdot\frac{1}{r_2}\nabla_{\Omega_2}\,f(\theta_{12}) = -\frac{1}{r_2}\frac{\partial f}{\partial\theta_{12}}$$

Now assembling:
$$\hat{r}_{12}\cdot\nabla_1\Big|_{\theta_{12}\text{ part}} = \frac{-\sin\alpha\sin\theta_{12}}{\sqrt{u}}\cdot\frac{1}{R\cos\alpha}\frac{\partial}{\partial\theta_{12}}$$

$$\hat{r}_{12}\cdot\nabla_2\Big|_{\theta_{12}\text{ part}} = \frac{\cos\alpha\sin\theta_{12}}{\sqrt{u}}\cdot\frac{-1}{R\sin\alpha}\frac{\partial}{\partial\theta_{12}}$$

Therefore the $\theta_{12}$ part of $G$:

$$G_{\theta_{12}} = -\frac{1}{2}\left[\frac{-\sin\alpha\sin\theta_{12}}{R\cos\alpha\sqrt{u}}\right]\frac{\partial}{\partial\theta_{12}} + \frac{1}{2}\left[\frac{-\cos\alpha\sin\theta_{12}}{R\sin\alpha\sqrt{u}}\right]\frac{\partial}{\partial\theta_{12}}$$

$$= \frac{\sin\theta_{12}}{2R\sqrt{u}}\left[\frac{\sin\alpha}{\cos\alpha} - \frac{\cos\alpha}{\sin\alpha}\right]\frac{\partial}{\partial\theta_{12}}$$

$$= \frac{\sin\theta_{12}}{2R\sqrt{u}}\cdot\frac{\sin^2\alpha - \cos^2\alpha}{\sin\alpha\cos\alpha}\frac{\partial}{\partial\theta_{12}}$$

$$= \frac{-\cos 2\alpha\sin\theta_{12}}{R\sin 2\alpha\sqrt{u}}\frac{\partial}{\partial\theta_{12}}$$

### 7.2 Complete Angular TC Operator

Combining the $\alpha$ and $\theta_{12}$ parts (dropping the hyperradial $\partial/\partial R$ piece):

$$G_{\text{ang}} = -\frac{\cos 2\alpha\cos\theta_{12}}{2R\sqrt{u}}\frac{\partial}{\partial\alpha} - \frac{\cos 2\alpha\sin\theta_{12}}{R\sin 2\alpha\sqrt{u}}\frac{\partial}{\partial\theta_{12}}$$

where $u = 1 - \sin 2\alpha\cos\theta_{12}$.

**Singularity analysis:** Both terms have a $1/\sqrt{u}$ factor, which diverges at the coalescence point $(\alpha = \pi/4, \theta_{12} = 0)$. However, at this point:
- $\cos 2\alpha|_{\alpha=\pi/4} = 0$
- $\sin\theta_{12}|_{\theta_{12}=0} = 0$

So both numerators vanish at the coalescence point! Let us check the limiting behavior.

Near coalescence, set $\alpha = \pi/4 + \delta\alpha$, $\theta_{12} = \delta\theta$:
$$u \approx 2\delta\alpha^2\cdot 1 + 1\cdot\delta\theta^2/2 + \ldots = 2(\delta\alpha)^2 + (\delta\theta)^2/2$$

(Using $\sin 2\alpha \approx 1 - 2(\delta\alpha)^2$ and $\cos\theta_{12} \approx 1 - (\delta\theta)^2/2$.)

So $\sqrt{u} \approx \sqrt{2(\delta\alpha)^2 + (\delta\theta)^2/2}$.

For the $\alpha$-term: $\frac{\cos 2\alpha}{\sqrt{u}} \approx \frac{-2\delta\alpha}{\sqrt{2(\delta\alpha)^2 + (\delta\theta)^2/2}}$

This is bounded: as $\delta\alpha \to 0$ with $\delta\theta$ fixed, the ratio $\to 0$; along $\delta\theta = 0$, the ratio $\to -\sqrt{2}$ (finite). **The $\alpha$-term coefficient is bounded.**

For the $\theta_{12}$-term: $\frac{\cos 2\alpha\sin\theta_{12}}{\sin 2\alpha\sqrt{u}} \approx \frac{(-2\delta\alpha)(\delta\theta)}{1\cdot\sqrt{2(\delta\alpha)^2 + (\delta\theta)^2/2}}$

Again bounded. Along $\delta\alpha = 0$: ratio $\to 0$. Along $\delta\theta = 0$: ratio $\to 0$. The maximum is at $\delta\alpha \sim \delta\theta$ where the ratio is $O(1)$. **The $\theta_{12}$-term coefficient is also bounded.**

**CRITICAL RESULT: The TC angular operator $G_{\text{ang}}$ has NO singularities.** The $1/\sqrt{u}$ divergence from $\nabla r_{12}$ is exactly cancelled by the zeros of $\cos 2\alpha$ and $\sin\theta_{12}$ at the coalescence point. This is the direct consequence of the Kato cusp cancellation: the singular parts of the TC transformation ($\nabla^2 J$ and $V_{ee}$) cancel each other, leaving only smooth operators.

---

## 8. TC Angular Hamiltonian

### 8.1 Complete Form

The TC angular eigenvalue problem at fixed $R$ is:

$$\left[\frac{\Lambda^2}{2} + R\,C_{\text{nuc}}(\alpha)\right]\Phi + G_{\text{ang}}\Phi = \mu_{\text{TC}}(R)\,\Phi$$

where:

$$C_{\text{nuc}}(\alpha) = -Z\left(\frac{1}{\cos\alpha} + \frac{1}{\sin\alpha}\right)$$

$$G_{\text{ang}} = -\frac{\cos 2\alpha\cos\theta_{12}}{2R\sqrt{u}}\frac{\partial}{\partial\alpha} - \frac{\cos 2\alpha\sin\theta_{12}}{R\sin 2\alpha\sqrt{u}}\frac{\partial}{\partial\theta_{12}}$$

Note: $G_{\text{ang}}$ scales as $1/R$ (like the charge function), so the angular problem has the standard structure $[\Lambda^2/2 + R\cdot(\text{potential})] + (1/R)\cdot(\text{derivative}) = \mu$.

Wait — checking dimensions: $G_{\text{ang}}$ has a $1/R$ prefactor and involves $\partial/\partial\alpha$ and $\partial/\partial\theta_{12}$, which are dimensionless (angular derivatives). So $G_{\text{ang}} \sim 1/R$, which scales the SAME as the charge function terms. Correct.

In the standard normalization where the angular equation is:
$$\left[\frac{\Lambda^2}{2} + R\,C(\alpha, \theta_{12})\right]\Phi = \mu(R)\Phi$$

the TC version has:
$$\left[\frac{\Lambda^2}{2} + R\,C_{\text{nuc}}(\alpha) + \tilde{G}_{\text{ang}}\right]\Phi = \mu_{\text{TC}}(R)\Phi$$

where $\tilde{G}_{\text{ang}} = R\cdot G_{\text{ang}}$ absorbs the $R$ dependence:

$$\tilde{G}_{\text{ang}} = -\frac{\cos 2\alpha\cos\theta_{12}}{2\sqrt{u}}\frac{\partial}{\partial\alpha} - \frac{\cos 2\alpha\sin\theta_{12}}{\sin 2\alpha\sqrt{u}}\frac{\partial}{\partial\theta_{12}}$$

This operator is **R-independent** (like the original $V_{ee}$ coupling) and has smooth, bounded coefficients.

### 8.2 The $V_{ee}$ Replacement

The original angular Hamiltonian has:
$$H_{\text{ang}} = \frac{\Lambda^2}{2} + R\,C_{\text{nuc}}(\alpha) + R\cdot\frac{1}{\sqrt{u}} \cdot P_l(\text{multipole expansion})$$

where the last term is the $V_{ee}$ contribution (singular at $u = 0$).

The TC angular Hamiltonian replaces this with:
$$H_{\text{ang}}^{\text{TC}} = \frac{\Lambda^2}{2} + R\,C_{\text{nuc}}(\alpha) + \tilde{G}_{\text{ang}}$$

where $\tilde{G}_{\text{ang}}$ has **no singularity**. The $1/\sqrt{u}$ from $V_{ee}$ is gone.

**This is a dramatic simplification.** The Neumann expansion of $1/\sqrt{u}$ that creates the V_ee coupling between all $l$-channels (via Gaunt integrals, the $\sum_k (r_</r_>)^k P_k$ series) is completely eliminated.

---

## 9. Selection Rules for TC Matrix Elements

### 9.1 Expansion of the TC Coefficients

The TC operator involves two coefficient functions:

$$A(\alpha, \theta_{12}) = \frac{\cos 2\alpha\cos\theta_{12}}{2\sqrt{u}} \qquad (\text{multiplies } \partial/\partial\alpha)$$

$$B(\alpha, \theta_{12}) = \frac{\cos 2\alpha\sin\theta_{12}}{\sin 2\alpha\sqrt{u}} \qquad (\text{multiplies } \partial/\partial\theta_{12})$$

To compute matrix elements in the basis $\phi_{\nu,l} \propto \chi_l(\alpha) P_l(\cos\theta_{12})$, we need the Legendre expansion of these coefficients in $\theta_{12}$.

**For $A$:** The $\theta_{12}$-dependence is $\cos\theta_{12}/\sqrt{1 - \sin 2\alpha\cos\theta_{12}}$. Using the Neumann-type expansion:

$$\frac{\cos\theta_{12}}{\sqrt{1 - t\cos\theta_{12}}} = \sum_{l=0}^{\infty} a_l(t) P_l(\cos\theta_{12})$$

where $t = \sin 2\alpha$. This can be derived from the standard generating function $(1 - 2tx + t^2)^{-1/2} = \sum_l t^l P_l(x)$ with appropriate variable substitution. However, the Neumann expansion of $1/\sqrt{1 - t\cos\theta_{12}}$ is NOT the standard Legendre generating function (which has $-2tx + t^2$ under the square root, not $-t\cos\theta$).

The expansion of $1/\sqrt{1 - t\cos\theta}$ for $|t| < 1$ is related to the Gegenbauer expansion:

$$\frac{1}{\sqrt{1 - t\cos\theta}} = \sum_{k=0}^{\infty} t^k \frac{(2k)!}{4^k(k!)^2} P_k(\cos\theta) \cdot \ldots$$

Actually, from the standard multipole expansion used in the code:
$$\frac{1}{\sqrt{1 - \sin 2\alpha\cos\theta_{12}}} = \sum_k \frac{(\min(\sin\alpha, \cos\alpha))^k}{(\max(\sin\alpha, \cos\alpha))^{k+1}} P_k(\cos\theta_{12})$$

This is the Neumann expansion where $r_< = R\min(\sin\alpha, \cos\alpha)$ and $r_> = R\max(\sin\alpha, \cos\alpha)$.

Multiplying by $\cos\theta_{12}$: since $\cos\theta_{12} \cdot P_l = \frac{l+1}{2l+1}P_{l+1} + \frac{l}{2l+1}P_{l-1}$, the $\cos\theta_{12}$ multiplication couples $l \to l\pm 1$. So the $\alpha$-derivative term couples channels with $\Delta l = \pm 1$ (from the additional $\cos\theta_{12}$ factor) **on top of** whatever coupling the Neumann expansion produces.

**For $B$:** The $\theta_{12}$-dependence is $\sin\theta_{12}/\sqrt{u}$, and $\partial/\partial\theta_{12}$ acting on $P_l(\cos\theta_{12})$ gives $-\sin\theta_{12} P_l'(\cos\theta_{12})$. So the combined integrand involves:

$$B \cdot \frac{\partial P_{l'}}{\partial\theta_{12}} = \frac{\cos 2\alpha}{\sin 2\alpha\sqrt{u}}\sin\theta_{12}\cdot(-\sin\theta_{12})P_{l'}'(\cos\theta_{12})$$

$$= -\frac{\cos 2\alpha}{\sin 2\alpha\sqrt{u}}\sin^2\theta_{12}\, P_{l'}'(\cos\theta_{12})$$

Using $\sin^2\theta P_l'(\cos\theta) = l P_{l-1}(\cos\theta) - l\cos\theta P_l(\cos\theta)$ (from the recurrence), this couples channels via $\Delta l = 0, \pm 1$.

### 9.2 Selection Rules Summary

The TC gradient operator $\tilde{G}_{\text{ang}}$ has the following selection rules:

| Term | $\Delta l$ | $\Delta k$ | Note |
|:-----|:----------:|:----------:|:-----|
| $A\cdot\partial_\alpha$ | $\pm 1$ (from $\cos\theta_{12}$) + $0, \pm 2, \ldots$ (from $1/\sqrt{u}$ expansion) | All (from $\partial_\alpha$ on Gegenbauer) | Couples adjacent $l$-channels |
| $B\cdot\partial_{\theta_{12}}$ | $0, \pm 1$ (from $\sin^2\theta P'_l$) + $0, \pm 2, \ldots$ (from $1/\sqrt{u}$) | All (from $\alpha$-dependent coefficient) | Couples nearby $l$-channels |

The **key difference from $V_{ee}$**: the original $V_{ee}$ coupling has selection rule $\Delta l = 0, \pm 2, \pm 4, \ldots$ (even, from the Gaunt integral). The TC gradient operator introduces **odd** $\Delta l$ coupling ($\Delta l = \pm 1$) because of the $\cos\theta_{12}$ factor.

However, the Neumann expansion $1/\sqrt{u}$ in the TC coefficients ALSO couples all $l$-channels, similar to the original $V_{ee}$. The difference is that the **integrand is smooth** (no $1/\sqrt{u}$ left uncompensated), so the matrix elements converge exponentially in $l$.

### 9.3 Non-Hermiticity

The TC angular Hamiltonian $H_{\text{ang}}^{\text{TC}}$ is **non-Hermitian** because $\tilde{G}_{\text{ang}}$ is a first-order differential operator:

$$\langle\phi_i|\tilde{G}|\phi_j\rangle \neq \langle\phi_j|\tilde{G}|\phi_i\rangle^*$$

In fact, by integration by parts:
$$\langle\phi_i|\tilde{G}|\phi_j\rangle + \langle\phi_j|\tilde{G}|\phi_i\rangle^* = \text{boundary terms} + \text{divergence terms}$$

The antisymmetric (non-Hermitian) part comes from the first-derivative operator, which is anti-Hermitian up to boundary terms. This requires a non-symmetric eigensolver (e.g., `numpy.linalg.eig` instead of `scipy.linalg.eigh`).

---

## 10. Matrix Element Formulas

### 10.1 Basis Functions

The angular basis for the TC problem is the same Gegenbauer $\times$ Legendre basis:
$$\phi_{k,l}(\alpha, \theta_{12}) = N_{k,l}\,(\sin\alpha\cos\alpha)^{l+1}\,C_k^{l+1}(\cos 2\alpha)\, P_l(\cos\theta_{12})$$

with $k = 0, 2, 4, \ldots$ (singlet) and normalization $N_{k,l}$ fixed by $\langle\phi_{k,l}|\phi_{k,l}\rangle = 1$ in the $S^5$ measure:

$$d\Omega_5 = \sin^2\alpha\cos^2\alpha\sin\theta_{12}\,d\alpha\,d\theta_{12}\,d(\text{Euler})$$

For $L=0$, after integrating out Euler angles, the effective measure is:
$$d\mu = \sin^2\alpha\cos^2\alpha\sin\theta_{12}\,d\alpha\,d\theta_{12}$$

### 10.2 Matrix Elements of $\tilde{G}_{\text{ang}}$

The matrix element $\langle\phi_{k,l}|\tilde{G}_{\text{ang}}|\phi_{k',l'}\rangle$ splits into two terms:

**$\alpha$-derivative term:**
$$T_\alpha = -\int_0^{\pi/4}d\alpha\int_0^{\pi}d\theta_{12}\,\sin^2\alpha\cos^2\alpha\sin\theta_{12}\,\phi_{k,l}^*\,\frac{\cos 2\alpha\cos\theta_{12}}{2\sqrt{u}}\,\frac{\partial\phi_{k',l'}}{\partial\alpha}$$

**$\theta_{12}$-derivative term:**
$$T_\theta = -\int_0^{\pi/4}d\alpha\int_0^{\pi}d\theta_{12}\,\sin^2\alpha\cos^2\alpha\sin\theta_{12}\,\phi_{k,l}^*\,\frac{\cos 2\alpha\sin\theta_{12}}{\sin 2\alpha\sqrt{u}}\,\frac{\partial\phi_{k',l'}}{\partial\theta_{12}}$$

Since $\frac{\partial\phi_{k',l'}}{\partial\theta_{12}} = N_{k',l'}\chi_{k',l'}(\alpha)\cdot(-\sin\theta_{12})P_{l'}'(\cos\theta_{12})$, the $\theta_{12}$ integral in $T_\theta$ involves:

$$\int_0^\pi \frac{\sin^2\theta_{12}\sin\theta_{12}}{\sqrt{u}}\,P_l(\cos\theta_{12})\,P_{l'}'(\cos\theta_{12})\,d\theta_{12}$$

where $u = 1 - t\cos\theta_{12}$ with $t = \sin 2\alpha$.

Substituting $x = \cos\theta_{12}$:

$$\int_{-1}^{1}\frac{(1-x^2)}{(1-tx)^{1/2}}\,P_l(x)\,P_{l'}'(x)\,dx$$

This integral can be evaluated analytically by expanding $(1-tx)^{-1/2}$ in Legendre polynomials and using the orthogonality relations. Similarly, the $\alpha$ integrals involve Gegenbauer polynomials with smooth weight functions.

### 10.3 Practical Evaluation Strategy

For numerical implementation, the matrix elements can be computed by Gauss-Legendre quadrature (as already used in `AlgebraicAngularSolver._precompute_coupling`):

1. Evaluate $\phi_{k,l}(\alpha_i, \theta_j)$, $\partial_\alpha\phi_{k',l'}(\alpha_i, \theta_j)$, and $\partial_{\theta}\phi_{k',l'}(\alpha_i, \theta_j)$ on a 2D quadrature grid.
2. Evaluate the smooth coefficient functions $A(\alpha_i, \theta_j)$ and $B(\alpha_i, \theta_j)$.
3. Accumulate the weighted sums.

The existing infrastructure in `algebraic_angular.py` already splits the $\alpha$ integration at $\pi/4$ and uses Gauss-Legendre quadrature. The $\theta_{12}$ integration can use the same approach (the Legendre polynomial orthogonality gives exact results for low-order integrands, and quadrature handles the $1/\sqrt{u}$ weight).

**Key advantage:** Since the coefficients $A$ and $B$ are smooth (Section 7.2 singularity analysis), the quadrature converges rapidly. This contrasts with the original $V_{ee}$ matrix elements where the $1/\sqrt{u}$ singularity requires careful handling.

---

## 11. The Full TC Energy Shift

For completeness, the eigenvalue $E$ of $H_{\text{TC}}$ relates to the physical energy as:

$$E_{\text{physical}} = E_{\text{TC}} - \frac{1}{4}$$

(subtracting the constant shift from Term 3).

---

## 12. Summary Table

| Term | Expression | Singularity Status | Selection Rules | Implementation Notes |
|:-----|:-----------|:-------------------|:----------------|:--------------------|
| **Original $V_{ee}$** | $+1/r_{12} = 1/(R\sqrt{u})$ | **Singular** ($1/\sqrt{u}$ at coalescence) | $\Delta l = 0, \pm 2, \pm 4, \ldots$ via Gaunt | Neumann multipole expansion, Gaunt integrals |
| **$-\frac{1}{2}\sum_i\nabla_i^2 J$** (Laplacian) | $-1/r_{12}$ | **Singular** (cancels $V_{ee}$ exactly) | Same as $V_{ee}$ | Cancels $V_{ee}$; net contribution = 0 |
| **$\frac{1}{2}[[T,J],J]$** (double commutator) | $+1/4$ (constant) | **None** | None (scalar shift) | Trivial: shift all eigenvalues by $+1/4$ |
| **$-\sum_i(\nabla_i J)\cdot\nabla_i$** (gradient) | $\tilde{G}_{\text{ang}} = -\frac{\cos 2\alpha\cos\theta_{12}}{2\sqrt{u}}\partial_\alpha - \frac{\cos 2\alpha\sin\theta_{12}}{\sin 2\alpha\sqrt{u}}\partial_{\theta_{12}}$ | **Smooth** (numerators vanish at coalescence, cancelling $1/\sqrt{u}$) | $\Delta l = \pm 1, \pm 3, \ldots$ (odd from $\cos\theta_{12}$) plus $\Delta l = 0, \pm 1, \ldots$ (from $\sin^2\theta P'_l$) | Non-Hermitian first-order operator; 2D quadrature over smooth integrand; uses existing Gauss-Legendre infrastructure |

## 13. Assessment of Implementation Difficulty

### Key Simplifications

1. **$V_{ee}$ completely cancels** — the most complex part of the angular coupling (Neumann expansion, Gaunt integrals) is eliminated.
2. **Double commutator is trivial** — just $+1/4$ energy shift.
3. **TC operator is smooth** — no singularities in the angular matrix elements, so standard quadrature converges rapidly.
4. **Same basis** — the existing Gegenbauer $\times$ Legendre basis is used unchanged.
5. **R-independent TC coupling** — the operator $\tilde{G}_{\text{ang}}$ does not depend on $R$, so the TC correction matrices are precomputed once (like the existing $V_{ee}$ coupling).

### Implementation Tasks

1. **Compute $\partial\phi_{k,l}/\partial\alpha$**: Requires derivative of Gegenbauer polynomials. The recurrence $\frac{d}{dx}C_k^\lambda(x) = 2\lambda C_{k-1}^{\lambda+1}(x)$ gives this analytically. Combined with the chain rule for the $(\sin\alpha\cos\alpha)^{l+1}$ envelope. Moderate effort.

2. **Compute $\partial\phi_{k,l}/\partial\theta_{12}$**: This is just $\phi_{k,l}\cdot(-\sin\theta_{12})P_l'(\cos\theta_{12})/P_l(\cos\theta_{12})$, but more precisely involves the Legendre derivative. Standard recurrence: $(1-x^2)P_l'(x) = -lxP_l + lP_{l-1}$. Low effort.

3. **Evaluate smooth coefficients $A(\alpha, \theta_{12})$ and $B(\alpha, \theta_{12})$ on quadrature grid**: Direct evaluation of the formulas. The $1/\sqrt{u}$ near coalescence is cancelled analytically; numerically, use L'Hopital or Taylor expansion near $u \approx 0$ to avoid 0/0. Low effort.

4. **Assemble TC correction matrix**: 2D quadrature (existing infrastructure) with smooth integrand. Low effort.

5. **Non-Hermitian eigensolver**: Replace `scipy.linalg.eigh` with `numpy.linalg.eig`. The eigenvalues should still be real (provably for exact Jastrow). Check for imaginary parts as a diagnostic. Low effort.

### Overall Difficulty: MODERATE

The main new work is computing the $\alpha$-derivatives of the Gegenbauer basis functions and setting up the 2D quadrature for the TC matrix elements. The theoretical derivation (this document) is the hard part; the implementation is straightforward given the existing infrastructure in `algebraic_angular.py`.

### Expected Outcome

At $l_{\max} = 2$, the standard solver gives 0.24% error (0.10% with Schwartz post-correction). The TC solver should converge **exponentially** in $l_{\max}$ because:

1. The $V_{ee}$ singularity is removed, so all matrix elements are smooth.
2. The TC eigenfunctions are smooth (cusp-free), so they are well-represented by the angular basis at low $l_{\max}$.
3. The convergence bottleneck (partial-wave expansion of $1/r_{12}$ cusp) is eliminated.

Target: TC at $l_{\max} = 2$ should achieve $< 0.05\%$ error, competitive with standard at $l_{\max} = 5$ + Schwartz correction.
