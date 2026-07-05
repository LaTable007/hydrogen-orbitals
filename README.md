# Hydrogen Atom Orbitals

Interactive 3D visualization of hydrogen atom orbitals built with Rust and [three-d](https://github.com/asny/three-d).

---

## Demos

### Probability density |ψ|²

> 60 000 points sampled by rejection sampling — viridis colormap.

https://github.com/LaTable007/hydrogen-orbitals/releases/download/v1.0/demo_density.mp4

---

### Wavefunction phase arg(ψ)

> Phase rotates at ω = 1/(2n²) a.u. — full HSV colour wheel. The phase winds m times around the Z axis each period.

https://github.com/LaTable007/hydrogen-orbitals/releases/download/v1.0/demo_phase.mp4

---

### Probability current **j**

> Particles physically orbit around Z. Inner particles orbit faster: v_φ = m / ρ_cyl (differential rotation). m > 0 → CCW, m < 0 → CW, m = 0 → static.

https://github.com/LaTable007/hydrogen-orbitals/releases/download/v1.0/demo_current.mp4

---

### Full overview

https://github.com/LaTable007/hydrogen-orbitals/releases/download/v1.0/demo.mp4

---

## Features

- **60 000 points** sampled via rejection sampling from the exact quantum probability density
- **18 orbitals** from 1s to 4f with quantum numbers n, l, m
- **3 visualization modes** switchable at runtime:
  - **Probability density** |ψ|² — viridis colormap
  - **Wavefunction phase** arg(ψ) = mφ — HSV colour wheel showing angular momentum winding
  - **Probability current** **j** — particles physically orbit around Z axis
- **Time animation** — Play/Pause, speed slider, period scales as n²
- **Angular clipping** — two independent sliders (azimuthal φ and polar θ) to cut into the orbital and observe its interior structure
- **XYZ axes** (red/green/blue) scaled to each orbital

## Physics

Uses the **complex eigenstates** ψ_{nlm} = R_{nl}(r) · Y_l^m(θ,φ):

| Quantity | Formula |
|----------|---------|
| Density | \|ψ\|² = R²·N²·P_l^\|m\|(cosθ)² — axisymmetric for m ≠ 0 |
| Phase | arg(ψ(t)) = mφ + t/(2n²) + π·H(−R·N·P) |
| Current | j_φ = m · \|ψ\|² / (r sinθ) — purely azimuthal |
| Velocity field | v_φ = j_φ / \|ψ\|² = m / ρ_cyl |

## Controls

| Action | Result |
|--------|--------|
| Click + drag | Rotate |
| Scroll | Zoom |
| φ slider | Azimuthal cut [0°, φ] |
| θ slider | Polar cut [0°, θ] |
| ▶ / ⏸ | Toggle time animation |
| Speed slider | 0.1× – 8× real-time |

## Build & Run

```bash
cargo run --release
```

Requires Rust 1.86+ and OpenGL 3.3.
