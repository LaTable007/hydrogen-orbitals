# Hydrogen Atom Orbitals

Interactive 3D visualization of hydrogen atom orbitals built with Rust and [three-d](https://github.com/asny/three-d).

https://github.com/LaTable007/hydrogen-orbitals/releases/download/v1.0/demo.mp4

## Features

- **60 000 points** sampled via rejection sampling from the exact quantum probability density
- **3 visualization modes** switchable at runtime:
  - **Probability density** |ψ|² — viridis colormap
  - **Wavefunction phase** arg(ψ) = mφ — full HSV colour wheel showing angular momentum winding
  - **Probability current** **j** — azimuthal current j_φ = m|ψ|²/(r sinθ), shown as elongated ellipsoids
- **Angular clipping** — two independent sliders (azimuthal φ and polar θ) to cut into the orbital and observe its interior structure
- **XYZ axes** (red/green/blue) scaled to each orbital
- All 18 orbitals from 1s to 4f with quantum numbers n, l, m

## Physics

The visualization uses the **complex eigenstates** ψ_{nlm} = R_{nl}(r) · Y_l^m(θ,φ):

- **Density**: |ψ|² = R²·N²·P_l^|m|(cosθ)² — axially symmetric for m ≠ 0
- **Phase**: arg(ψ) = m·φ + π·H(−R·N·P) — the phase winds m times around the Z axis
- **Current**: j_φ = m · |ψ|² / (r sinθ) — purely azimuthal, zero for m = 0

## Controls

| Action | Result |
|--------|--------|
| Click + drag | Rotate |
| Scroll | Zoom |
| φ slider | Azimuthal cut [0°, φ] |
| θ slider | Polar cut [0°, θ] |

## Build & Run

```bash
cargo run --release
```

Requires Rust 1.86+ and a machine with OpenGL 3.3 support.
