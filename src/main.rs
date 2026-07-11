use three_d::*;
use rand::Rng;
use std::f64::consts::PI;

// ── Maths ─────────────────────────────────────────────────────────────────────

fn factorial(n: u32) -> f64 {
    (1..=n).map(|x| x as f64).product::<f64>().max(1.0)
}

fn associated_legendre(l: i32, m_abs: u32, x: f64) -> f64 {
    let sin_t = (1.0 - x * x).max(0.0).sqrt();
    let mut pmm = 1.0_f64;
    for i in 1..=m_abs {
        pmm *= -((2 * i - 1) as f64) * sin_t;
    }
    if l == m_abs as i32 {
        return pmm;
    }
    let mut pmm1 = x * (2 * m_abs + 1) as f64 * pmm;
    if l == m_abs as i32 + 1 {
        return pmm1;
    }
    let mut pll = 0.0;
    for ll in (m_abs as i32 + 2)..=l {
        pll = ((2 * ll - 1) as f64 * x * pmm1
            - (ll as f64 + m_abs as f64 - 1.0) * pmm)
            / (ll as f64 - m_abs as f64);
        pmm = pmm1;
        pmm1 = pll;
    }
    pll
}

fn laguerre_associated(n: u32, l: u32, x: f64) -> f64 {
    let alpha = (2 * l + 1) as f64;
    let k = n - l - 1;
    if k == 0 {
        return 1.0;
    }
    let mut l0 = 1.0_f64;
    let mut l1 = 1.0 + alpha - x;
    for j in 2..=k {
        let tmp = (((2 * j - 1) as f64 + alpha) - x) * l1 - ((j - 1) as f64 + alpha) * l0;
        l0 = l1;
        l1 = tmp / j as f64;
    }
    l1
}

/// Eigenstate ψ_{nlm} with its normalisation constants precomputed once,
/// instead of re-deriving factorials for each of the millions of Monte-Carlo
/// samples.
///
/// |ψ|² = R² · N² · P²(cosθ)   — no φ dependence (axisymmetric for m≠0)
/// arg(ψ)  = m·φ  +  π·H(−R·N·P)   where H is the Heaviside function
/// j_φ  = m · |ψ|² / (r sinθ)   →  |j| = |m| · |ψ|² / ρ_cyl
#[derive(Clone, Copy)]
struct Orbital {
    n: u32,
    l: i32,
    m: i32,
    m_abs: u32,
    /// true → real linear combination (ψ_{l,m} ± ψ_{l,−m})/√2, i.e. the
    /// chemist's p_x, d_xy, … lobes, instead of the complex L_z eigenstate.
    real: bool,
    inv_n: f64,
    norm_r: f64, // radial normalisation
    norm_a: f64, // angular normalisation of Y_lm
}

impl Orbital {
    fn new(n: u32, l: i32, m: i32, real: bool) -> Self {
        let l_u = l as u32;
        let m_abs = m.unsigned_abs();
        // Convention moderne des polynômes de Laguerre (celle de la récurrence
        // dans laguerre_associated): le dénominateur est (n+l)!, pas ((n+l)!)³.
        // Vérif: R_10 = 2e^{-r}, R_21 = r·e^{-r/2}/√24.
        let norm_r_sq = (2.0 / n as f64).powi(3) * factorial(n - l_u - 1)
            / (2.0 * n as f64 * factorial(n + l_u));
        let norm_a = (((2 * l + 1) as f64 / (4.0 * PI))
            * (factorial(l_u - m_abs) / factorial(l_u + m_abs)))
            .sqrt();
        Self {
            n,
            l,
            m,
            m_abs,
            real,
            inv_n: 1.0 / n as f64,
            norm_r: norm_r_sq.sqrt(),
            norm_a,
        }
    }

    /// Normalised radial part R_{nl}(r).
    fn radial(&self, r: f64) -> f64 {
        let rho = 2.0 * r * self.inv_n;
        self.norm_r
            * (-rho / 2.0).exp()
            * rho.powi(self.l)
            * laguerre_associated(self.n, self.l as u32, rho)
    }

    /// Polar part Θ(θ) = N·P_lm(cosθ) — Y_lm without the e^{imφ} factor.
    fn angular(&self, cos_theta: f64) -> f64 {
        self.norm_a * associated_legendre(self.l, self.m_abs, cos_theta)
    }

    /// Azimuthal factor of the real form: √2·cos(|m|φ) for m>0, √2·sin(|m|φ) for m<0.
    fn azimuthal_real(&self, phi: f64) -> f64 {
        use std::f64::consts::SQRT_2;
        if self.m > 0 {
            SQRT_2 * (self.m_abs as f64 * phi).cos()
        } else {
            SQRT_2 * (self.m_abs as f64 * phi).sin()
        }
    }

    /// Signed real amplitude (for the complex eigenstate: without e^{imφ}).
    fn amplitude(&self, x: f64, y: f64, z: f64) -> f64 {
        let r = (x * x + y * y + z * z).sqrt();
        if r < 1e-10 {
            return 0.0;
        }
        let cos_theta = (z / r).clamp(-1.0, 1.0);
        let mut amp = self.radial(r) * self.angular(cos_theta);
        if self.real && self.m != 0 {
            amp *= self.azimuthal_real(y.atan2(x));
        }
        amp
    }

    /// Fast path for rejection sampling: density only, no phase/current.
    fn density(&self, x: f64, y: f64, z: f64) -> f64 {
        let amp = self.amplitude(x, y, z);
        amp * amp
    }

    fn values(&self, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
        let amp = self.amplitude(x, y, z);
        let density = amp * amp;

        // Phase:  arg(ψ) = m·φ + (π if amp < 0) — real orbitals are… real,
        // so their phase is just the sign of the lobe (0 or π).
        let phase = if self.real {
            if amp >= 0.0 { 0.0 } else { PI }
        } else {
            let phi = y.atan2(x); // [-π, π]
            let phase_base = if amp >= 0.0 { 0.0 } else { PI };
            ((phase_base + self.m as f64 * phi) % (2.0 * PI) + 2.0 * PI) % (2.0 * PI)
        };

        // Probability current: j_φ = m·|ψ|²/(r·sinθ) = m·density/ρ_cyl.
        // A real wavefunction carries no current.
        let j_mag = if self.real {
            0.0
        } else {
            let rho_cyl = (x * x + y * y).sqrt().max(1e-10);
            self.m_abs as f64 * density / rho_cyl
        };

        (density, phase, j_mag)
    }
}

// ── Colours ───────────────────────────────────────────────────────────────────

fn viridis(t: f32) -> Srgba {
    let t = t.clamp(0.0, 1.0);
    let (r, g, b): (f32, f32, f32) = if t < 0.25 {
        let s = t / 0.25;
        (0.267 + s * -0.036, 0.005 + s * 0.317, 0.329 + s * 0.217)
    } else if t < 0.5 {
        let s = (t - 0.25) / 0.25;
        (0.231 + s * -0.102, 0.322 + s * 0.247, 0.546 + s * 0.003)
    } else if t < 0.75 {
        let s = (t - 0.5) / 0.25;
        (0.129 + s * 0.240, 0.569 + s * 0.219, 0.549 + s * -0.165)
    } else {
        let s = (t - 0.75) / 0.25;
        (0.369 + s * 0.624, 0.788 + s * 0.118, 0.384 + s * -0.239)
    };
    Srgba::new((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8, 220)
}

fn phase_color(phase: f32) -> Srgba {
    // Full HSV colour wheel: hue = phase / 2π
    let h = (phase / (2.0 * PI as f32) * 360.0) % 360.0;
    let c = 1.0_f32;
    let x = c * (1.0 - ((h / 60.0) % 2.0 - 1.0).abs());
    let (r, g, b): (f32, f32, f32) = if h < 60.0 {
        (c, x, 0.0)
    } else if h < 120.0 {
        (x, c, 0.0)
    } else if h < 180.0 {
        (0.0, c, x)
    } else if h < 240.0 {
        (0.0, x, c)
    } else if h < 300.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };
    Srgba::new((r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8, 230)
}

// ── Data ──────────────────────────────────────────────────────────────────────

#[derive(Clone, Copy, PartialEq)]
enum VisuMode {
    Density,
    WavePhase,
    ProbCurrent,
}

struct RawPoint {
    pos: Vec3,
    phi: f32,     // [0, 2π] — azimuthal angle for clipping & animation
    theta: f32,   // [0, π]  — polar angle for clipping
    rho_cyl: f32, // √(x²+y²) — cylindrical radius, drives per-particle orbital speed
    density: f32, // normalised [0,1] with γ-correction
    phase: f32,   // [0, 2π]
    j_norm: f32,  // |j| normalised [0,1]
}

fn generate_raw_points(orb: Orbital, num_points: usize) -> Vec<RawPoint> {
    let scale = (orb.n as f64).powi(2) * 3.5;
    let threads = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4);

    // Estimate max density for rejection sampling (parallel)
    let max_density = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut md = 1e-40_f64;
                    for _ in 0..(100_000 / threads + 1) {
                        let x = rng.gen_range(-scale..scale);
                        let y = rng.gen_range(-scale..scale);
                        let z = rng.gen_range(-scale..scale);
                        md = md.max(orb.density(x, y, z));
                    }
                    md
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap())
            .fold(1e-40_f64, f64::max)
    });

    // Rejection sampling, split across all cores
    let buf: Vec<(Vec3, f32, f32, f32, f64, f64, f64)> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|i| {
                let target = num_points / threads + usize::from(i < num_points % threads);
                s.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut part = Vec::with_capacity(target);
                    let mut attempts = 0usize;
                    while part.len() < target && attempts < target * 300 {
                        attempts += 1;
                        let x = rng.gen_range(-scale..scale);
                        let y = rng.gen_range(-scale..scale);
                        let z = rng.gen_range(-scale..scale);
                        let density = orb.density(x, y, z);
                        if rng.r#gen::<f64>() < density / max_density {
                            let (_, phase, j_mag) = orb.values(x, y, z);
                            let r = (x * x + y * y + z * z).sqrt().max(1e-10);
                            let rho_cyl = (x * x + y * y).sqrt() as f32;
                            let phi_raw = y.atan2(x);
                            let phi =
                                (if phi_raw < 0.0 { phi_raw + 2.0 * PI } else { phi_raw }) as f32;
                            let theta = (z / r).clamp(-1.0, 1.0).acos() as f32;
                            part.push((
                                Vec3::new(x as f32, y as f32, z as f32),
                                phi, rho_cyl, theta, density, phase, j_mag,
                            ));
                        }
                    }
                    part
                })
            })
            .collect();
        let mut all = Vec::with_capacity(num_points);
        for h in handles {
            all.extend(h.join().unwrap());
        }
        all
    });

    // Normalise density and j over the sampled set
    let max_d = buf.iter().map(|p| p.4).fold(0.0_f64, f64::max).max(1e-40);
    let max_j = buf.iter().map(|p| p.6).fold(0.0_f64, f64::max).max(1e-40);

    buf.into_iter()
        .map(|(pos, phi, rho_cyl, theta, d, phase, j)| RawPoint {
            pos,
            phi,
            rho_cyl,
            theta,
            density: (d / max_d).powf(0.35) as f32,
            phase: phase as f32,
            j_norm: (j / max_j).powf(0.45) as f32,
        })
        .collect()
}

// ── Superposition 1s + 2p_z ───────────────────────────────────────────────────
//
// ψ(t) = (ψ_1s·e^{-iE₁t} + ψ_2p·e^{-iE₂t})/√2 — les deux fonctions sont réelles:
//   |ψ(t)|² = ½(a² + b²) + a·b·cos(ω₁₂·t)   avec a = ψ_1s(x), b = ψ_2p(x)
// On échantillonne une seule fois l'enveloppe env = ½(a²+b²) + |ab| (le max
// temporel de la densité), puis chaque frame pondère les points existants par
// w = |ψ(t)|²/env ∈ [0,1] — aucun ré-échantillonnage pendant l'animation.

struct SupPoint {
    pos: Vec3,
    phi: f32,   // pour la coupe angulaire
    theta: f32, // pour la coupe angulaire
    a: f32,     // amplitude ψ_1s au point
    b: f32,     // amplitude ψ_2p_z au point
    env: f32,   // densité max sur une période
}

fn generate_sup_points(num_points: usize) -> Vec<SupPoint> {
    let o1 = Orbital::new(1, 0, 0, false);
    let o2 = Orbital::new(2, 1, 0, false);
    let scale = 14.0; // extension du 2p: 3.5·n² avec n=2
    let threads = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(4);

    let env_of = move |x: f64, y: f64, z: f64| -> (f64, f64, f64) {
        let a = o1.amplitude(x, y, z);
        let b = o2.amplitude(x, y, z);
        (a, b, 0.5 * (a * a + b * b) + (a * b).abs())
    };

    let max_env = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut md = 1e-40_f64;
                    for _ in 0..(100_000 / threads + 1) {
                        let x = rng.gen_range(-scale..scale);
                        let y = rng.gen_range(-scale..scale);
                        let z = rng.gen_range(-scale..scale);
                        md = md.max(env_of(x, y, z).2);
                    }
                    md
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().unwrap())
            .fold(1e-40_f64, f64::max)
    });

    let pts: Vec<SupPoint> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|i| {
                let target = num_points / threads + usize::from(i < num_points % threads);
                s.spawn(move || {
                    let mut rng = rand::thread_rng();
                    let mut part = Vec::with_capacity(target);
                    let mut attempts = 0usize;
                    while part.len() < target && attempts < target * 300 {
                        attempts += 1;
                        let x = rng.gen_range(-scale..scale);
                        let y = rng.gen_range(-scale..scale);
                        let z = rng.gen_range(-scale..scale);
                        let (a, b, env) = env_of(x, y, z);
                        if rng.r#gen::<f64>() < env / max_env {
                            let r = (x * x + y * y + z * z).sqrt().max(1e-10);
                            let phi_raw = y.atan2(x);
                            let phi =
                                (if phi_raw < 0.0 { phi_raw + 2.0 * PI } else { phi_raw }) as f32;
                            let theta = (z / r).clamp(-1.0, 1.0).acos() as f32;
                            part.push(SupPoint {
                                pos: Vec3::new(x as f32, y as f32, z as f32),
                                phi,
                                theta,
                                a: a as f32,
                                b: b as f32,
                                env: env as f32,
                            });
                        }
                    }
                    part
                })
            })
            .collect();
        let mut all = Vec::with_capacity(num_points);
        for h in handles {
            all.extend(h.join().unwrap());
        }
        all
    });
    pts
}

/// Une période du battement (fréquence de Bohr ω₁₂ = E₂−E₁ = 3/8 u.a.)
/// dure T_SUP secondes d'affichage à vitesse 1.
const T_SUP: f64 = 6.0;

fn build_sup_instances(
    pts: &[SupPoint],
    point_size: f32,
    clip_phi: f32,
    clip_theta: f32,
    anim_time: f64,
) -> Instances {
    let cosw = (2.0 * PI * anim_time / T_SUP).cos() as f32;
    let env_max = pts.iter().map(|p| p.env).fold(1e-30_f32, f32::max);

    let mut transforms = Vec::with_capacity(pts.len());
    let mut colors = Vec::with_capacity(pts.len());
    for p in pts {
        if clip_phi > 0.01 && p.phi <= clip_phi {
            continue;
        }
        if clip_theta > 0.01 && p.theta <= clip_theta {
            continue;
        }
        let dens = 0.5 * (p.a * p.a + p.b * p.b) + p.a * p.b * cosw;
        let w = (dens / p.env).clamp(0.0, 1.0);
        if w < 0.01 {
            continue;
        }
        // Rayon ∝ w^⅓: le volume affiché de chaque point suit |ψ(t)|²,
        // donc la « masse » visible du nuage oscille fidèlement.
        transforms.push(Mat4::from_translation(p.pos) * Mat4::from_scale(point_size * w.cbrt()));
        colors.push(viridis((dens / env_max).powf(0.35)));
    }
    Instances {
        transformations: transforms,
        colors: Some(colors),
        ..Default::default()
    }
}

// ── Instanced geometry ────────────────────────────────────────────────────────

fn build_instances(
    raw: &[RawPoint],
    point_size: f32,
    clip_phi: f32,   // radians
    clip_theta: f32, // radians
    mode: VisuMode,
    // Time evolution: ψ(t) = ψ₀ · e^{+it/(2n²)}  →  phase offset = ω·t, ω = 2π/(T·n²)
    anim_time: f64,  // display seconds elapsed
    n_val: u32,
    m_sign: f32, // ±1 = current direction, 0 = no current (m=0 or real orbital)
) -> Instances {
    let scale = Mat4::from_scale(point_size);
    let mut transforms = Vec::with_capacity(raw.len());
    let mut colors = Vec::with_capacity(raw.len());

    // T_BASE: display seconds for one phase revolution at n=1.
    // Physically ω_n = 1/(2n²) a.u.; here we rescale so n=1 completes in T_BASE seconds.
    // Also used for per-particle orbital speed in ProbCurrent mode.
    const T_BASE: f64 = 4.0;
    let omega: f32 = (2.0 * PI / (T_BASE * (n_val as f64).powi(2))) as f32;
    let phase_offset = (omega as f64 * anim_time) as f32; // total phase advance

    for p in raw {
        if clip_theta > 0.01 && p.theta <= clip_theta {
            continue;
        }

        let (t, c) = match mode {
            VisuMode::Density => {
                if clip_phi > 0.01 && p.phi <= clip_phi {
                    continue;
                }
                (Mat4::from_translation(p.pos) * scale, viridis(p.density))
            }
            VisuMode::WavePhase => {
                if clip_phi > 0.01 && p.phi <= clip_phi {
                    continue;
                }
                // Phase rotates: arg(ψ(t)) = arg(ψ₀) + ω·t
                let phase_t = (p.phase + phase_offset).rem_euclid(2.0 * PI as f32);
                (Mat4::from_translation(p.pos) * scale, phase_color(phase_t))
            }
            VisuMode::ProbCurrent => {
                // Particles physically orbit around Z.
                // Velocity field: v_φ = j_φ / |ψ|² = m / (r sinθ) = m / ρ_cyl  (a.u.)
                // Display: each particle revolves at  ω_i = m_sign · (2π/T_BASE) / ρ_cyl
                // → inner particles orbit faster (differential rotation, like a galaxy).
                // Only the SPEED divisor is clamped near the axis; the displayed
                // position keeps the true ρ_cyl, otherwise near-axis points get
                // pushed onto an artificial ρ=0.3 tube (ugly once the cloud is
                // sliced open).
                const DISPLAY_SPEED: f32 = 2.0 * std::f32::consts::PI / T_BASE as f32;
                let omega_i = m_sign * DISPLAY_SPEED / p.rho_cyl.max(0.3);
                let phi_t = (p.phi + omega_i * anim_time as f32).rem_euclid(2.0 * PI as f32);

                // La coupe φ s'applique à la position AFFICHÉE: la tranche
                // reste fixe dans l'espace et les particules la traversent.
                if clip_phi > 0.01 && phi_t <= clip_phi {
                    continue;
                }

                let rho = p.rho_cyl;
                let new_pos = Vec3::new(rho * phi_t.cos(), rho * phi_t.sin(), p.pos.z);

                // Ellipsoid stays aligned with the local azimuthal direction at phi_t
                let sp = phi_t.sin();
                let cp = phi_t.cos();
                let a = point_size * (1.0 + p.j_norm * 5.0);
                let b = point_size;
                let t = Mat4::new(
                    -sp * a,  cp * a, 0.0, 0.0,
                     cp * b,  sp * b, 0.0, 0.0,
                    0.0, 0.0, b, 0.0,
                    new_pos.x, new_pos.y, new_pos.z, 1.0,
                );
                (t, viridis(p.j_norm))
            }
        };
        transforms.push(t);
        colors.push(c);
    }

    Instances {
        transformations: transforms,
        colors: Some(colors),
        ..Default::default()
    }
}

/// XYZ axes.  Returns [positive_arrows, negative_stubs].
/// All primitives are along the +X axis in the CpuMesh convention
/// (cylinder and arrow both run from x=0 to x=1).
fn build_axes(context: &Context, length: f32) -> [Gm<InstancedMesh, ColorMaterial>; 2] {
    let r = (length / 55.0).max(0.04);

    // Rotations mapping +X to each positive axis direction
    let rx = Mat4::identity();
    let ry = Mat4::from_angle_z(Deg(90.0_f32));  // +X → +Y
    let rz = Mat4::from_angle_y(Deg(-90.0_f32)); // +X → +Z

    // Non-uniform scale: length along axis, r in cross-section
    let sc = |rot: Mat4| rot * Mat4::from_nonuniform_scale(length, r, r);

    let pos_inst = Instances {
        transformations: vec![sc(rx), sc(ry), sc(rz)],
        colors: Some(vec![
            Srgba::new(220, 60, 60, 255),  // X red
            Srgba::new(60, 200, 60, 255),  // Y green
            Srgba::new(60, 110, 230, 255), // Z blue
        ]),
        ..Default::default()
    };
    let pos_arrows = Gm::new(
        InstancedMesh::new(context, &pos_inst, &CpuMesh::arrow(0.85, 1.0, 12)),
        ColorMaterial { color: Srgba::WHITE, ..Default::default() },
    );

    // Rotations mapping +X to each negative axis direction
    let rnx = Mat4::from_angle_y(Deg(180.0_f32)); // +X → -X
    let rny = Mat4::from_angle_z(Deg(-90.0_f32)); // +X → -Y
    let rnz = Mat4::from_angle_y(Deg(90.0_f32));  // +X → -Z

    let rl = r * 0.4;
    let scn = |rot: Mat4| rot * Mat4::from_nonuniform_scale(length, rl, rl);

    let neg_inst = Instances {
        transformations: vec![scn(rnx), scn(rny), scn(rnz)],
        colors: Some(vec![
            Srgba::new(140, 40, 40, 180),
            Srgba::new(40, 130, 40, 180),
            Srgba::new(40, 60, 150, 180),
        ]),
        ..Default::default()
    };
    let neg_lines = Gm::new(
        InstancedMesh::new(context, &neg_inst, &CpuMesh::cylinder(8)),
        ColorMaterial { color: Srgba::WHITE, ..Default::default() },
    );

    [pos_arrows, neg_lines]
}

// ── Orbital naming ────────────────────────────────────────────────────────────

/// Lettre spectroscopique associée à l (s, p, d, f, g, h, i, ...).
fn subshell_letter(l: i32) -> char {
    match l {
        0 => 's',
        1 => 'p',
        2 => 'd',
        3 => 'f',
        _ => (b'g' + (l - 4) as u8) as char,
    }
}

/// Suffixe chimique des combinaisons réelles usuelles (p_x, d_xy, …).
fn real_orbital_suffix(l: i32, m: i32) -> Option<&'static str> {
    match (l, m) {
        (1, 0) => Some("z"),
        (1, 1) => Some("x"),
        (1, -1) => Some("y"),
        (2, 0) => Some("z²"),
        (2, 1) => Some("xz"),
        (2, -1) => Some("yz"),
        (2, 2) => Some("x²−y²"),
        (2, -2) => Some("xy"),
        (3, 0) => Some("z³"),
        (3, 1) => Some("xz²"),
        (3, -1) => Some("yz²"),
        (3, 2) => Some("z(x²−y²)"),
        (3, -2) => Some("xyz"),
        (3, 3) => Some("x(x²−3y²)"),
        (3, -3) => Some("y(3x²−y²)"),
        _ => None,
    }
}

/// Construit le nom correspondant à un triplet (n, l, m):
/// "2p (m=+1)" pour l'état propre complexe, "2p_x" pour la forme réelle.
fn orbital_label(n: u32, l: i32, m: i32, real: bool) -> String {
    if l == 0 {
        return format!("{n}s");
    }
    if real {
        match real_orbital_suffix(l, m) {
            Some(s) => format!("{n}{}_{s}", subshell_letter(l)),
            None => {
                let f = if m >= 0 { "cos" } else { "sin" };
                format!("{n}{} · {f}({}φ)", subshell_letter(l), m.abs())
            }
        }
    } else {
        format!("{n}{} (m={:+})", subshell_letter(l), m)
    }
}

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let window = Window::new(WindowSettings {
        title: "Orbitales de l'hydrogène".to_string(),
        max_size: Some((1400, 900)),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();

    // Vue « physicien »: l'axe de quantification Z est vertical, vers le haut.
    // Le up (0,0,1) est aussi respecté par OrbitControl (rotation à Z fixe,
    // façon tourne-disque), donc Z reste vertical quand on tourne la vue.
    let mut camera = Camera::new_perspective(
        window.viewport(),
        vec3(0.0, -47.0, 17.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 0.0, 1.0),
        degrees(45.0),
        0.1,
        2000.0,
    );
    let mut control = OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 500.0);
    let mut gui = GUI::new(&context);

    // ── State ─────────────────────────────────────────────────────────────────
    let num_points = 60_000usize;
    let mut n_val: u32 = 1;
    let mut l_val: i32 = 0;
    let mut m_val: i32 = 0;
    let mut real_mode = false; // orbitales réelles (p_x, d_xy, …)
    let mut superposition = false; // état (1s + 2p_z)/√2, dipôle oscillant
    let mut prev_sup = superposition;
    let mut prev_n = n_val;
    let mut prev_l = l_val;
    let mut prev_m = m_val;
    let mut prev_real = real_mode;
    let mut clip_phi_deg = 0.0_f32;
    let mut clip_theta_deg = 0.0_f32;
    let mut visu_mode = VisuMode::Density;

    // Previous-frame values to detect changes
    let mut prev_phi = -1.0_f32;
    let mut prev_theta = -1.0_f32;
    let mut prev_mode = VisuMode::Density;

    // Animation
    let mut anim_time: f64 = 0.0;
    let mut is_animating = false;
    let mut time_speed: f32 = 1.0;

    // Flags
    let mut regenerate = true; // true on first frame to initialise
    let mut rebuild = false;

    // Created once; instance buffers are updated in place afterwards
    // (recreating the InstancedMesh every frame was the main bottleneck).
    let empty_inst = Instances { transformations: vec![], ..Default::default() };
    let mut point_mesh = Gm::new(
        InstancedMesh::new(&context, &empty_inst, &CpuMesh::sphere(4)),
        ColorMaterial { color: Srgba::WHITE, ..Default::default() },
    );
    let mut axes = build_axes(&context, 10.0);
    let mut raw_points: Vec<RawPoint> = vec![];
    let mut sup_points: Vec<SupPoint> = vec![];

    window.render_loop(move |mut frame_input| {
        let mut redraw = frame_input.first_frame;
        redraw |= camera.set_viewport(frame_input.viewport);
        // ── GUI ───────────────────────────────────────────────────────────────
        // Traité AVANT le contrôle caméra: les panneaux consomment leurs
        // événements souris, sinon glisser un slider fait aussi tourner la vue.
        // Taille logique de la fenêtre, pour projeter les étiquettes X/Y/Z.
        let vp_w = frame_input.viewport.width as f32 / frame_input.device_pixel_ratio;
        let vp_h = frame_input.viewport.height as f32 / frame_input.device_pixel_ratio;
        redraw |= gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |ctx| {
                use egui::*;
                SidePanel::left("panel").min_width(180.0).show(ctx, |ui| {
                    ui.heading("Orbitale H");
                    ui.separator();

                    // Grisé quand la superposition est active (l'état n'est
                    // plus une orbitale propre unique)
                    ui.add_enabled_ui(!superposition, |ui| {
                    // Nombres quantiques n, l, m
                    ui.label("Nombres quantiques:");
                    ui.horizontal(|ui| {
                        ui.label("n");
                        ui.add(Slider::new(&mut n_val, 1..=6));
                    });
                    // l est contraint à [0, n-1]
                    if l_val > (n_val as i32 - 1) {
                        l_val = n_val as i32 - 1;
                    }
                    if l_val < 0 {
                        l_val = 0;
                    }
                    ui.horizontal(|ui| {
                        ui.label("l");
                        ui.add(Slider::new(&mut l_val, 0..=(n_val as i32 - 1)));
                    });
                    // m est contraint à [-l, l]
                    if m_val > l_val {
                        m_val = l_val;
                    }
                    if m_val < -l_val {
                        m_val = -l_val;
                    }
                    ui.horizontal(|ui| {
                        ui.label("m");
                        ui.add(Slider::new(&mut m_val, -l_val..=l_val));
                    });

                    ui.add_space(4.0);
                    ui.label(format!(
                        "Orbitale correspondante : {}",
                        orbital_label(n_val, l_val, m_val, real_mode)
                    ));
                    ui.small(format!(
                        "E = {:.2} eV · <r> = {:.1} a₀",
                        -13.6 / (n_val * n_val) as f64,
                        (3.0 * (n_val * n_val) as f64 - (l_val * (l_val + 1)) as f64) / 2.0
                    ));

                    ui.add_space(4.0);
                    ui.label("Forme de l'orbitale:");
                    ui.radio_value(&mut real_mode, false, "Complexe (état propre de L_z)");
                    ui.radio_value(&mut real_mode, true, "Réelle (chimie: p_x, d_xy, ...)");
                    ui.small(if real_mode {
                        "Superposition des états +m et -m:\nles lobes pointent selon x/y\n(les orbitales des livres de chimie)"
                    } else {
                        "Contient e^{imφ}: densité en anneau\nautour de Z, la phase tourne avec φ"
                    });
                    }); // fin zone grisée (superposition)

                    ui.separator();

                    // Superposition temporelle de deux états propres
                    ui.label("Superposition d'états:");
                    ui.checkbox(&mut superposition, "1s + 2p_z (dipôle oscillant)");
                    if superposition {
                        ui.small("ψ(t) = (ψ_1s·e^{-iE1·t} + ψ_2p·e^{-iE2·t})/√2");
                        ui.small(
                            "Deux énergies différentes: l'interférence\n\
                             bat à la fréquence de Bohr ν = (E2-E1)/h.\n\
                             Le nuage bascule entre +z et -z:\n\
                             un dipôle électrique oscillant.",
                        );
                        ui.small(
                            "C'est lui qui rayonne le photon Lyman-α\n\
                             (10,2 eV, λ = 121 nm) de la transition\n\
                             2p -> 1s.",
                        );
                    }

                    ui.separator();

                    // Visualisation mode
                    ui.add_enabled_ui(!superposition, |ui| {
                    ui.label("Mode:");
                    ui.radio_value(&mut visu_mode, VisuMode::Density,    "Densité |ψ|²");
                    ui.radio_value(&mut visu_mode, VisuMode::WavePhase,  "Fonction d'onde ψ");
                    ui.radio_value(&mut visu_mode, VisuMode::ProbCurrent,"Courant j");
                    }); // fin zone grisée (superposition)
                    if superposition {
                        ui.small("(superposition: densité |ψ(t)|², la taille\ndes points suit la densité instantanée)");
                    }

                    // Mode-dependent legend
                    if !superposition {
                    ui.add_space(4.0);
                    match visu_mode {
                        VisuMode::Density => {
                            ui.small("Couleur : densité de probabilité");
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(68,1,84),"■");
                                ui.label("faible");
                                ui.colored_label(Color32::from_rgb(253,231,37),"■");
                                ui.label("élevée");
                            });
                            ui.small("(densité complexe: axisymétrique pour m≠0)");
                        }
                        VisuMode::WavePhase => {
                            ui.small("Couleur : phase arg(ψ) = mφ");
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(255,0,0),"■"); ui.label("0");
                                ui.colored_label(Color32::from_rgb(0,200,0),"■"); ui.label("2π/3");
                                ui.colored_label(Color32::from_rgb(0,0,255),"■"); ui.label("4π/3");
                            });
                            if m_val == 0 {
                                ui.small("m=0: signe de ψ (0 ou π)");
                            } else if real_mode {
                                ui.small("Orbitale réelle: phase 0 ou π (signe des lobes)");
                            } else {
                                ui.small(format!("Phase tourne {} fois autour de Z", m_val.abs()));
                            }
                        }
                        VisuMode::ProbCurrent => {
                            ui.small("Couleur + élongation : |j|");
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(68,1,84),"■"); ui.label("faible");
                                ui.colored_label(Color32::from_rgb(253,231,37),"■"); ui.label("fort");
                            });
                            ui.small("Ellipsoïdes le long de la direction azimutale");
                            if m_val == 0 {
                                ui.small("m=0 : j = 0 partout");
                            } else if real_mode {
                                ui.small("Orbitale réelle : j = 0 partout");
                            } else if m_val > 0 {
                                ui.small("m>0: sens trigonométrique");
                            } else {
                                ui.small("m<0: sens horaire");
                            }
                        }
                    }
                    } // fin légende (masquée en superposition)

                    ui.separator();

                    // Animation
                    ui.label("Animation temporelle:");
                    ui.small("ψ(t) = ψ₀·e^{it/2n²}  (phase tourne à ω=1/2n²)");
                    ui.horizontal(|ui| {
                        let label = if is_animating { "⏸ Pause" } else { "▶ Play" };
                        if ui.button(label).clicked() {
                            is_animating = !is_animating;
                        }
                        if ui.button("↺ Reset").clicked() {
                            anim_time = 0.0;
                            rebuild = true;
                        }
                    });
                    ui.add(
                        Slider::new(&mut time_speed, 0.1_f32..=8.0)
                            .suffix("×").text("Vitesse"),
                    );
                    if superposition {
                        let period = T_SUP as f32 / time_speed;
                        ui.small(format!("Période du battement ≈ {period:.1}s"));
                    } else if visu_mode == VisuMode::Density {
                        ui.small("(densité: pas de dépendance temporelle)");
                    } else {
                        let period = 4.0 * (n_val as f32).powi(2) / time_speed;
                        ui.small(format!("Période ≈ {period:.1}s pour cette orbitale"));
                    }

                    ui.separator();

                    // Angular clipping
                    ui.label("Coupe angulaire:");
                    ui.add(Slider::new(&mut clip_phi_deg, 0.0_f32..=355.0)
                        .suffix("°").text("φ azimutal"));
                    ui.add(Slider::new(&mut clip_theta_deg, 0.0_f32..=175.0)
                        .suffix("°").text("θ polaire"));
                    ui.small("Retire les points dans [0°, angle]");

                    ui.separator();
                    ui.small("Axes: X rouge · Y vert · Z bleu");
                    ui.separator();
                    ui.small("Clic+glisser: rotation");
                    ui.small("Molette: zoom");
                });

                // ── Panneau d'analyse (graphiques) ────────────────────────
                SidePanel::right("plots").min_width(310.0).show(ctx, |ui| {
                    ScrollArea::vertical().show(ui, |ui| {
                        let orb = Orbital::new(n_val, l_val, m_val, real_mode);
                        let show_density = visu_mode == VisuMode::Density;
                        let mean_r =
                            (3.0 * (n_val * n_val) as f64 - (l_val * (l_val + 1)) as f64) / 2.0;

                        ui.heading("Analyse");
                        ui.separator();

                        if superposition {
                            // ── Superposition: les deux parties radiales ──
                            ui.label("Parties radiales R(r)");
                            let o1 = Orbital::new(1, 0, 0, false);
                            let o2 = Orbital::new(2, 1, 0, false);
                            let pts1: egui_plot::PlotPoints = (0..=256)
                                .map(|i| {
                                    let r = 16.0 * i as f64 / 256.0;
                                    [r, o1.radial(r)]
                                })
                                .collect();
                            let pts2: egui_plot::PlotPoints = (0..=256)
                                .map(|i| {
                                    let r = 16.0 * i as f64 / 256.0;
                                    [r, o2.radial(r)]
                                })
                                .collect();
                            egui_plot::Plot::new("sup_radial")
                                .height(170.0)
                                .allow_drag(false)
                                .allow_zoom(false)
                                .allow_scroll(false)
                                .allow_boxed_zoom(false)
                                .legend(egui_plot::Legend::default())
                                .show(ui, |pui| {
                                    pui.line(
                                        egui_plot::Line::new(pts1)
                                            .color(Color32::from_rgb(120, 190, 255))
                                            .width(1.8)
                                            .name("1s"),
                                    );
                                    pui.line(
                                        egui_plot::Line::new(pts2)
                                            .color(Color32::from_rgb(255, 160, 60))
                                            .width(1.8)
                                            .name("2p"),
                                    );
                                });
                            ui.small(
                                "Le battement vit là où les deux courbes\n\
                                 se recouvrent (r = 1 à 6 a₀ environ).",
                            );

                            ui.separator();
                            ui.label("Dipôle oscillant");
                            ui.small(
                                "|ψ(t)|² = ½(ψ_1s² + ψ_2p²)\n\
                                 + ψ_1s·ψ_2p·cos(ωt)",
                            );
                            ui.small(
                                "ψ_2p change de signe entre z>0 et z<0:\n\
                                 le terme croisé renforce la densité d'un\n\
                                 côté, la creuse de l'autre, puis s'inverse.\n\
                                 Le barycentre de la charge oscille le\n\
                                 long de z: une antenne dipolaire\n\
                                 microscopique.",
                            );
                        } else {

                        // ── Partie radiale ────────────────────────────────
                        if show_density {
                            ui.label("Densité radiale P(r) = r²R²(r)");
                        } else {
                            ui.label("Partie radiale R(r)");
                        }
                        let r_max = (n_val as f64).powi(2) * 3.5;
                        let pts: egui_plot::PlotPoints = (0..=256)
                            .map(|i| {
                                let r = r_max * i as f64 / 256.0;
                                let rr = orb.radial(r);
                                [r, if show_density { r * r * rr * rr } else { rr }]
                            })
                            .collect();
                        egui_plot::Plot::new("radial_plot")
                            .height(170.0)
                            .allow_drag(false)
                            .allow_zoom(false)
                            .allow_scroll(false)
                            .allow_boxed_zoom(false)
                            .show(ui, |pui| {
                                pui.line(
                                    egui_plot::Line::new(pts)
                                        .color(Color32::from_rgb(120, 190, 255))
                                        .width(1.8),
                                );
                                pui.vline(
                                    egui_plot::VLine::new(mean_r)
                                        .color(Color32::from_rgb(250, 200, 80))
                                        .style(egui_plot::LineStyle::dashed_loose()),
                                );
                            });
                        ui.small(format!(
                            "r en a₀ · <r> = {mean_r:.1} a₀ (pointillé)"
                        ));
                        ui.small(format!(
                            "Nœuds radiaux: {}",
                            n_val - l_val as u32 - 1
                        ));

                        ui.separator();

                        // ── Partie angulaire: diagramme polaire, plan x–z ──
                        if show_density {
                            ui.label("Partie angulaire |Y(θ)|² (coupe x-z)");
                        } else {
                            ui.label("Partie angulaire Y(θ) (coupe x-z)");
                        }
                        let size = ui.available_width().min(290.0);
                        let (prect, _) =
                            ui.allocate_exact_size(vec2(size, size), Sense::hover());
                        let painter = ui.painter_at(prect);
                        painter.rect_filled(prect, 4.0, Color32::from_rgb(13, 13, 24));
                        let c = prect.center();
                        let max_r = size * 0.44;
                        const NS: usize = 240;
                        let vals: Vec<f32> = (0..=NS)
                            .map(|i| {
                                let th = 2.0 * PI * i as f64 / NS as f64;
                                orb.angular(th.cos()) as f32
                            })
                            .collect();
                        let vmax =
                            vals.iter().fold(0.0_f32, |a, v| a.max(v.abs())).max(1e-12);
                        // Axes de la coupe: z vertical, x horizontal
                        painter.line_segment(
                            [pos2(c.x, prect.top() + 6.0), pos2(c.x, prect.bottom() - 6.0)],
                            Stroke::new(0.5, Color32::from_gray(70)),
                        );
                        painter.line_segment(
                            [pos2(prect.left() + 6.0, c.y), pos2(prect.right() - 6.0, c.y)],
                            Stroke::new(0.5, Color32::from_gray(70)),
                        );
                        painter.text(
                            pos2(c.x + 6.0, prect.top() + 10.0),
                            Align2::LEFT_CENTER,
                            "z",
                            FontId::proportional(11.0),
                            Color32::from_gray(160),
                        );
                        painter.text(
                            pos2(prect.right() - 10.0, c.y - 9.0),
                            Align2::CENTER_CENTER,
                            "x",
                            FontId::proportional(11.0),
                            Color32::from_gray(160),
                        );
                        // Rayon de démonstration: une direction θ mesurée depuis z
                        let th_demo = 0.7_f32;
                        let demo_dir = vec2(th_demo.sin(), -th_demo.cos());
                        painter.extend(Shape::dashed_line(
                            &[c, c + max_r * demo_dir],
                            Stroke::new(0.6, Color32::from_gray(110)),
                            4.0,
                            4.0,
                        ));
                        let arc: Vec<Pos2> = (0..=12)
                            .map(|i| {
                                let t = th_demo * i as f32 / 12.0;
                                c + 22.0 * vec2(t.sin(), -t.cos())
                            })
                            .collect();
                        painter.add(Shape::line(arc, Stroke::new(0.8, Color32::from_gray(130))));
                        painter.text(
                            c + 34.0 * vec2((th_demo * 0.55).sin(), -(th_demo * 0.55).cos()),
                            Align2::CENTER_CENTER,
                            "θ",
                            FontId::proportional(11.0),
                            Color32::from_gray(170),
                        );
                        // Lobes remplis: éventail de triangles depuis le centre
                        // (valide car r(θ) est une fonction radiale)
                        let radius = |v: f32| {
                            let t = (v / vmax).abs();
                            max_r * if show_density { t * t } else { t }
                        };
                        for i in 0..NS {
                            let th0 = 2.0 * PI * i as f64 / NS as f64;
                            let th1 = 2.0 * PI * (i + 1) as f64 / NS as f64;
                            let p0 = c + radius(vals[i])
                                * vec2(th0.sin() as f32, -(th0.cos() as f32));
                            let p1 = c + radius(vals[i + 1])
                                * vec2(th1.sin() as f32, -(th1.cos() as f32));
                            let (line_col, fill_col) = if show_density {
                                (
                                    Color32::from_rgb(180, 220, 100),
                                    Color32::from_rgba_unmultiplied(180, 220, 100, 45),
                                )
                            } else if vals[i] + vals[i + 1] >= 0.0 {
                                (
                                    Color32::from_rgb(255, 105, 90), // lobe Y > 0
                                    Color32::from_rgba_unmultiplied(255, 105, 90, 45),
                                )
                            } else {
                                (
                                    Color32::from_rgb(90, 140, 255), // lobe Y < 0
                                    Color32::from_rgba_unmultiplied(90, 140, 255, 45),
                                )
                            };
                            painter.add(Shape::convex_polygon(
                                vec![c, p0, p1],
                                fill_col,
                                Stroke::NONE,
                            ));
                            painter.line_segment([p0, p1], Stroke::new(1.6, line_col));
                        }
                        if !show_density {
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(255, 105, 90), "■");
                                ui.label("Y > 0");
                                ui.colored_label(Color32::from_rgb(90, 140, 255), "■");
                                ui.label("Y < 0");
                            });
                        }
                        ui.small(
                            "Lecture: dans chaque direction θ (angle depuis z),\n\
                             la distance courbe-centre donne le poids de\n\
                             l'orbitale dans cette direction.",
                        );
                        ui.small(if real_mode && m_val != 0 {
                            "Nuage 3D = cette coupe × R(r), modulée\n\
                             en azimut par cos/sin(mφ) (forme réelle)."
                        } else {
                            "Nuage 3D = cette coupe × R(r), identique\n\
                             dans tous les plans contenant z."
                        });
                        let mut nodes = format!(
                            "Nœuds angulaires: {} cône(s)",
                            l_val - m_val.abs()
                        );
                        if real_mode && m_val != 0 {
                            nodes += &format!(" + {} plan(s)", m_val.abs());
                        }
                        ui.small(nodes);

                        } // fin état propre (le diagramme d'énergie est commun)

                        ui.separator();

                        // ── Niveaux d'énergie ─────────────────────────────
                        ui.label("Niveaux E_n = -13,6/n² eV (cliquer)");
                        let w = ui.available_width().min(290.0);
                        let (erect, eresp) =
                            ui.allocate_exact_size(vec2(w, 190.0), Sense::click());
                        let painter = ui.painter_at(erect);
                        painter.rect_filled(erect, 4.0, Color32::from_rgb(13, 13, 24));
                        let y_of =
                            |e: f32| erect.top() + 8.0 + (erect.height() - 16.0) * (e / -14.0);
                        let y0 = y_of(0.0);
                        painter.line_segment(
                            [pos2(erect.left() + 42.0, y0), pos2(erect.right() - 8.0, y0)],
                            Stroke::new(0.5, Color32::from_gray(90)),
                        );
                        painter.text(
                            pos2(erect.left() + 38.0, y0),
                            Align2::RIGHT_CENTER,
                            "E=0",
                            FontId::proportional(9.0),
                            Color32::from_gray(120),
                        );
                        for nn in 1..=6u32 {
                            let y = y_of(-13.6 / (nn * nn) as f32);
                            let sel = if superposition { nn <= 2 } else { nn == n_val };
                            let col = if sel {
                                Color32::from_rgb(255, 210, 80)
                            } else {
                                Color32::from_gray(150)
                            };
                            painter.line_segment(
                                [pos2(erect.left() + 42.0, y), pos2(erect.right() - 8.0, y)],
                                Stroke::new(if sel { 2.5 } else { 1.0 }, col),
                            );
                            painter.text(
                                pos2(erect.left() + 38.0, y),
                                Align2::RIGHT_CENTER,
                                format!("n={nn}"),
                                FontId::proportional(9.0),
                                col,
                            );
                        }
                        if superposition {
                            // Flèche de transition entre les deux niveaux
                            let x = erect.left() + 90.0;
                            let (y1, y2) = (y_of(-13.6), y_of(-3.4));
                            painter.line_segment(
                                [pos2(x, y2), pos2(x, y1)],
                                Stroke::new(1.2, Color32::from_rgb(255, 160, 60)),
                            );
                            painter.text(
                                pos2(x + 5.0, (y1 + y2) * 0.5),
                                Align2::LEFT_CENTER,
                                "10,2 eV",
                                FontId::proportional(10.0),
                                Color32::from_rgb(255, 160, 60),
                            );
                        }
                        if !superposition && eresp.clicked() {
                            if let Some(p) = eresp.interact_pointer_pos() {
                                let mut best: Option<(f32, u32)> = None;
                                for nn in 1..=6u32 {
                                    let d = (y_of(-13.6 / (nn * nn) as f32) - p.y).abs();
                                    if d < 10.0 && best.is_none_or(|(bd, _)| d < bd) {
                                        best = Some((d, nn));
                                    }
                                }
                                if let Some((_, nn)) = best {
                                    n_val = nn;
                                }
                            }
                        }
                        if superposition {
                            ui.small(
                                "Transition 2p -> 1s: ΔE = 10,2 eV,\n\
                                 photon Lyman-α (λ = 121 nm)",
                            );
                        } else {
                            ui.small(format!(
                                "E{} = {:.2} eV · dégénérescence: {} orbitales",
                                n_val,
                                -13.6 / (n_val * n_val) as f64,
                                n_val * n_val
                            ));
                        }
                    });
                });

                // ── Noms des axes 3D, projetés aux extrémités des flèches ──
                {
                    use egui::*;
                    let painter = ctx.layer_painter(LayerId::new(
                        Order::Background,
                        Id::new("axis_labels"),
                    ));
                    let tip = (n_val as f32).powi(2) * 4.2 * 1.07;
                    let mvp = camera.projection() * camera.view();
                    for (dir, name, col) in [
                        (three_d::vec3(1.0_f32, 0.0, 0.0), "X", Color32::from_rgb(235, 95, 95)),
                        (three_d::vec3(0.0, 1.0, 0.0), "Y", Color32::from_rgb(95, 220, 95)),
                        (three_d::vec3(0.0, 0.0, 1.0), "Z", Color32::from_rgb(110, 155, 245)),
                    ] {
                        let p = dir * tip;
                        let clip = mvp * three_d::vec4(p.x, p.y, p.z, 1.0);
                        if clip.w > 0.1 {
                            let sx = (clip.x / clip.w * 0.5 + 0.5) * vp_w;
                            let sy = (1.0 - (clip.y / clip.w * 0.5 + 0.5)) * vp_h;
                            painter.text(
                                pos2(sx, sy),
                                Align2::CENTER_CENTER,
                                name,
                                FontId::proportional(14.0),
                                col,
                            );
                        }
                    }
                }
            },
        );

        // La caméra ne reçoit que les événements restants
        redraw |= control.handle_events(&mut camera, &mut frame_input.events);

        // ── Detect GUI-driven changes ─────────────────────────────────────────
        // n peut aussi changer via un clic sur le diagramme d'énergie:
        // on re-contraint l et m ici avant de régénérer.
        l_val = l_val.clamp(0, n_val as i32 - 1);
        m_val = m_val.clamp(-l_val, l_val);
        if n_val != prev_n || l_val != prev_l || m_val != prev_m || real_mode != prev_real {
            prev_n = n_val;
            prev_l = l_val;
            prev_m = m_val;
            prev_real = real_mode;
            regenerate = true;
        }
        if superposition != prev_sup {
            prev_sup = superposition;
            regenerate = true;
            if superposition {
                // Démarre à cos(ωt)=1: densité déportée vers +z
                anim_time = 0.0;
                is_animating = true;
            }
        }
        if visu_mode != prev_mode {
            // Le mode Courant se met en mouvement sans passer par Play
            if visu_mode == VisuMode::ProbCurrent {
                is_animating = true;
            }
            prev_mode = visu_mode;
            rebuild = true;
        }
        if (clip_phi_deg - prev_phi).abs() > 0.4 || (clip_theta_deg - prev_theta).abs() > 0.4 {
            prev_phi = clip_phi_deg;
            prev_theta = clip_theta_deg;
            rebuild = true;
        }

        // ── Advance animation clock ───────────────────────────────────────────
        // La superposition anime la densité elle-même, quel que soit le mode.
        let animates_visually =
            is_animating && (superposition || visu_mode != VisuMode::Density);
        if animates_visually {
            anim_time += frame_input.elapsed_time / 1000.0 * time_speed as f64;
            rebuild = true; // rebuild colours every frame
            redraw = true;
        }

        // ── Full regeneration (new orbital) ───────────────────────────────────
        if regenerate {
            if superposition {
                sup_points = generate_sup_points(num_points);
                let inst = build_sup_instances(
                    &sup_points, 0.30,
                    clip_phi_deg.to_radians(), clip_theta_deg.to_radians(),
                    anim_time,
                );
                point_mesh.geometry.set_instances(&inst);
                axes = build_axes(&context, 16.8); // extension du 2p (n=2)
            } else {
                let orb = Orbital::new(n_val, l_val, m_val, real_mode);
                let ps = 0.15 * n_val as f32;
                // Les axes dépassent du nuage (échantillonné jusqu'à 3.5·n²)
                let al = (n_val as f32).powi(2) * 4.2;
                raw_points = generate_raw_points(orb, num_points);
                let m_sign = if real_mode { 0.0 } else { m_val.signum() as f32 };
                let inst = build_instances(
                    &raw_points, ps,
                    clip_phi_deg.to_radians(), clip_theta_deg.to_radians(),
                    visu_mode, anim_time, n_val, m_sign,
                );
                point_mesh.geometry.set_instances(&inst);
                axes = build_axes(&context, al);
            }
            regenerate = false;
            rebuild = false;
            redraw = true;
        }

        // ── Fast rebuild (clip / mode / animation change, no resampling) ──────
        if rebuild {
            let inst = if superposition {
                build_sup_instances(
                    &sup_points, 0.30,
                    clip_phi_deg.to_radians(), clip_theta_deg.to_radians(),
                    anim_time,
                )
            } else {
                let ps = 0.15 * n_val as f32;
                let m_sign = if real_mode { 0.0 } else { m_val.signum() as f32 };
                build_instances(
                    &raw_points, ps,
                    clip_phi_deg.to_radians(), clip_theta_deg.to_radians(),
                    visu_mode, anim_time, n_val, m_sign,
                )
            };
            point_mesh.geometry.set_instances(&inst);
            rebuild = false;
            redraw = true;
        }

        // ── Render ────────────────────────────────────────────────────────────
        if redraw {
            let screen = frame_input.screen();
            screen.clear(ClearState::color_and_depth(0.02, 0.02, 0.05, 1.0, 1.0));
            screen.render(&camera, [&point_mesh, &axes[0], &axes[1]], &[]);
            screen.write(|| gui.render()).unwrap();
        }

        FrameOutput { swap_buffers: redraw, ..Default::default() }
    });
}
