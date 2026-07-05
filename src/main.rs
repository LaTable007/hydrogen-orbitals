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

/// Returns (density, phase ∈ [0,2π], |j_φ|) for the *complex* eigenstate ψ_{nlm}.
///
/// |ψ|² = R² · N² · P²(cosθ)   — no φ dependence (axisymmetric for m≠0)
/// arg(ψ)  = m·φ  +  π·H(−R·N·P)   where H is the Heaviside function
/// j_φ  = m · |ψ|² / (r sinθ)   →  |j| = |m| · |ψ|² / ρ_cyl
fn hydrogen_values(n: u32, l: i32, m: i32, x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let r = (x * x + y * y + z * z).sqrt();
    if r < 1e-10 {
        return (0.0, 0.0, 0.0);
    }
    let theta = (z / r).clamp(-1.0, 1.0).acos();
    let phi = y.atan2(x); // [-π, π]
    let rho = 2.0 * r / n as f64;
    let l_u = l as u32;
    let m_abs = m.unsigned_abs();

    // Radial
    let norm_r_sq = (2.0 / n as f64).powi(3) * factorial(n - l_u - 1)
        / (2.0 * n as f64 * factorial(n + l_u).powi(3));
    let radial = norm_r_sq.sqrt()
        * (-rho / 2.0).exp()
        * rho.powi(l)
        * laguerre_associated(n, l_u, rho);

    // Angular amplitude (real, without e^{imφ})
    let norm_a = (((2 * l + 1) as f64 / (4.0 * PI))
        * (factorial(l as u32 - m_abs) / factorial(l as u32 + m_abs)))
        .sqrt();
    let plm = associated_legendre(l, m_abs, theta.cos());
    let amp = radial * norm_a * plm; // signed real amplitude

    let density = amp * amp;

    // Phase:  arg(ψ) = m·φ + (π if amp < 0)
    let phase_base = if amp >= 0.0 { 0.0 } else { PI };
    let phase = ((phase_base + m as f64 * phi) % (2.0 * PI) + 2.0 * PI) % (2.0 * PI);

    // Probability current: j_φ = m·|ψ|²/(r·sinθ) = m·density/ρ_cyl
    let rho_cyl = (x * x + y * y).sqrt().max(1e-10);
    let j_mag = m_abs as f64 * density / rho_cyl;

    (density, phase, j_mag)
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
    phi: f32,     // [0, 2π] — used for azimuthal clipping
    theta: f32,   // [0, π]  — used for polar clipping
    density: f32, // normalised to [0,1] with γ-correction
    phase: f32,   // [0, 2π]
    j_norm: f32,  // |j| normalised [0,1]
}

fn generate_raw_points(n: u32, l: i32, m: i32, num_points: usize) -> Vec<RawPoint> {
    let mut rng = rand::thread_rng();
    let scale = (n as f64).powi(2) * 3.5;

    // Estimate max density for rejection sampling
    let mut max_density = 1e-40_f64;
    for _ in 0..100_000 {
        let x = rng.gen_range(-scale..scale);
        let y = rng.gen_range(-scale..scale);
        let z = rng.gen_range(-scale..scale);
        let (d, _, _) = hydrogen_values(n, l, m, x, y, z);
        if d > max_density {
            max_density = d;
        }
    }

    // Rejection sampling
    let mut buf: Vec<(Vec3, f32, f32, f64, f64, f64)> = Vec::with_capacity(num_points);
    let mut attempts = 0usize;
    while buf.len() < num_points && attempts < num_points * 300 {
        attempts += 1;
        let x = rng.gen_range(-scale..scale);
        let y = rng.gen_range(-scale..scale);
        let z = rng.gen_range(-scale..scale);
        let (density, phase, j_mag) = hydrogen_values(n, l, m, x, y, z);
        if rng.r#gen::<f64>() < density / max_density {
            let r = (x * x + y * y + z * z).sqrt().max(1e-10);
            let phi_raw = y.atan2(x);
            let phi = (if phi_raw < 0.0 { phi_raw + 2.0 * PI } else { phi_raw }) as f32;
            let theta = (z / r).clamp(-1.0, 1.0).acos() as f32;
            buf.push((Vec3::new(x as f32, y as f32, z as f32), phi, theta, density, phase, j_mag));
        }
    }

    // Normalise density and j over the sampled set
    let max_d = buf.iter().map(|p| p.3).fold(0.0_f64, f64::max).max(1e-40);
    let max_j = buf.iter().map(|p| p.5).fold(0.0_f64, f64::max).max(1e-40);

    buf.into_iter()
        .map(|(pos, phi, theta, d, phase, j)| RawPoint {
            pos,
            phi,
            theta,
            density: (d / max_d).powf(0.35) as f32,
            phase: phase as f32,
            j_norm: (j / max_j).powf(0.45) as f32,
        })
        .collect()
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
    m_val: i32,
) -> Instances {
    let scale = Mat4::from_scale(point_size);
    let mut transforms = Vec::new();
    let mut colors = Vec::new();

    // Display angular frequency: one full phase revolution every T_base·n² seconds.
    // Physically: ω_n = 1/(2n²) a.u.; we rescale to human timescale.
    const T_BASE: f64 = 4.0; // seconds per revolution for n=1
    let omega: f32 = (2.0 * PI / (T_BASE * (n_val as f64).powi(2))) as f32;
    let phase_offset = (omega as f64 * anim_time) as f32; // total phase advance

    let m_sign = m_val.signum() as f32; // +1, 0, or -1

    for p in raw {
        if clip_phi > 0.01 && p.phi <= clip_phi {
            continue;
        }
        if clip_theta > 0.01 && p.theta <= clip_theta {
            continue;
        }

        let (t, c) = match mode {
            VisuMode::Density => (
                Mat4::from_translation(p.pos) * scale,
                viridis(p.density),
            ),
            VisuMode::WavePhase => {
                // Phase rotates: arg(ψ(t)) = arg(ψ₀) + ω·t
                let phase_t = (p.phase + phase_offset).rem_euclid(2.0 * PI as f32);
                (Mat4::from_translation(p.pos) * scale, phase_color(phase_t))
            }
            VisuMode::ProbCurrent => {
                // |j| is time-independent (stationary state), but we animate a
                // "tracer" spotlight that sweeps in the direction of actual current flow:
                //   m > 0 → CCW,  m < 0 → CW.
                // Bright spot at φ_spot = phase_offset / |m_sign|, moving in +m direction.
                let flow_mod = if m_val == 0 {
                    1.0_f32
                } else {
                    // cos argument = 0 at φ_spot; spot moves CCW for m>0, CW for m<0
                    0.35 + 0.65 * (m_sign * p.phi - phase_offset).cos().abs()
                };

                let sp = p.phi.sin();
                let cp = p.phi.cos();
                let a = point_size * (1.0 + p.j_norm * 5.0);
                let b = point_size;
                let t = Mat4::new(
                    -sp * a,  cp * a, 0.0, 0.0,
                     cp * b,  sp * b, 0.0, 0.0,
                    0.0, 0.0, b, 0.0,
                    p.pos.x, p.pos.y, p.pos.z, 1.0,
                );
                (t, viridis(p.j_norm * flow_mod))
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

fn build_point_mesh(context: &Context, inst: &Instances) -> Gm<InstancedMesh, ColorMaterial> {
    Gm::new(
        InstancedMesh::new(context, inst, &CpuMesh::sphere(4)),
        ColorMaterial { color: Srgba::WHITE, ..Default::default() },
    )
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

// ── Orbital table ─────────────────────────────────────────────────────────────

const ORBITALS: &[(u32, i32, i32, &str)] = &[
    (1, 0, 0,  "1s"),
    (2, 0, 0,  "2s"),
    (2, 1, -1, "2p (m=-1)"),
    (2, 1, 0,  "2p (m=0)"),
    (2, 1, 1,  "2p (m=+1)"),
    (3, 0, 0,  "3s"),
    (3, 1, -1, "3p (m=-1)"),
    (3, 1, 0,  "3p (m=0)"),
    (3, 1, 1,  "3p (m=+1)"),
    (3, 2, -2, "3d (m=-2)"),
    (3, 2, -1, "3d (m=-1)"),
    (3, 2, 0,  "3d (m=0)"),
    (3, 2, 1,  "3d (m=+1)"),
    (3, 2, 2,  "3d (m=+2)"),
    (4, 0, 0,  "4s"),
    (4, 1, 0,  "4p (m=0)"),
    (4, 2, 0,  "4d (m=0)"),
    (4, 3, 0,  "4f (m=0)"),
];

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() {
    let window = Window::new(WindowSettings {
        title: "Orbitales de l'hydrogène".to_string(),
        max_size: Some((1400, 900)),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();

    let mut camera = Camera::new_perspective(
        window.viewport(),
        vec3(0.0, 0.0, 50.0),
        vec3(0.0, 0.0, 0.0),
        vec3(0.0, 1.0, 0.0),
        degrees(45.0),
        0.1,
        2000.0,
    );
    let mut control = OrbitControl::new(vec3(0.0, 0.0, 0.0), 1.0, 500.0);
    let mut gui = GUI::new(&context);

    // ── State ─────────────────────────────────────────────────────────────────
    let num_points = 60_000usize;
    let mut orbital_idx = 0usize;
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

    // Placeholders (replaced immediately by regenerate=true)
    let empty_inst = Instances { transformations: vec![], ..Default::default() };
    let empty_cpu = CpuMesh::sphere(1);
    let mut point_mesh = Gm::new(
        InstancedMesh::new(&context, &empty_inst, &empty_cpu),
        ColorMaterial { color: Srgba::WHITE, ..Default::default() },
    );
    let mut axes = build_axes(&context, 10.0);
    let mut raw_points: Vec<RawPoint> = vec![];

    window.render_loop(move |mut frame_input| {
        let mut redraw = frame_input.first_frame;
        redraw |= camera.set_viewport(frame_input.viewport);
        redraw |= control.handle_events(&mut camera, &mut frame_input.events);

        // ── Advance animation clock ───────────────────────────────────────────
        let animates_visually = is_animating && visu_mode != VisuMode::Density;
        if animates_visually {
            anim_time += frame_input.elapsed_time / 1000.0 * time_speed as f64;
            rebuild = true; // rebuild colours every frame
            redraw = true;
        }

        // ── Full regeneration (new orbital) ───────────────────────────────────
        if regenerate {
            let (n, l, m, _) = ORBITALS[orbital_idx];
            let ps = 0.15 * n as f32;
            let al = (n as f32).powi(2) * 2.8;
            raw_points = generate_raw_points(n, l, m, num_points);
            let inst = build_instances(
                &raw_points, ps,
                clip_phi_deg.to_radians(), clip_theta_deg.to_radians(),
                visu_mode, anim_time, n, m,
            );
            point_mesh = build_point_mesh(&context, &inst);
            axes = build_axes(&context, al);
            regenerate = false;
            rebuild = false;
            redraw = true;
        }

        // ── Fast rebuild (clip / mode / animation change, no resampling) ──────
        if rebuild {
            let (n, _, m, _) = ORBITALS[orbital_idx];
            let ps = 0.15 * n as f32;
            let inst = build_instances(
                &raw_points, ps,
                clip_phi_deg.to_radians(), clip_theta_deg.to_radians(),
                visu_mode, anim_time, n, m,
            );
            point_mesh = build_point_mesh(&context, &inst);
            rebuild = false;
            redraw = true;
        }

        // ── GUI ───────────────────────────────────────────────────────────────
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |ctx| {
                use egui::*;
                SidePanel::left("panel").min_width(180.0).show(ctx, |ui| {
                    ui.heading("Orbitale H");
                    ui.separator();

                    // Orbital selector
                    ui.label("Orbitale:");
                    ScrollArea::vertical().max_height(220.0).show(ui, |ui| {
                        for (i, &(_, _, _, label)) in ORBITALS.iter().enumerate() {
                            if ui.selectable_label(i == orbital_idx, label).clicked()
                                && i != orbital_idx
                            {
                                orbital_idx = i;
                                regenerate = true;
                            }
                        }
                    });

                    ui.separator();

                    // Visualisation mode
                    ui.label("Mode:");
                    ui.radio_value(&mut visu_mode, VisuMode::Density,    "Densité |ψ|²");
                    ui.radio_value(&mut visu_mode, VisuMode::WavePhase,  "Fonction d'onde ψ");
                    ui.radio_value(&mut visu_mode, VisuMode::ProbCurrent,"Courant j");

                    // Mode-dependent legend
                    ui.add_space(4.0);
                    match visu_mode {
                        VisuMode::Density => {
                            ui.small("Couleur → densité de probabilité");
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(68,1,84),"■");
                                ui.label("faible");
                                ui.colored_label(Color32::from_rgb(253,231,37),"■");
                                ui.label("élevée");
                            });
                            ui.small("(densité complexe: axisymétrique pour m≠0)");
                        }
                        VisuMode::WavePhase => {
                            ui.small("Couleur → phase arg(ψ) = mφ");
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(255,0,0),"■"); ui.label("0");
                                ui.colored_label(Color32::from_rgb(0,200,0),"■"); ui.label("2π/3");
                                ui.colored_label(Color32::from_rgb(0,0,255),"■"); ui.label("4π/3");
                            });
                            let (_, _, m, _) = ORBITALS[orbital_idx];
                            if m == 0 {
                                ui.small("m=0: signe de ψ (0 ou π)");
                            } else {
                                ui.small(format!("Phase tourne {} fois autour de Z", m.abs()));
                            }
                        }
                        VisuMode::ProbCurrent => {
                            let (_, _, m, _) = ORBITALS[orbital_idx];
                            ui.small("Couleur + élongation → |j|");
                            ui.horizontal(|ui| {
                                ui.colored_label(Color32::from_rgb(68,1,84),"■"); ui.label("faible");
                                ui.colored_label(Color32::from_rgb(253,231,37),"■"); ui.label("fort");
                            });
                            ui.small("Ellipsoïde orienté selon φ̂ (azimutal)");
                            if m == 0 {
                                ui.small("m=0 → j = 0 partout");
                            } else if m > 0 {
                                ui.small("m>0: sens trigonométrique");
                            } else {
                                ui.small("m<0: sens horaire");
                            }
                        }
                    }

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
                    if visu_mode == VisuMode::Density {
                        ui.small("(densité: pas de dépendance temporelle)");
                    } else {
                        let (n, _, _, _) = ORBITALS[orbital_idx];
                        let period = 4.0 * (n as f32).powi(2) / time_speed;
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
            },
        );

        // ── Detect GUI-driven changes ─────────────────────────────────────────
        if !regenerate {
            let changed = (clip_phi_deg - prev_phi).abs() > 0.4
                || (clip_theta_deg - prev_theta).abs() > 0.4
                || visu_mode != prev_mode;
            if changed {
                rebuild = true;
                prev_phi = clip_phi_deg;
                prev_theta = clip_theta_deg;
                prev_mode = visu_mode;
            }
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
