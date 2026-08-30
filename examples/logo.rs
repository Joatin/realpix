//! Draws the project logo — and does it with the library, so the tessellation in the mark
//! is a real HEALPix grid rather than a picture of one.
//!
//! Run with `cargo run --example logo`; it writes `assets/logo.svg`.
//!
//! The sphere is drawn in orthographic projection: a cell's corners come from
//! [`vertices`](realpix::nested::Layer::vertices), each edge is subdivided along its great
//! circle so the curvature survives the projection, and every edge shared with a different
//! base cell is stroked more brightly to bring out the twelve HEALPix faces. The amber
//! patch is an actual `cone_coverage` result.

use realpix::nested;
use std::fmt::Write as _;

/// Depth of the drawn grid. Twelve base cells of 16 each: bold enough to read at favicon
/// size, fine enough to show the tessellation.
const DEPTH: u8 = 2;
/// Depth the highlighted cone is resolved at. Drawing the search four levels below the
/// grid keeps its edge round, and is a fair picture of how the crate is actually used: the
/// coverage does not have to share a depth with anything else.
const CONE_DEPTH: u8 = 4;
/// Samples per cell edge. Cell edges are not great circles, but over a depth-2 edge the
/// difference is far below a pixel and the curvature is what matters.
const STEPS: usize = 6;
/// Samples for the twelve base-cell boundaries, which are long and the most prominent
/// lines in the mark, so they get a finer curve.
const BOLD_STEPS: usize = 12;

const SIZE: f64 = 512.0;
const R: f64 = 232.0;

/// Direction the sphere is viewed from, as (lon, lat) in radians. Tilted so that the
/// north polar cap, the equatorial belt and a hint of the south are all in view.
const VIEW: (f64, f64) = (0.85, 0.42);

/// Centre and angular radius of the highlighted cone search.
const CONE: (f64, f64, f64) = (1.40, 0.58, 0.25);

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = dot(v, v).sqrt();
    [v[0] / n, v[1] / n, v[2] / n]
}

/// Great-circle interpolation, so a subdivided edge stays on the sphere.
fn slerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    let c = dot(a, b).clamp(-1.0, 1.0);
    let omega = c.acos();
    if omega < 1e-9 {
        return a;
    }
    let (s, sa, sb) = (omega.sin(), ((1.0 - t) * omega).sin(), (t * omega).sin());
    unit([
        (a[0] * sa + b[0] * sb) / s,
        (a[1] * sa + b[1] * sb) / s,
        (a[2] * sa + b[2] * sb) / s,
    ])
}

struct View {
    axis: [f64; 3],
    right: [f64; 3],
    up: [f64; 3],
}

impl View {
    fn new(lon: f64, lat: f64) -> View {
        let axis = realpix::lonlat_to_vec(lon, lat);
        let right = unit(cross([0.0, 0.0, 1.0], axis));
        View {
            axis,
            right,
            up: cross(axis, right),
        }
    }

    /// Projects a direction to canvas coordinates. Points on the far side are pushed out
    /// to the limb rather than folded across it, so a cell straddling the edge still draws
    /// as a closed shape; the surrounding clip path trims it.
    fn project(&self, p: [f64; 3]) -> (f64, f64) {
        let (mut x, mut y) = (dot(p, self.right), dot(p, self.up));
        if dot(p, self.axis) < 0.0 {
            let n = (x * x + y * y).sqrt().max(1e-12);
            x /= n;
            y /= n;
        }
        (SIZE / 2.0 + R * x, SIZE / 2.0 - R * y)
    }
}

/// The four edges of a cell, each as the neighbour it is shared with. `vertices` returns
/// the corners as [N, W, S, E], so the edge N→W is the one the north-west neighbour is on,
/// and so round.
const EDGE_NEIGHBOUR: [realpix::Direction; 4] = [
    realpix::Direction::NW,
    realpix::Direction::SW,
    realpix::Direction::SE,
    realpix::Direction::NE,
];

fn main() -> std::io::Result<()> {
    let layer = nested::get(DEPTH);
    let view = View::new(VIEW.0, VIEW.1);

    let cone_center = realpix::lonlat_to_vec(CONE.0, CONE.1);
    let cone_layer = nested::get(CONE_DEPTH);
    let covered = cone_layer.cone_coverage_cells(cone_center, CONE.2);

    let mut cells = String::new();
    let mut faint = String::new();
    let mut bold = String::new();
    let mut cone_fill = String::new();
    let mut cone_edge = String::new();

    for cell in layer.iter() {
        let corners = layer.vertices(cell);
        // Cull the far hemisphere, keeping cells that merely straddle the limb.
        if dot(layer.center_vec(cell), view.axis) < -0.15 {
            continue;
        }

        // The filled cell, with its edges subdivided along their great circles.
        let mut path = String::new();
        for edge in 0..4 {
            let (a, b) = (corners[edge], corners[(edge + 1) % 4]);
            for step in 0..STEPS {
                let (x, y) = view.project(slerp(a, b, step as f64 / STEPS as f64));
                let _ = write!(
                    path,
                    "{}{x:.1} {y:.1} ",
                    if edge == 0 && step == 0 { "M" } else { "L" }
                );
            }
        }
        path.push('Z');

        // Shade by latitude so the sphere reads as curved, and light the side facing the
        // viewer. Cells inside the cone take the accent instead.
        let lit = dot(layer.center_vec(cell), view.axis).max(0.0);
        let z = layer.center_vec(cell)[2];
        let t = (0.22 + 0.78 * lit) * (0.70 + 0.30 * (1.0 - z.abs()));
        let fill = format!(
            "rgb({},{},{})",
            (16.0 + 54.0 * t) as u32,
            (30.0 + 100.0 * t) as u32,
            (70.0 + 152.0 * t) as u32
        );
        let _ = writeln!(cells, r#"    <path d="{path}" fill="{fill}"/>"#);

        // Stroke each edge once, brightly where it borders another base cell.
        for (edge, direction) in EDGE_NEIGHBOUR.iter().enumerate() {
            let neighbour = layer.neighbour(cell, *direction);
            let base_edge = match neighbour {
                None => true,
                Some(n) => layer.parent(n, 0) != layer.parent(cell, 0),
            };
            // Interior edges are shared by two cells; draw them from the lower index only.
            if !base_edge && neighbour.is_some_and(|n| n < cell) {
                continue;
            }
            let (a, b) = (corners[edge], corners[(edge + 1) % 4]);
            let n = if base_edge { BOLD_STEPS } else { STEPS };
            let mut seg = String::new();
            for step in 0..=n {
                let (x, y) = view.project(slerp(a, b, step as f64 / n as f64));
                let _ = write!(seg, "{}{x:.1} {y:.1} ", if step == 0 { "M" } else { "L" });
            }
            let _ = writeln!(
                if base_edge { &mut bold } else { &mut faint },
                r#"    <path d="{}"/>"#,
                seg.trim_end()
            );
        }
    }

    for &cell in &covered {
        if dot(cone_layer.center_vec(cell), view.axis) < -0.15 {
            continue;
        }
        let corners = cone_layer.vertices(cell);
        let mut path = String::new();
        for edge in 0..4 {
            let (a, b) = (corners[edge], corners[(edge + 1) % 4]);
            for step in 0..4 {
                let (x, y) = view.project(slerp(a, b, step as f64 / 4.0));
                let _ = write!(
                    path,
                    "{}{x:.1} {y:.1} ",
                    if edge == 0 && step == 0 { "M" } else { "L" }
                );
            }
        }
        path.push('Z');
        let lit = dot(cone_layer.center_vec(cell), view.axis).max(0.0);
        let t = 0.45 + 0.55 * lit;
        let _ = writeln!(
            cone_fill,
            r#"    <path d="{path}" fill="rgb({},{},{})"/>"#,
            (138.0 + 107.0 * t) as u32,
            (72.0 + 90.0 * t) as u32,
            (24.0 + 44.0 * t) as u32
        );

        // Outline only where the coverage stops, which is the edge of the search itself.
        for (edge, direction) in EDGE_NEIGHBOUR.iter().enumerate() {
            let outside = match cone_layer.neighbour(cell, *direction) {
                None => true,
                Some(n) => covered.binary_search(&n).is_err(),
            };
            if !outside {
                continue;
            }
            let (a, b) = (corners[edge], corners[(edge + 1) % 4]);
            let mut seg = String::new();
            for step in 0..=4 {
                let (x, y) = view.project(slerp(a, b, step as f64 / 4.0));
                let _ = write!(seg, "{}{x:.1} {y:.1} ", if step == 0 { "M" } else { "L" });
            }
            let _ = writeln!(cone_edge, r#"    <path d="{}"/>"#, seg.trim_end());
        }
    }

    let svg = format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {SIZE:.0} {SIZE:.0}" width="{SIZE:.0}" height="{SIZE:.0}" role="img" aria-label="REALPix">
  <title>REALPix</title>
  <defs>
    <clipPath id="sphere">
      <circle cx="{c:.0}" cy="{c:.0}" r="{R:.0}"/>
    </clipPath>
    <radialGradient id="limb" cx="35%" cy="30%" r="78%">
      <stop offset="0%" stop-color="#1b2a5e"/>
      <stop offset="62%" stop-color="#0d1430"/>
      <stop offset="100%" stop-color="#05070f"/>
    </radialGradient>
    <radialGradient id="glow" cx="34%" cy="28%" r="74%">
      <stop offset="0%" stop-color="#ffffff" stop-opacity="0.18"/>
      <stop offset="46%" stop-color="#ffffff" stop-opacity="0.02"/>
      <stop offset="78%" stop-color="#000010" stop-opacity="0.30"/>
      <stop offset="100%" stop-color="#000008" stop-opacity="0.70"/>
    </radialGradient>
  </defs>

  <circle cx="{c:.0}" cy="{c:.0}" r="{R:.0}" fill="url(#limb)"/>

  <g clip-path="url(#sphere)">
{cells}
    <g fill="none" stroke="#8fbcff" stroke-opacity="0.38" stroke-width="1.1" stroke-linecap="round">
{faint}    </g>
{cone_fill}    <g fill="none" stroke="#ffd9a0" stroke-opacity="0.95" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
{cone_edge}    </g>
    <g fill="none" stroke="#dbeafe" stroke-opacity="0.92" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round">
{bold}    </g>
    <circle cx="{c:.0}" cy="{c:.0}" r="{R:.0}" fill="url(#glow)"/>
  </g>

  <circle cx="{c:.0}" cy="{c:.0}" r="{R:.0}" fill="none" stroke="#9ec5ff" stroke-opacity="0.55" stroke-width="2"/>
</svg>
"##,
        c = SIZE / 2.0
    );

    std::fs::create_dir_all("assets")?;
    std::fs::write("assets/logo.svg", &svg)?;
    println!(
        "assets/logo.svg written ({} bytes, {} cells in the cone)",
        svg.len(),
        covered.len()
    );
    Ok(())
}
