//! Measures the largest centre-to-boundary angular distance over a whole layer.
//! Used to derive `nested::MAX_CENTER_TO_VERTEX`; not part of the public API.
use realpix::nested;

fn ang(a: [f64; 3], b: [f64; 3]) -> f64 {
    let c = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
        .sqrt()
        .atan2(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
}

fn main() {
    println!("depth  nside  max_vertex        nside*max   max_edge(sampled)  nside*max");
    for depth in 0..=9u8 {
        let layer = nested::get(depth);
        let nside = layer.nside() as f64;
        let (mut mv, mut me) = (0.0f64, 0.0f64);
        for cell in layer.iter() {
            let c = layer.center_vec(cell);
            for v in layer.vertices(cell) {
                mv = mv.max(ang(c, v));
            }
            // Sample the four edges as well: HEALPix cell edges are not great circles, so
            // a bound built only from corners could in principle be too small.
            if depth <= 6 {
                let vs = layer.vertices(cell);
                for k in 0..4 {
                    let (a, b) = (vs[k], vs[(k + 1) % 4]);
                    for s in 1..8 {
                        let t = s as f64 / 8.0;
                        let p = [
                            a[0] + (b[0] - a[0]) * t,
                            a[1] + (b[1] - a[1]) * t,
                            a[2] + (b[2] - a[2]) * t,
                        ];
                        let n = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                        me = me.max(ang(c, [p[0] / n, p[1] / n, p[2] / n]));
                    }
                }
            }
        }
        println!(
            "{depth:5}  {:5}  {mv:.15}  {:.6}   {me:.15}  {:.6}",
            layer.nside(),
            mv * nside,
            me * nside
        );
    }
}
