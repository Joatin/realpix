//! Keeps the benchmark suite honest: every public function must be benchmarked.
//!
//! `benches/healpix.rs` names each benchmark `module::function`. This test derives the
//! same set of names from the source and fails if one of them is missing, so a new public
//! function cannot be added without a benchmark to go with it.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// The module a source file contributes its public functions to.
fn module_of(path: &Path) -> String {
    let relative = path.strip_prefix(root().join("src")).unwrap();
    let parts: Vec<&str> = relative
        .components()
        .map(|c| c.as_os_str().to_str().unwrap())
        .collect();
    // `nested/mod.rs` and `nested/cone.rs` both belong to `nested`.
    if parts.len() > 1 {
        parts[0].to_string()
    } else {
        parts[0].trim_end_matches(".rs").to_string()
    }
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

/// Every public function in the crate, as `module::function`.
fn public_functions() -> BTreeSet<String> {
    let mut files = Vec::new();
    rust_files(&root().join("src"), &mut files);
    files.sort();

    let mut out = BTreeSet::new();
    for path in &files {
        // `lib.rs` only re-exports; the definitions live in the other files.
        if path.file_name().unwrap() == "lib.rs" {
            continue;
        }
        let module = module_of(path);
        let text = fs::read_to_string(path).unwrap();
        let mut in_test_module = false;
        for line in text.lines() {
            if line.trim_start().starts_with("mod tests") {
                in_test_module = true;
            }
            if in_test_module {
                continue;
            }
            let trimmed = line.trim_start();
            let Some(rest) = trimmed.strip_prefix("pub ") else {
                continue;
            };
            let rest = rest.strip_prefix("const ").unwrap_or(rest);
            let Some(rest) = rest.strip_prefix("fn ") else {
                continue;
            };
            let name: String = rest
                .chars()
                .take_while(|c| c.is_alphanumeric() || *c == '_')
                .collect();
            if !name.is_empty() {
                out.insert(format!("{module}::{name}"));
            }
        }
    }
    out
}

#[test]
fn every_public_function_is_benchmarked() {
    let benches = fs::read_to_string(root().join("benches/healpix.rs")).unwrap();
    let functions = public_functions();

    // A sanity check on the parsing itself: if this ever collapses to a handful, the
    // extraction has broken and the coverage check below would pass vacuously.
    assert!(
        functions.len() > 80,
        "only found {} public functions — the parser is broken",
        functions.len()
    );

    let missing: Vec<&String> = functions
        .iter()
        .filter(|name| !benches.contains(name.as_str()))
        .collect();

    assert!(
        missing.is_empty(),
        "{} public function(s) have no benchmark in benches/healpix.rs:\n{}",
        missing.len(),
        missing
            .iter()
            .map(|m| format!("  {m}"))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// The reverse: a benchmark naming a function that no longer exists is a stale benchmark.
#[test]
fn no_benchmark_names_a_function_that_is_gone() {
    let benches = fs::read_to_string(root().join("benches/healpix.rs")).unwrap();
    let functions = public_functions();

    let mut stale = Vec::new();
    for line in benches.lines() {
        // Benchmark ids appear as string literals containing `::`.
        for chunk in line.split('"').skip(1).step_by(2) {
            let id = chunk.split('/').next().unwrap_or(chunk);
            let id = id.split_whitespace().next().unwrap_or(id);
            if id.contains("::") && !functions.contains(id) {
                stale.push(id.to_string());
            }
        }
    }
    stale.sort();
    stale.dedup();
    assert!(
        stale.is_empty(),
        "benchmark id(s) name a function that is not public any more:\n{}",
        stale
            .iter()
            .map(|s| format!("  {s}"))
            .collect::<Vec<_>>()
            .join("\n")
    );
}
