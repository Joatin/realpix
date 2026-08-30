//! Merging of the index ranges emitted by a coverage query.

use core::ops::Range;

/// Merges ranges into the largest possible runs before handing them to the caller.
///
/// Ranges must be pushed in non-decreasing order of `start`: each new range is either
/// contiguous with (or contained in) the pending one, or entirely beyond it. Both cone
/// searches emit their cells in increasing index order, which satisfies this.
pub(crate) struct Merger<F: FnMut(Range<u64>)> {
    pending: Option<Range<u64>>,
    sink: F,
}

impl<F: FnMut(Range<u64>)> Merger<F> {
    #[inline]
    pub(crate) fn new(sink: F) -> Self {
        Self {
            pending: None,
            sink,
        }
    }

    #[inline]
    pub(crate) fn push(&mut self, range: Range<u64>) {
        if let Some(pending) = self.pending.as_mut()
            && range.start <= pending.end
        {
            pending.end = pending.end.max(range.end);
            return;
        }
        if let Some(previous) = self.pending.replace(range) {
            (self.sink)(previous);
        }
    }

    #[inline]
    pub(crate) fn flush(mut self) {
        if let Some(pending) = self.pending.take() {
            (self.sink)(pending);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec::Vec;

    fn merged(input: &[Range<u64>]) -> Vec<Range<u64>> {
        let mut out = Vec::new();
        let mut m = Merger::new(|r| out.push(r));
        for r in input {
            m.push(r.clone());
        }
        m.flush();
        out
    }

    #[test]
    fn adjacent_and_overlapping_runs_are_joined() {
        let one = |r: Range<u64>| alloc::vec![r];
        assert_eq!(merged(&[0..2, 2..5, 5..6]), one(0..6));
        assert_eq!(merged(&[0..4, 1..3]), one(0..4));
        assert_eq!(merged(&[0..2, 3..5]), alloc::vec![0..2, 3..5]);
        assert!(merged(&[]).is_empty());
    }
}
