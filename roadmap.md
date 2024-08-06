# Roadmap

- [✅] Port patch
    - [✅] Add DiffMatchPatch object instead of @This() (which won't work)
- [✅] Port match.
- [✅] Port test coverage.
- [✅] Unicode-aware `diffLineMode`.
    - [✅] Coverage for all corner cases of preventing diff splits which aren't
          on valid UTF-8 boundaries.
    - [✅] Convert `line_array` to encode UTF-8 byte sequences and store `u21` keys
    - [✅] Make the inner function accept a stream iterator, one which delivers the
          entire string with boundaries (where applicable) at the end.
- [ ] Refactor: the port currently treats Diffs and Patches as raw ArrayLists,
      these should be proper Zig objects, with member functions, and probably
      come in an Unmanaged and normal form.
    - [?] Diff struct becomes Edit.  Patch also needs a name, because Diff and
          Patch should be the names of the user-facing structs.  The name for
          what a patch is in classic diff/patch programs is Hunk, so that's a
          justifiable choice.
    - [ ] DiffList and PatchList remain same, used internally.
    - [ ] New Diff struct, and DiffUnmanaged.
        - [ ] Namespaces subsequent operations on diffs.
    - [ ] Same for Patch and PatchUnmanaged.  These are little more than the
          relevant DiffList, a DiffMatchPatch instance, and some decl functions,
          plus the Allocator for managed versions.
- [ ] Enhancements
    - [ ] Extend Bitap algorithm to handle larger patches. The algorithm takes
          `m * n` space, where `m` is unique bytes in the pattern and `n` is the
          pattern length, so I think the idea of doing it up to 2048 bytes/bits
          was optimistic on my part.  But comptime-gated function specializations
          for 64 (status quo), 128, and 256 bytes, would mean a lot less frobbing
          and munging the patches internally.  Performance implications expected
          to be positive, if not hugely so.  The algorithm is also amenable to
          SIMD acceleration, although I'm not going to do that.
    - [ ] `diffsForRegion`: provides every diff pertaining to a specific
          region of `before`.  Needs to also include how much overlap, if
          any, the diff includes.  Should have "borrow" and "copy"
          versions.  Signature being `diffsForRegion(diffs: DiffList, start: usize,`
          `end: usize, <borrow or copy enum>) ?DiffList`.
    - [ ] Implement a delta format which doesn't suck so badly.  I have copious
          notes on this.
    - [ ] I'd also like to break compatibility with the 'Unidiff' format, in a
          less dramatic way.  It's mostly the compulsive percent-encoding of
          everything which doesn't fit in a URI, it's Googley [derogatory] and
          a UTF-8 native patch format has no need for this.  This would be a
          separate, sucks-less text format, differentiable by its header,
          decodes into a Patch in the same basic way.  The legacy form is
          already ported and should be kept.
    - [ ] Add `Diff.differs() bool`, which checks if there are any differences
          between the before and after text.
    - [✅] Diff stream
        - [✅] Use Unicode characters and codepoint indices - 32.
        - [✅] Implement line diff as a stream.
        - [✅] Also gives word diff, token diff, etc.
- [ ] Histogram?
    - [ ] Imara diff has an optimized histogram:
          https://github.com/pascalkuthe/imara-diff
    - [ ] Calculating the histogram while hashing the lines would be
          straightforward, this could be comptime-gated, but probably
          just a second copy of the munge function is fine.
    - [ ] This one is getting into overkill territory perhaps.
- [ ] POSIX-diff compatible patch output?
    - [ ] This one seems pretty worthwhile to me.  It would need to call line
          mode without refining further, but everything else is fairly simple.
- [ ] Delta functions?  They aren't used internally.  I favor ignoring the
      legacy version and implementing a better one.

Covers the bases.
