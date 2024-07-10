# Roadmap

- [✅] Port patch
    - [✅] Add DiffMatchPatch object instead of @This() (which won't work)
- [✅] Port match.
- [✅] Port test coverage.
- [ ] Unicode-aware `diffLineMode`.
    - [ ] Coverage for all corner cases of preventing diff splits which aren't
          on valid UTF-8 boundaries.
    - [ ] Convert `line_array` to encode UTF-8 byte sequences and store `u21` keys
    - [ ] Make the inner function accept a stream iterator, one which delivers the
          entire string with boundaries (where applicable) at the end.
- [ ] Refactor: the port currently treats Diffs and Patches as raw ArrayLists,
      these should be proper Zig objects, with member functions, and probably
      come in an Unmanaged and normal form.
    - [?] Diff struct becomes Edit.
    - [ ] DiffList and PatchList remain same, used internally.
    - [ ] New Diff struct, and DiffUnmanaged.
        - [ ] Namespaces subsequent operations on diffs.
- [ ] Enhancements
    - [ ] Extend Bitap algorithm to handle larger patches. The algorithm takes
          `m * n` space, where `m` is unique bytes in the pattern and `n` is the
          pattern length, so I think the idea of doing it up to 2048 bytes/bits
          was optimistic on my part.
    - [ ] `diffsForRegion`: provides every diff pertaining to a specific
          region of `before`.  Needs to also include how much overlap, if
          any, the diff includes.  Should have "borrow" and "copy"
          versions.
    - [ ] Implement a delta function which doesn't suck so badly.
    - [ ] Diff stream
        - [ ] Use Unicode characters and codepoint indices - 32.
        - [ ] Implement line diff as a stream.
        - [ ] Also gives word diff, token diff, etc.
- [ ] Histogram?
    - [ ] Imara diff has an optimized histogram:
          https://github.com/pascalkuthe/imara-diff
- [ ] POSIX-diff compatible patch output?
    - [ ] This one seems pretty worthwhile to me.  It would need to call line
          mode without refining further, but everything else is fairly simple.
- [ ] Delta functions?  They aren't used internally.  I favor ignoring the
      legacy version and implementing a better one.

Covers the bases.
