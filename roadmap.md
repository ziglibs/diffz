# Roadmap

- [✅] Port patch
    - [ ] Add DiffMatchPatch object instead of @This() (which won't work)
- [✅] Port match
- [ ] Port test coverage
- [ ] Diff stream
    - [ ] Use Unicode characters and codepoint indices - 32
    - [ ] Implement line diff as a stream
    - [ ] Also gives word diff, token diff, etc.
- [ ] Refactor:
    - [ ] Diff struct becomes Edit
    - [ ] DiffList stays
    - [ ] New Diff struct, and DiffUnmanaged
        - [ ] Namespaces subsequent operations on diffs
- [ ] Histogram?
    - [ ] Imara diff has an optimized histogram:
          https://github.com/pascalkuthe/imara-diff
- [ ] POSIX-diff compatible patch output?
    - [ ] This one seems pretty worthwhile to me.
- [ ] Delta functions?  They aren't used internally.

Covers the bases.
