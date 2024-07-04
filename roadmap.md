# Roadmap

- [ ] Port patch
- [ ] Port match
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
- [ ] Delta functions?  They aren't used internally.

Covers the bases.
