# diffz

An implementation of Google's diff-match-patch.

Currently implemented:

- [x] Diff
- [ ] Match
- [ ] Patch

## Installation

> [!NOTE]
> The minimum supported Zig version is `0.14.0`.

```bash
# Initialize a `zig build` project if you haven't already
zig init
# Add the `diffz` package to your `build.zig.zon`
zig fetch --save git+https://github.com/ziglibs/diffz.git
```

You can then import `diffz` in your `build.zig` with:

```zig
const diffz = b.dependency("diffz", .{});
const exe = b.addExecutable(...);
exe.root_module.addImport("diffz", diffz.module("diffz"));
```

## License

This library is based off of https://github.com/google/diff-match-patch, which is licensed under the [Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0). This library itself is licensed under the MIT License, see `LICENSE`.
