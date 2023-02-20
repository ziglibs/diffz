const std = @import("std");

pub fn build(b: *std.Build) void {
    b.addModule(.{
        .name = "diffz",
        .source_file = .{ .path = "DiffMatchPatch.zig" },
    });
}
