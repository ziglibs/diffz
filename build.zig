const std = @import("std");

pub fn build(b: *std.Build) void {
    _ = b.addModule("diffz", .{
        .root_source_file = .{ .path = "DiffMatchPatch.zig" },
    });
}
