const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    _ = b.addModule("diffz", .{
        .root_source_file = b.path("DiffMatchPatch.zig"),
        .target = target,
        .optimize = optimize,
    });

    const lib = b.addStaticLibrary(.{
        .name = "diffz",
        .root_source_file = b.path("DiffMatchPatch.zig"),
        .target = target,
        .optimize = optimize,
    });

    // This declares intent for the library to be installed into the standard
    // location when the user invokes the "install" step (the default step when
    // running `zig build`).
    b.installArtifact(lib);

    // Run tests
    const tests = b.addTest(.{
        .name = "tests",
        .root_source_file = b.path("DiffMatchPatch.zig"),
        .target = target,
        .optimize = optimize,
    });
    const step_tests = b.addRunArtifact(tests);
    step_tests.has_side_effects = true;

    b.step("test", "Run diffz tests").dependOn(&step_tests.step);

    // Adds a step to generate code coverage
    const cov_step = b.step("cov", "Generate coverage (kcov must be installed)");

    const cov_run = b.addSystemCommand(&.{
        "kcov",
        "--clean",
        "--include-pattern=DiffMatchPatch.zig",
        "--exclude-line=unreachable,expect(false)",
        "kcov-output",
    });
    cov_run.addArtifactArg(tests);
    cov_step.dependOn(&cov_run.step);
    _ = cov_run.captureStdOut();
    _ = cov_run.captureStdErr();
}
