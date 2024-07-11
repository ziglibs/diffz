const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    _ = b.addModule("diffz", .{
        .root_source_file = b.path("DiffMatchPatch.zig"),
        .target = target,
        .optimize = optimize,
    });

    // Run tests
    const tests = b.addTest(.{
        .name = "tests",
        .root_source_file = b.path("DiffMatchPatch.zig"),
        .target = target,
        .optimize = optimize,
    });
    const step_tests = b.addRunArtifact(tests);

    b.step("test", "Run diffz tests").dependOn(&step_tests.step);

    const addOutputDirectoryArg = comptime if (@import("builtin").zig_version.order(.{ .major = 0, .minor = 13, .patch = 0 }) == .lt)
        std.Build.Step.Run.addOutputFileArg
    else
        std.Build.Step.Run.addOutputDirectoryArg;

    const run_kcov = b.addSystemCommand(&.{
        "kcov",
        "--clean",
        "--exclude-line=unreachable,expect(false)",
    });
    run_kcov.addPrefixedDirectoryArg("--include-pattern=", b.path("."));
    const coverage_output = addOutputDirectoryArg(run_kcov, ".");
    run_kcov.addArtifactArg(tests);
    run_kcov.enableTestRunnerMode();

    const install_coverage = b.addInstallDirectory(.{
        .source_dir = coverage_output,
        .install_dir = .{ .custom = "coverage" },
        .install_subdir = "",
    });

    const coverage_step = b.step("coverage", "Generate coverage (kcov must be installed)");
    coverage_step.dependOn(&install_coverage.step);
}
