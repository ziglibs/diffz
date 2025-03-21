const std = @import("std");
const builtin = @import("builtin");

const minimum_zig_version = std.SemanticVersion.parse("0.14.0") catch unreachable;

pub fn build(b: *std.Build) void {
    if (comptime (builtin.zig_version.order(minimum_zig_version) == .lt)) {
        @compileError(std.fmt.comptimePrint(
            \\Your Zig version does not meet the minimum build requirement:
            \\  required Zig version: {[minimum_zig_version]}
            \\  actual   Zig version: {[current_version]}
            \\
        , .{
            .current_version = builtin.zig_version,
            .minimum_zig_version = minimum_zig_version,
        }));
    }

    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const diffz_module = b.addModule("diffz", .{
        .root_source_file = b.path("DiffMatchPatch.zig"),
        .target = target,
        .optimize = optimize,
    });

    const tests = b.addTest(.{ .root_module = diffz_module });
    const run_tests = b.addRunArtifact(tests);

    const test_step = b.step("test", "Run all the tests");
    test_step.dependOn(&run_tests.step);

    const kcov_bin = b.findProgram(&.{"kcov"}, &.{}) catch "kcov";

    const run_kcov = b.addSystemCommand(&.{
        kcov_bin,
        "--clean",
        "--exclude-line=unreachable,expect(false)",
    });
    run_kcov.addPrefixedDirectoryArg("--include-pattern=", b.path("."));
    const coverage_output = run_kcov.addOutputDirectoryArg(".");
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
