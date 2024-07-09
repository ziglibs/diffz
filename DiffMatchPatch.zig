const DiffMatchPatch = @This();

const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const Allocator = std.mem.Allocator;
const ArrayListUnmanaged = std.ArrayListUnmanaged;
const DiffList = ArrayListUnmanaged(Diff);
const PatchList = ArrayListUnmanaged(Patch);

pub const DiffError = error{
    OutOfMemory,
    BadPatchString,
};

//| XXX This boolean is entirely for calming the compiler down while working

const XXX = false;

//| Fields

/// Number of milliseconds to map a diff before giving up (0 for infinity).
diff_timeout: u64 = 1000,
/// Cost of an empty edit operation in terms of edit characters.
diff_edit_cost: u16 = 4,

/// At what point is no match declared (0.0 = perfection, 1.0 = very loose).
/// This defaults to 0.05, on the premise that the library will mostly be
/// used in cases where failure is better than a bad patch application.
match_threshold: f64 = 0.05,

/// How far to search for a match (0 = exact location, 1000+ = broad match).
/// A match this many characters away from the expected location will add
/// 1.0 to the score (0.0 is a perfect match).
match_distance: u32 = 1000,

/// The number of bits in a usize.
match_max_bits: u8 = @bitSizeOf(usize),

/// When deleting a large block of text (over ~64 characters), how close
/// do the contents have to be to match the expected contents. (0.0 =
/// perfection, 1.0 = very loose).  Note that Match_Threshold controls
/// how closely the end points of a delete need to match.
patch_delete_threshold: f32 = 0.5,

/// Chunk size for context length.
patch_margin: u8 = 4,

//| Allocation Management Helpers

/// Deinit an `ArrayListUnmanaged(Diff)` and the allocated slices of
/// text in each `Diff`.
pub fn deinitDiffList(allocator: Allocator, diffs: *DiffList) void {
    defer diffs.deinit(allocator);
    for (diffs.items) |d| {
        allocator.free(d.text);
    }
}

/// Free a range of Diffs inside a list.  Used during cleanups and
/// edits.
fn freeRangeDiffList(
    allocator: Allocator,
    diffs: *DiffList,
    start: usize,
    len: usize,
) void {
    const after_range = start + len;
    const range = diffs.items[start..after_range];
    for (range) |d| {
        allocator.free(d.text);
    }
}

pub fn deinitPatchList(allocator: Allocator, patches: *PatchList) void {
    defer patches.deinit(allocator);
    for (patches.items) |*a_patch| {
        deinitDiffList(allocator, &a_patch.diffs);
    }
}

/// DMP with default configuration options
pub const default = DiffMatchPatch{};

/// Represents a single edit operation.
/// TODO rename this Edit
pub const Diff = struct {
    pub const Operation = enum {
        insert,
        delete,
        equal,
    };

    operation: Operation,
    text: []const u8,

    pub fn format(value: Diff, _: anytype, _: anytype, writer: anytype) !void {
        try writer.print("({s}, \"{s}\")", .{
            switch (value.operation) {
                .equal => "=",
                .insert => "+",
                .delete => "-",
            },
            value.text,
        });
    }

    pub fn init(operation: Operation, text: []const u8) Diff {
        return .{ .operation = operation, .text = text };
    }
    pub fn eql(a: Diff, b: Diff) bool {
        return a.operation == b.operation and std.mem.eql(u8, a.text, b.text);
    }

    pub fn clone(self: Diff, allocator: Allocator) !Diff {
        return Diff{
            .operation = self.operation,
            .text = try allocator.dupe(u8, self.text),
        };
    }
};

pub const Patch = struct {
    /// Diffs to be applied
    diffs: DiffList = DiffList{},
    /// Start of patch in before text
    start1: usize = 0,
    length1: usize = 0,
    /// Start of patch in after text
    start2: usize = 0,
    length2: usize = 0,

    /// Make a clone of the Patch, including all Diffs.
    pub fn clone(patch: Patch, allocator: Allocator) !Patch {
        var new_diffs = DiffList{};
        new_diffs.initCapacity(allocator, patch.diffs.items.len);
        for (patch.diffs) |a_diff| {
            try new_diffs.append(try a_diff.clone(allocator));
        }
        return Patch{
            .diffs = new_diffs,
            .start1 = patch.start1,
            .length1 = patch.length1,
            .start2 = patch.start2,
            .length2 = patch.length2,
        };
    }

    pub fn deinit(patch: *Patch, allocator: Allocator) void {
        deinitDiffList(allocator, &patch.diffs);
    }

    /// Emit patch in Unidiff format, as specifified here:
    /// https://github.com/google/diff-match-patch/wiki/Unidiff
    /// This is similar to GNU Unidiff format, but not identical.
    /// Header: @@ -382,8 +481,9 @@
    /// Indices are printed as 1-based, not 0-based.
    /// @return The GNU diff string.
    pub fn asText(patch: Patch, allocator: Allocator) ![]const u8 {
        var text_array = std.ArrayList(u8).init(allocator);
        defer text_array.deinit();
        const writer = text_array.writer();
        try patch.writeText(writer);
        return text_array.toOwnedSlice();
    }

    const format = std.fmt.format;

    /// Stream textual patch representation to Writer.  See `asText`
    /// for more information.
    pub fn writeText(patch: Patch, writer: anytype) !void {
        // Write header.
        _ = try writer.write(PATCH_HEAD);
        // Stream coordinates
        if (patch.length1 == 0) {
            try format(writer, "{d},0", .{patch.start1});
        } else if (patch.length1 == 1) {
            try format(writer, "{d}", .{patch.start1 + 1});
        } else {
            try format(writer, "{d},{d}", .{ patch.start1 + 1, patch.length1 });
        }
        _ = try writer.write(" +");
        if (patch.length2 == 0) {
            try std.fmt.format(writer, "{d},0", .{patch.start2});
        } else if (patch.length2 == 1) {
            _ = try format(writer, "{d}", .{patch.start2 + 1});
        } else {
            try format(writer, "{d},{d}", .{ patch.start2 + 1, patch.length2 });
        }
        _ = try writer.write(PATCH_TAIL);
        // Escape the body of the patch with %xx notation.
        for (patch.diffs.items) |a_diff| {
            switch (a_diff.operation) {
                .insert => try writer.writeByte('+'),
                .delete => try writer.writeByte('-'),
                .equal => try writer.writeByte(' '),
            }
            _ = try writeUriEncoded(writer, a_diff.text);
            try writer.writeByte('\n');
        }
        return;
    }
};

const PATCH_HEAD = "@@ -";
const PATCH_TAIL = " @@\n";

/// Find the differences between two texts.
/// @param before Old string to be diffed.
/// @param after New string to be diffed.
/// @param checklines Speedup flag.  If false, then don't run a
///     line-level diff first to identify the changed areas.
///     If true, then run a faster slightly less optimal diff.
/// @return List of Diff objects.
pub fn diff(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    /// If false, then don't run a line-level diff first
    /// to identify the changed areas. If true, then run
    /// a faster slightly less optimal diff.
    check_lines: bool,
) DiffError!DiffList {
    const deadline = if (dmp.diff_timeout == 0)
        std.math.maxInt(u64)
    else
        @as(u64, @intCast(std.time.milliTimestamp())) + dmp.diff_timeout;
    return dmp.diffInternal(allocator, before, after, check_lines, deadline);
}

fn diffInternal(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    check_lines: bool,
    deadline: u64,
) DiffError!DiffList {
    // Check for equality (speedup).
    if (std.mem.eql(u8, before, after)) {
        var diffs = DiffList{};
        errdefer deinitDiffList(allocator, &diffs);
        if (before.len != 0) {
            try diffs.ensureUnusedCapacity(allocator, 1);
            diffs.appendAssumeCapacity(Diff.init(
                .equal,
                try allocator.dupe(u8, before),
            ));
        }
        return diffs;
    }

    // Trim off common prefix (speedup).
    var common_length = diffCommonPrefix(before, after);
    const common_prefix = before[0..common_length];
    var trimmed_before = before[common_length..];
    var trimmed_after = after[common_length..];

    // Trim off common suffix (speedup).
    common_length = diffCommonSuffix(trimmed_before, trimmed_after);
    const common_suffix = trimmed_before[trimmed_before.len - common_length ..];
    trimmed_before = trimmed_before[0 .. trimmed_before.len - common_length];
    trimmed_after = trimmed_after[0 .. trimmed_after.len - common_length];

    // Compute the diff on the middle block.
    var diffs = try dmp.diffCompute(allocator, trimmed_before, trimmed_after, check_lines, deadline);
    errdefer deinitDiffList(allocator, &diffs);

    // Restore the prefix and suffix.
    if (common_prefix.len != 0) {
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.insertAssumeCapacity(0, Diff.init(
            .equal,
            try allocator.dupe(u8, common_prefix),
        ));
    }
    if (common_suffix.len != 0) {
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.appendAssumeCapacity(Diff.init(
            .equal,
            try allocator.dupe(u8, common_suffix),
        ));
    }

    try diffCleanupMerge(allocator, &diffs);
    return diffs;
}

/// Test if a byte is a UTF-8 follow byte
inline fn is_follow(byte: u8) bool {
    return byte & 0b1100_0000 == 0b1000_0000;
}

/// Find a common prefix which respects UTF-8 code point boundaries.
fn diffCommonPrefix(before: []const u8, after: []const u8) usize {
    const n = @min(before.len, after.len);
    var i: usize = 0;
    while (i < n) : (i += 1) {
        var b = before[i];
        const a = after[i];
        if (a != b) {
            if (is_follow(a) and is_follow(b)) {
                // We've clipped a codepoint, back out
                if (i == 0) return i; // Malformed UTF-8 is always possible
                i -= 1;
                // We'll track `before` since they must be the same:
                b = before[i];
                assert(b == after[i]);
                while (i != 0 and is_follow(b)) {
                    i -= 1;
                    b = before[i];
                    assert(b == after[i]);
                }
                // Now we're either at zero, or at the lead:
                return i;
            } else {
                return i;
            }
        }
    }

    return n;
}

/// Find a common suffix which respects UTF-8 code point boundaries
fn diffCommonSuffix(before: []const u8, after: []const u8) usize {
    const n = @min(before.len, after.len);
    var i: usize = 1;
    var was_follow = false;
    while (i <= n) : (i += 1) {
        var b = before[before.len - i];
        const a = after[after.len - i];
        if (a != b) {
            if (was_follow) {
                // Means we're at at least 2:
                assert(i > 1);
                // We just saw an identical follow byte, so we back
                // out forward:
                i -= 1;
                b = before[before.len - i];
                assert(b == after[after.len - i]);
                while (i > 1 and is_follow(b)) {
                    i -= 1;
                    b = before[before.len - i];
                    assert(b == after[after.len - i]);
                } // Either at one, or no more follow bytes:
                return i - 1;
            } else {
                return i - 1;
            }
        } else {
            was_follow = is_follow(b); // no need to check twice
        }
    }

    return n;
}

/// Find the differences between two texts.  Assumes that the texts do not
/// have any common prefix or suffix.
/// @param before Old string to be diffed.
/// @param after New string to be diffed.
/// @param checklines Speedup flag.  If false, then don't run a
///     line-level diff first to identify the changed areas.
///     If true, then run a faster slightly less optimal diff.
/// @param deadline Time when the diff should be complete by.
/// @return List of Diff objects.
fn diffCompute(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    check_lines: bool,
    deadline: u64,
) DiffError!DiffList {
    if (before.len == 0) {
        // Just add some text (speedup).
        var diffs = DiffList{};
        errdefer deinitDiffList(allocator, &diffs);
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.appendAssumeCapacity(Diff.init(
            .insert,
            try allocator.dupe(u8, after),
        ));
        return diffs;
    }

    if (after.len == 0) {
        // Just delete some text (speedup).
        var diffs = DiffList{};
        errdefer deinitDiffList(allocator, &diffs);
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.appendAssumeCapacity(Diff.init(
            .delete,
            try allocator.dupe(u8, before),
        ));
        return diffs;
    }

    const long_text = if (before.len > after.len) before else after;
    const short_text = if (before.len > after.len) after else before;

    if (std.mem.indexOf(u8, long_text, short_text)) |index| {
        // Shorter text is inside the longer text (speedup).
        var diffs = DiffList{};
        errdefer deinitDiffList(allocator, &diffs);
        const op: Diff.Operation = if (before.len > after.len)
            .delete
        else
            .insert;
        try diffs.ensureUnusedCapacity(allocator, 3);
        diffs.appendAssumeCapacity(Diff.init(
            op,
            try allocator.dupe(u8, long_text[0..index]),
        ));
        diffs.appendAssumeCapacity(Diff.init(
            .equal,
            try allocator.dupe(u8, short_text),
        ));
        diffs.appendAssumeCapacity(Diff.init(
            op,
            try allocator.dupe(u8, long_text[index + short_text.len ..]),
        ));
        return diffs;
    }

    if (short_text.len == 1) {
        // Single character string.
        // After the previous speedup, the character can't be an equality.
        var diffs = DiffList{};
        errdefer deinitDiffList(allocator, &diffs);
        try diffs.ensureUnusedCapacity(allocator, 2);
        diffs.appendAssumeCapacity(Diff.init(
            .delete,
            try allocator.dupe(u8, before),
        ));
        diffs.appendAssumeCapacity(Diff.init(
            .insert,
            try allocator.dupe(u8, after),
        ));
        return diffs;
    }

    // Check to see if the problem can be split in two.
    if (try dmp.diffHalfMatch(allocator, before, after)) |half_match| {
        // A half-match was found, sort out the return data.
        defer half_match.deinit(allocator);
        // Send both pairs off for separate processing.
        var diffs = try dmp.diffInternal(
            allocator,
            half_match.prefix_before,
            half_match.prefix_after,
            check_lines,
            deadline,
        );
        errdefer deinitDiffList(allocator, &diffs);
        var diffs_b = try dmp.diffInternal(
            allocator,
            half_match.suffix_before,
            half_match.suffix_after,
            check_lines,
            deadline,
        );
        defer diffs_b.deinit(allocator);
        // we have to deinit regardless, so deinitDiffList would be
        // a double free:
        errdefer {
            for (diffs_b.items) |d| {
                allocator.free(d.text);
            }
        }

        // Merge the results.
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.appendAssumeCapacity(
            Diff.init(.equal, try allocator.dupe(
                u8,
                half_match.common_middle,
            )),
        );
        try diffs.appendSlice(allocator, diffs_b.items);
        return diffs;
    }

    if (check_lines and before.len > 100 and after.len > 100) {
        return dmp.diffLineMode(allocator, before, after, deadline);
    }
    return dmp.diffBisect(allocator, before, after, deadline);
}

const HalfMatchResult = struct {
    prefix_before: []const u8,
    suffix_before: []const u8,
    prefix_after: []const u8,
    suffix_after: []const u8,
    common_middle: []const u8,

    // TODO maybe check for empty slice here for fewer copies,
    // as in, maybe we can transfer ownership and replace with "".
    pub fn deinit(hmr: HalfMatchResult, alloc: Allocator) void {
        alloc.free(hmr.prefix_before);
        alloc.free(hmr.suffix_before);
        alloc.free(hmr.prefix_after);
        alloc.free(hmr.suffix_after);
        alloc.free(hmr.common_middle);
    }
};

/// Do the two texts share a Substring which is at least half the length of
/// the longer text?
/// This speedup can produce non-minimal diffs.
/// @param before First string.
/// @param after Second string.
/// @return Five element String array, containing the prefix of text1, the
///     suffix of text1, the prefix of text2, the suffix of text2 and the
///     common middle.  Or null if there was no match.
fn diffHalfMatch(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
) DiffError!?HalfMatchResult {
    if (dmp.diff_timeout == 0) {
        // Don't risk returning a non-optimal diff if we have unlimited time.
        return null;
    }
    const long_text = if (before.len > after.len) before else after;
    const short_text = if (before.len > after.len) after else before;

    if (long_text.len < 4 or short_text.len * 2 < long_text.len) {
        return null; // Pointless.
    }

    // First check if the second quarter is the seed for a half-match.
    const half_match_1 = try dmp.diffHalfMatchInternal(allocator, long_text, short_text, (long_text.len + 3) / 4);
    errdefer {
        if (half_match_1) |h_m| h_m.deinit(allocator);
    }
    // Check again based on the third quarter.
    const half_match_2 = try dmp.diffHalfMatchInternal(allocator, long_text, short_text, (long_text.len + 1) / 2);
    errdefer {
        if (half_match_2) |h_m| h_m.deinit(allocator);
    }

    var half_match: ?HalfMatchResult = null;
    if (half_match_1 == null and half_match_2 == null) {
        return null;
    } else if (half_match_2 == null) {
        half_match = half_match_1.?;
    } else if (half_match_1 == null) {
        half_match = half_match_2.?;
    } else {
        // Both matched. Select the longest.
        half_match = half: {
            if (half_match_1.?.common_middle.len > half_match_2.?.common_middle.len) {
                half_match_2.?.deinit(allocator);
                break :half half_match_1;
            } else {
                half_match_1.?.deinit(allocator);
                break :half half_match_2;
            }
        };
    }

    // A half-match was found, sort out the return data.
    if (before.len > after.len) {
        return half_match.?;
    } else {
        // Transfers ownership of all memory to new, permuted, half_match.
        const half_match_yes = half_match.?;
        return .{
            .prefix_before = half_match_yes.prefix_after,
            .suffix_before = half_match_yes.suffix_after,
            .prefix_after = half_match_yes.prefix_before,
            .suffix_after = half_match_yes.suffix_before,
            .common_middle = half_match_yes.common_middle,
        };
    }
}

/// Does a Substring of shorttext exist within longtext such that the
/// Substring is at least half the length of longtext?
/// @param longtext Longer string.
/// @param shorttext Shorter string.
/// @param i Start index of quarter length Substring within longtext.
/// @return Five element string array, containing the prefix of longtext, the
///     suffix of longtext, the prefix of shorttext, the suffix of shorttext
///     and the common middle.  Or null if there was no match.
fn diffHalfMatchInternal(
    _: DiffMatchPatch,
    allocator: std.mem.Allocator,
    long_text: []const u8,
    short_text: []const u8,
    i: usize,
) DiffError!?HalfMatchResult {
    // Start with a 1/4 length Substring at position i as a seed.
    const seed = long_text[i .. i + long_text.len / 4];
    var j: isize = -1;

    var best_common = std.ArrayListUnmanaged(u8){};
    defer best_common.deinit(allocator);
    var best_long_text_a: []const u8 = "";
    var best_long_text_b: []const u8 = "";
    var best_short_text_a: []const u8 = "";
    var best_short_text_b: []const u8 = "";

    while (j < short_text.len and b: {
        j = @as(isize, @intCast(std.mem.indexOf(u8, short_text[@as(usize, @intCast(j + 1))..], seed) orelse break :b false)) + j + 1;
        break :b true;
    }) {
        const prefix_length = diffCommonPrefix(long_text[i..], short_text[@as(usize, @intCast(j))..]);
        const suffix_length = diffCommonSuffix(long_text[0..i], short_text[0..@as(usize, @intCast(j))]);
        if (best_common.items.len < suffix_length + prefix_length) {
            best_common.items.len = 0;
            const a = short_text[@as(usize, @intCast(j - @as(isize, @intCast(suffix_length)))) .. @as(usize, @intCast(j - @as(isize, @intCast(suffix_length)))) + suffix_length];
            try best_common.appendSlice(allocator, a);
            const b = short_text[@as(usize, @intCast(j)) .. @as(usize, @intCast(j)) + prefix_length];
            try best_common.appendSlice(allocator, b);

            best_long_text_a = long_text[0 .. i - suffix_length];
            best_long_text_b = long_text[i + prefix_length ..];
            best_short_text_a = short_text[0..@as(usize, @intCast(j - @as(isize, @intCast(suffix_length))))];
            best_short_text_b = short_text[@as(usize, @intCast(j + @as(isize, @intCast(prefix_length))))..];
        }
    }
    if (best_common.items.len * 2 >= long_text.len) {
        const prefix_before = try allocator.dupe(u8, best_long_text_a);
        errdefer allocator.free(prefix_before);
        const suffix_before = try allocator.dupe(u8, best_long_text_b);
        errdefer allocator.free(suffix_before);
        const prefix_after = try allocator.dupe(u8, best_short_text_a);
        errdefer allocator.free(prefix_after);
        const suffix_after = try allocator.dupe(u8, best_short_text_b);
        const best_common_text = try best_common.toOwnedSlice(allocator);
        return .{
            .prefix_before = prefix_before,
            .suffix_before = suffix_before,
            .prefix_after = prefix_after,
            .suffix_after = suffix_after,
            .common_middle = best_common_text,
        };
    } else {
        return null;
    }
}

/// Find the 'middle snake' of a diff, split the problem in two
/// and return the recursively constructed diff.
/// See Myers 1986 paper: An O(ND) Difference Algorithm and Its Variations.
/// @param before Old string to be diffed.
/// @param after New string to be diffed.
/// @param deadline Time at which to bail if not yet complete.
/// @return List of Diff objects.
fn diffBisect(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    deadline: u64,
) DiffError!DiffList {
    const before_length: isize = @intCast(before.len);
    const after_length: isize = @intCast(after.len);
    const max_d: isize = @intCast((before.len + after.len + 1) / 2);
    const v_offset = max_d;
    const v_length = 2 * max_d;

    var v1 = try ArrayListUnmanaged(isize).initCapacity(allocator, @as(usize, @intCast(v_length)));
    defer v1.deinit(allocator);
    v1.items.len = @intCast(v_length);
    var v2 = try ArrayListUnmanaged(isize).initCapacity(allocator, @as(usize, @intCast(v_length)));
    defer v2.deinit(allocator);
    v2.items.len = @intCast(v_length);

    var x: usize = 0;
    while (x < v_length) : (x += 1) {
        v1.items[x] = -1;
        v2.items[x] = -1;
    }
    v1.items[@intCast(v_offset + 1)] = 0;
    v2.items[@intCast(v_offset + 1)] = 0;
    const delta = before_length - after_length;
    // If the total number of characters is odd, then the front path will
    // collide with the reverse path.
    const front = (@mod(delta, 2) != 0);
    // Offsets for start and end of k loop.
    // Prevents mapping of space beyond the grid.
    var k1start: isize = 0;
    var k1end: isize = 0;
    var k2start: isize = 0;
    var k2end: isize = 0;

    var d: isize = 0;
    while (d < max_d) : (d += 1) {
        // Bail out if deadline is reached.
        if (@as(u64, @intCast(std.time.milliTimestamp())) > deadline) {
            break;
        }

        // Walk the front path one step.
        var k1 = -d + k1start;
        while (k1 <= d - k1end) : (k1 += 2) {
            const k1_offset = v_offset + k1;
            var x1: isize = 0;
            if (k1 == -d or (k1 != d and
                v1.items[@intCast(k1_offset - 1)] < v1.items[@intCast(k1_offset + 1)]))
            {
                x1 = v1.items[@intCast(k1_offset + 1)];
            } else {
                x1 = v1.items[@intCast(k1_offset - 1)] + 1;
            }
            var y1 = x1 - k1;
            while (x1 < before_length and y1 < after_length) {
                const match, const d1 = equalForward(before, after, x1, y1);
                if (match) {
                    x1 += d1;
                    y1 += d1;
                } else {
                    break;
                }
            }
            v1.items[@intCast(k1_offset)] = x1;
            if (x1 > before_length) {
                // Ran off the right of the graph.
                k1end += 2;
            } else if (y1 > after_length) {
                // Ran off the bottom of the graph.
                k1start += 2;
            } else if (front) {
                const k2_offset = v_offset + delta - k1;
                if (k2_offset >= 0 and k2_offset < v_length and v2.items[@intCast(k2_offset)] != -1) {
                    // Mirror x2 onto top-left coordinate system.
                    const x2 = before_length - v2.items[@intCast(k2_offset)];
                    if (x1 >= x2) {
                        // Overlap detected.
                        return dmp.diffBisectSplit(allocator, before, after, x1, y1, deadline);
                    }
                }
            }
        }

        // Walk the reverse path one step.
        var k2: isize = -d + k2start;
        while (k2 <= d - k2end) : (k2 += 2) {
            const k2_offset = v_offset + k2;
            var x2: isize = 0;
            if (k2 == -d or (k2 != d and
                v2.items[@intCast(k2_offset - 1)] < v2.items[@intCast(k2_offset + 1)]))
            {
                x2 = v2.items[@intCast(k2_offset + 1)];
            } else {
                x2 = v2.items[@intCast(k2_offset - 1)] + 1;
            }
            var y2: isize = x2 - k2;
            while (x2 < before_length and y2 < after_length) {
                const match, const d1 = equalBackward(
                    before,
                    after,
                    before_length - x2 - 1,
                    after_length - y2 - 1,
                );
                if (match) {
                    x2 += d1;
                    y2 += d1;
                } else {
                    break;
                }
            }
            v2.items[@intCast(k2_offset)] = x2;
            if (x2 > before_length) {
                // Ran off the left of the graph.
                k2end += 2;
            } else if (y2 > after_length) {
                // Ran off the top of the graph.
                k2start += 2;
            } else if (!front) {
                const k1_offset = v_offset + delta - k2;
                if (k1_offset >= 0 and k1_offset < v_length and v1.items[@intCast(k1_offset)] != -1) {
                    const x1 = v1.items[@intCast(k1_offset)];
                    const y1 = v_offset + x1 - k1_offset;
                    // Mirror x2 onto top-left coordinate system.
                    x2 = before_length - v2.items[@intCast(k2_offset)];
                    if (x1 >= x2) {
                        // Overlap detected.
                        return dmp.diffBisectSplit(allocator, before, after, x1, y1, deadline);
                    }
                }
            }
        }
    }
    // Diff took too long and hit the deadline or
    // number of diffs equals number of characters, no commonality at all.
    var diffs = DiffList{};
    errdefer deinitDiffList(allocator, &diffs);
    try diffs.ensureUnusedCapacity(allocator, 2);
    diffs.appendAssumeCapacity(Diff.init(
        .delete,
        try allocator.dupe(u8, before),
    ));
    diffs.appendAssumeCapacity(Diff.init(
        .insert,
        try allocator.dupe(u8, after),
    ));
    return diffs;
}

/// Match up to a full character in the forward direction.  Note the
/// goal here: we aren't validating Unicode, we're making sure we don't
/// split code unit sequences.  We might get non-minimal diffs on bad
/// UTF-8, but that's fine.
fn equalForward(
    before: []const u8,
    after: []const u8,
    b_i: isize,
    a_i: isize,
) struct { bool, isize } {
    const b_u: usize = @intCast(b_i);
    const a_u: usize = @intCast(a_i);
    const b1c = before[b_u];
    const a1c = after[a_u];
    if (b1c == a1c) {
        // how many codeunits might we expect?
        // ASCII is easy:
        if (b1c < 0x80) {
            return .{ true, 1 };
        } else {
            switch (b1c) {
                0xc2...0xdf => {
                    // two bytes
                    if (b_u + 1 >= before.len or a_u + 1 >= after.len) {
                        // it's a match ¯\_(ツ)_/¯
                        return .{ true, 1 };
                    } // length is unused for false results
                    return .{ before[b_u + 1] == after[a_u + 1], 2 };
                },
                0xe0...0xef => {
                    // three bytes
                    if (b_u + 2 >= before.len or a_u + 2 >= after.len) {
                        return .{ true, 1 };
                    }
                    const m2 = before[b_u + 1] == after[a_u + 1];
                    const m3 = before[b_u + 2] == after[a_u + 2];
                    return .{ m2 and m3, 3 };
                },
                0xf0...0xf4 => {
                    // four bytes
                    if (b_u + 3 >= before.len or a_u + 3 >= after.len) {
                        return .{ true, 1 };
                    }
                    const m = same: {
                        const m2 = before[b_u + 1] == after[a_u + 1];
                        const m3 = before[b_u + 2] == after[a_u + 2];
                        const m4 = before[b_u + 3] == after[a_u + 3];
                        break :same m2 and m3 and m4;
                    };
                    return .{ m, 4 };
                }, // follow byte or invalid high, doesn't matter, match
                else => return .{ true, 1 },
            }
        }
    } else {
        return .{ false, 0 };
    }
}

/// Match characters backward, avoiding splitting two valid codeunits with a
/// common suffix.  Once again, we are not interested in validating the text,
/// just in preventing a spurious diff which truncates Unicode.
fn equalBackward(
    before: []const u8,
    after: []const u8,
    b_i: isize,
    a_i: isize,
) struct { bool, isize } {
    const b_u: usize = @intCast(b_i);
    const a_u: usize = @intCast(a_i);
    const b1c = before[b_u];
    const a1c = after[a_u];
    if (b1c == a1c) {
        // how many codeunits might we expect?
        // different jam here! We have to match back to a lead:
        switch (b1c) {
            // follow byte might be a code unit sequence
            0x80...0xbf => {
                // I'd rather double the offsets then deal with
                // casting.  Feel free to optimize...
                var off: usize = 1;
                var offi: isize = @intCast(off);
                while (off < 4 and b_i - offi >= 0 and a_i - offi >= 0) {
                    const b = before[b_u - off];
                    if (b != after[b_u - off]) {
                        // whole thing is a fail
                        return .{ false, 0 }; // here the offset doesn't matter
                    }
                    // check for lead byte
                    // since we presume well-formedness, any lead will do
                    if (0xc1 < b and b < 0xf5) {
                        return .{ true, offi + 1 };
                    }
                    off += 1;
                    offi += 1;
                } // since we didn't spot a plausible character, match 1
                return .{ true, 1 };
            }, // ASCII, malformed, don't care,
            else => return .{ true, 1 },
        }
    } else {
        return .{ false, 0 };
    }
}

/// Given the location of the 'middle snake', split the diff in two parts
/// and recurse.
/// @param text1 Old string to be diffed.
/// @param text2 New string to be diffed.
/// @param x Index of split point in text1.
/// @param y Index of split point in text2.
/// @param deadline Time at which to bail if not yet complete.
/// @return LinkedList of Diff objects.
fn diffBisectSplit(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    text1: []const u8,
    text2: []const u8,
    x: isize,
    y: isize,
    deadline: u64,
) DiffError!DiffList {
    const text1a = text1[0..@intCast(x)];
    const text2a = text2[0..@intCast(y)];
    const text1b = text1[@intCast(x)..];
    const text2b = text2[@intCast(y)..];

    // Compute both diffs serially.
    var diffs = try dmp.diffInternal(allocator, text1a, text2a, false, deadline);
    errdefer deinitDiffList(allocator, &diffs);
    var diffs_b = try dmp.diffInternal(allocator, text1b, text2b, false, deadline);
    // Free the list, but not the contents:
    defer diffs_b.deinit(allocator);
    errdefer {
        for (diffs_b.items) |d| {
            allocator.free(d.text);
        }
    }
    try diffs.appendSlice(allocator, diffs_b.items);
    return diffs;
}

/// Do a quick line-level diff on both strings, then rediff the parts for
/// greater accuracy.
/// This speedup can produce non-minimal diffs.
/// @param text1 Old string to be diffed.
/// @param text2 New string to be diffed.
/// @param deadline Time when the diff should be complete by.
/// @return List of Diff objects.
fn diffLineMode(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    text1_in: []const u8,
    text2_in: []const u8,
    deadline: u64,
) DiffError!DiffList {
    // Scan the text on a line-by-line basis first.
    var a = try diffLinesToChars(allocator, text1_in, text2_in);
    defer a.deinit(allocator);
    const text1 = a.chars_1;
    const text2 = a.chars_2;
    const line_array = a.line_array;

    var diffs: DiffList = try dmp.diffInternal(allocator, text1, text2, false, deadline);
    errdefer diffs.deinit(allocator);
    // Convert the diff back to original text.
    try diffCharsToLines(allocator, diffs.items, line_array.items);
    // Eliminate freak matches (e.g. blank lines)
    try diffCleanupSemantic(allocator, &diffs);

    // Rediff any replacement blocks, this time character-by-character.
    // Add a dummy entry at the end.
    try diffs.append(allocator, Diff.init(.equal, ""));

    var pointer: usize = 0;
    var count_delete: usize = 0;
    var count_insert: usize = 0;
    var text_delete = ArrayListUnmanaged(u8){};
    var text_insert = ArrayListUnmanaged(u8){};
    defer {
        text_delete.deinit(allocator);
        text_insert.deinit(allocator);
    }

    while (pointer < diffs.items.len) {
        switch (diffs.items[pointer].operation) {
            .insert => {
                count_insert += 1;
                try text_insert.appendSlice(allocator, diffs.items[pointer].text);
            },
            .delete => {
                count_delete += 1;
                try text_delete.appendSlice(allocator, diffs.items[pointer].text);
            },
            .equal => {
                // Upon reaching an equality, check for prior redundancies.
                if (count_delete >= 1 and count_insert >= 1) {
                    // Delete the offending records and add the merged ones.
                    freeRangeDiffList(
                        allocator,
                        &diffs,
                        pointer - count_delete - count_insert,
                        count_delete + count_insert,
                    );
                    try diffs.replaceRange(
                        allocator,
                        pointer - count_delete - count_insert,
                        count_delete + count_insert,
                        &.{},
                    );
                    pointer = pointer - count_delete - count_insert;
                    var sub_diff = try dmp.diffInternal(allocator, text_delete.items, text_insert.items, false, deadline);
                    defer sub_diff.deinit(allocator);
                    try diffs.insertSlice(allocator, pointer, sub_diff.items);
                    pointer = pointer + sub_diff.items.len;
                }
                count_insert = 0;
                count_delete = 0;
                text_delete.items.len = 0;
                text_insert.items.len = 0;
            },
        }
        pointer += 1;
    }
    diffs.items.len -= 1; // Remove the dummy entry at the end.

    return diffs;
}

const LinesToCharsResult = struct {
    chars_1: []const u8,
    chars_2: []const u8,
    line_array: ArrayListUnmanaged([]const u8),

    pub fn deinit(self: *LinesToCharsResult, allocator: Allocator) void {
        allocator.free(self.chars_1);
        allocator.free(self.chars_2);
        self.line_array.deinit(allocator);
    }
};

/// Split two texts into a list of strings.  Reduce the texts to a string of
/// hashes where each Unicode character represents one line.
/// @param text1 First string.
/// @param text2 Second string.
/// @return Three element Object array, containing the encoded text1, the
///     encoded text2 and the List of unique strings.  The zeroth element
///     of the List of unique strings is intentionally blank.
fn diffLinesToChars(
    allocator: std.mem.Allocator,
    text1: []const u8,
    text2: []const u8,
) DiffError!LinesToCharsResult {
    var line_array = ArrayListUnmanaged([]const u8){};
    errdefer line_array.deinit(allocator);
    var line_hash = std.StringHashMapUnmanaged(usize){};
    defer line_hash.deinit(allocator);
    // e.g. line_array[4] == "Hello\n"
    // e.g. line_hash.get("Hello\n") == 4

    // "\x00" is a valid character, but various debuggers don't like it.
    // So we'll insert a junk entry to avoid generating a null character.
    try line_array.append(allocator, "");

    // Allocate 2/3rds of the space for text1, the rest for text2.
    const chars1 = try diffLinesToCharsMunge(allocator, text1, &line_array, &line_hash, 170);
    const chars2 = try diffLinesToCharsMunge(allocator, text2, &line_array, &line_hash, 255);
    return .{ .chars_1 = chars1, .chars_2 = chars2, .line_array = line_array };
}

/// Split a text into a list of strings.  Reduce the texts to a string of
/// hashes where each Unicode character represents one line.
/// @param text String to encode.
/// @param lineArray List of unique strings.
/// @param lineHash Map of strings to indices.
/// @param maxLines Maximum length of lineArray.
/// @return Encoded string.
fn diffLinesToCharsMunge(
    allocator: std.mem.Allocator,
    text: []const u8,
    line_array: *ArrayListUnmanaged([]const u8),
    line_hash: *std.StringHashMapUnmanaged(usize),
    max_lines: usize,
) DiffError![]const u8 {
    var line_start: isize = 0;
    var line_end: isize = -1;
    var line: []const u8 = undefined;
    var chars = ArrayListUnmanaged(u8){};
    defer chars.deinit(allocator);
    // Walk the text, pulling out a Substring for each line.
    // TODO this can be handled with a Reader, avoiding all the manual splitting
    while (line_end < @as(isize, @intCast(text.len)) - 1) {
        line_end = b: {
            break :b @as(isize, @intCast(std.mem.indexOf(u8, text[@intCast(line_start)..], "\n") orelse
                break :b @intCast(text.len - 1))) + line_start;
        };
        line = text[@intCast(line_start) .. @as(usize, @intCast(line_start)) + @as(usize, @intCast(line_end + 1 - line_start))];

        if (line_hash.get(line)) |value| {
            try chars.append(allocator, @intCast(value));
        } else {
            if (line_array.items.len == max_lines) {
                // Bail out at 255 because char 256 == char 0.
                line = text[@intCast(line_start)..];
                line_end = @intCast(text.len);
            }
            try line_array.append(allocator, line);
            try line_hash.put(allocator, line, line_array.items.len - 1);
            try chars.append(allocator, @intCast(line_array.items.len - 1));
        }
        line_start = line_end + 1;
    }
    return try chars.toOwnedSlice(allocator);
}

/// Rehydrate the text in a diff from a string of line hashes to real lines
/// of text.
/// @param diffs List of Diff objects.
/// @param lineArray List of unique strings.
fn diffCharsToLines(
    allocator: std.mem.Allocator,
    diffs: []Diff,
    line_array: []const []const u8,
) DiffError!void {
    var text = ArrayListUnmanaged(u8){};
    defer text.deinit(allocator);

    for (diffs) |*d| {
        var j: usize = 0;
        while (j < d.text.len) : (j += 1) {
            try text.appendSlice(allocator, line_array[d.text[j]]);
        }
        allocator.free(d.text);
        d.text = try text.toOwnedSlice(allocator);
    }
}

/// Reorder and merge like edit sections.  Merge equalities.
/// Any edit section can move as long as it doesn't cross an equality.
/// @param diffs List of Diff objects.
fn diffCleanupMerge(allocator: std.mem.Allocator, diffs: *DiffList) DiffError!void {
    // Add a dummy entry at the end.
    try diffs.append(allocator, Diff.init(.equal, ""));
    var pointer: usize = 0;
    var count_delete: usize = 0;
    var count_insert: usize = 0;

    var text_delete = ArrayListUnmanaged(u8){};
    defer text_delete.deinit(allocator);

    var text_insert = ArrayListUnmanaged(u8){};
    defer text_insert.deinit(allocator);

    var common_length: usize = undefined;
    while (pointer < diffs.items.len) {
        switch (diffs.items[pointer].operation) {
            .insert => {
                count_insert += 1;
                try text_insert.appendSlice(allocator, diffs.items[pointer].text);
                pointer += 1;
            },
            .delete => {
                count_delete += 1;
                try text_delete.appendSlice(allocator, diffs.items[pointer].text);
                pointer += 1;
            },
            .equal => {
                // Upon reaching an equality, check for prior redundancies.
                if (count_delete + count_insert > 1) {
                    if (count_delete != 0 and count_insert != 0) {
                        // Factor out any common prefixes.
                        common_length = diffCommonPrefix(text_insert.items, text_delete.items);
                        if (common_length != 0) {
                            if ((pointer - count_delete - count_insert) > 0 and
                                diffs.items[pointer - count_delete - count_insert - 1].operation == .equal)
                            { // The prefix is not at the start of the diffs
                                const ii = pointer - count_delete - count_insert - 1;
                                var nt = try allocator.alloc(u8, diffs.items[ii].text.len + common_length);
                                const ot = diffs.items[ii].text;
                                @memcpy(nt[0..ot.len], ot);
                                @memcpy(nt[ot.len..], text_insert.items[0..common_length]);
                                diffs.items[ii].text = nt;
                                allocator.free(ot);
                            } else {
                                try diffs.ensureUnusedCapacity(allocator, 1);
                                const text = try allocator.dupe(u8, text_insert.items[0..common_length]);
                                diffs.insertAssumeCapacity(0, Diff.init(.equal, text));
                                pointer += 1;
                            }
                            try text_insert.replaceRange(allocator, 0, common_length, &.{});
                            try text_delete.replaceRange(allocator, 0, common_length, &.{});
                        }
                        // Factor out any common suffixies.
                        // @ZigPort this seems very wrong
                        common_length = diffCommonSuffix(text_insert.items, text_delete.items);
                        if (common_length != 0) {
                            const old_text = diffs.items[pointer].text;
                            diffs.items[pointer].text = try std.mem.concat(allocator, u8, &.{
                                text_insert.items[text_insert.items.len - common_length ..],
                                old_text,
                            });
                            allocator.free(old_text);
                            text_insert.items.len -= common_length;
                            text_delete.items.len -= common_length;
                        }
                    }
                    // Delete the offending records and add the merged ones.
                    pointer -= count_delete + count_insert;
                    if (count_delete + count_insert > 0) {
                        freeRangeDiffList(allocator, diffs, pointer, count_delete + count_insert);
                        try diffs.replaceRange(allocator, pointer, count_delete + count_insert, &.{});
                    }

                    if (text_delete.items.len != 0) {
                        try diffs.ensureUnusedCapacity(allocator, 1);
                        diffs.insertAssumeCapacity(pointer, Diff.init(
                            .delete,
                            try allocator.dupe(u8, text_delete.items),
                        ));
                        pointer += 1;
                    }
                    if (text_insert.items.len != 0) {
                        try diffs.ensureUnusedCapacity(allocator, 1);
                        diffs.insertAssumeCapacity(pointer, Diff.init(
                            .insert,
                            try allocator.dupe(u8, text_insert.items),
                        ));
                        pointer += 1;
                    }
                    pointer += 1;
                } else if (pointer != 0 and diffs.items[pointer - 1].operation == .equal) {
                    // Merge this equality with the previous one.
                    // Diff texts are []const u8 so a realloc isn't practical here
                    var nt = try allocator.alloc(u8, diffs.items[pointer - 1].text.len + diffs.items[pointer].text.len);
                    const ot = diffs.items[pointer - 1].text;
                    defer (allocator.free(ot));
                    @memcpy(nt[0..ot.len], ot);
                    @memcpy(nt[ot.len..], diffs.items[pointer].text);
                    diffs.items[pointer - 1].text = nt;
                    const dead_diff = diffs.orderedRemove(pointer);
                    allocator.free(dead_diff.text);
                } else {
                    pointer += 1;
                }
                count_insert = 0;
                count_delete = 0;
                text_delete.items.len = 0;
                text_insert.items.len = 0;
            },
        }
    }
    if (diffs.items[diffs.items.len - 1].text.len == 0) {
        diffs.items.len -= 1;
    }

    // Second pass: look for single edits surrounded on both sides by
    // equalities which can be shifted sideways to eliminate an equality.
    // e.g: A<ins>BA</ins>C -> <ins>AB</ins>AC
    var changes = false;
    pointer = 1;
    // Intentionally ignore the first and last element (don't need checking).
    while (pointer < (diffs.items.len - 1)) {
        if (diffs.items[pointer - 1].operation == .equal and
            diffs.items[pointer + 1].operation == .equal)
        {
            // This is a single edit surrounded by equalities.
            if (std.mem.endsWith(u8, diffs.items[pointer].text, diffs.items[pointer - 1].text)) {
                const old_pt = diffs.items[pointer].text;
                const pt = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer - 1].text,
                    diffs.items[pointer].text[0 .. diffs.items[pointer].text.len -
                        diffs.items[pointer - 1].text.len],
                });
                allocator.free(old_pt);
                diffs.items[pointer].text = pt;
                const old_pt1t = diffs.items[pointer + 1].text;
                const p1t = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer - 1].text,
                    diffs.items[pointer + 1].text,
                });
                allocator.free(old_pt1t);
                diffs.items[pointer + 1].text = p1t;
                freeRangeDiffList(allocator, diffs, pointer - 1, 1);
                try diffs.replaceRange(allocator, pointer - 1, 1, &.{});
                changes = true;
            } else if (std.mem.startsWith(u8, diffs.items[pointer].text, diffs.items[pointer + 1].text)) {
                const old_ptm1 = diffs.items[pointer - 1].text;
                const pm1t = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer - 1].text,
                    diffs.items[pointer + 1].text,
                });
                allocator.free(old_ptm1);
                diffs.items[pointer - 1].text = pm1t;
                const old_pt = diffs.items[pointer].text;
                const pt = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer].text[diffs.items[pointer + 1].text.len..],
                    diffs.items[pointer + 1].text,
                });
                allocator.free(old_pt);
                diffs.items[pointer].text = pt;
                freeRangeDiffList(allocator, diffs, pointer + 1, 1);
                try diffs.replaceRange(allocator, pointer + 1, 1, &.{});
                changes = true;
            }
        }
        pointer += 1;
    }
    // If shifts were made, the diff needs reordering and another shift sweep.
    if (changes) {
        try diffCleanupMerge(allocator, diffs);
    }
}

/// Reduce the number of edits by eliminating semantically trivial
/// equalities.
/// @param diffs List of Diff objects.
fn diffCleanupSemantic(allocator: std.mem.Allocator, diffs: *DiffList) DiffError!void {
    var changes = false;
    // Stack of indices where equalities are found.
    var equalities = ArrayListUnmanaged(isize){};
    defer equalities.deinit(allocator);
    // Always equal to equalities[equalitiesLength-1][1]
    var last_equality: ?[]const u8 = null;
    var pointer: isize = 0; // Index of current position.
    // Number of characters that changed prior to the equality.
    var length_insertions1: usize = 0;
    var length_deletions1: usize = 0;
    // Number of characters that changed after the equality.
    var length_insertions2: usize = 0;
    var length_deletions2: usize = 0;
    while (pointer < diffs.items.len) {
        if (diffs.items[@intCast(pointer)].operation == .equal) { // Equality found.
            try equalities.append(allocator, pointer);
            length_insertions1 = length_insertions2;
            length_deletions1 = length_deletions2;
            length_insertions2 = 0;
            length_deletions2 = 0;
            last_equality = diffs.items[@intCast(pointer)].text;
        } else { // an insertion or deletion
            if (diffs.items[@intCast(pointer)].operation == .insert) {
                length_insertions2 += diffs.items[@intCast(pointer)].text.len;
            } else {
                length_deletions2 += diffs.items[@intCast(pointer)].text.len;
            }
            // Eliminate an equality that is smaller or equal to the edits on both
            // sides of it.
            if (last_equality != null and
                (last_equality.?.len <= @max(length_insertions1, length_deletions1)) and
                (last_equality.?.len <= @max(length_insertions2, length_deletions2)))
            {
                // Duplicate record.
                try diffs.ensureUnusedCapacity(allocator, 1);
                diffs.insertAssumeCapacity(
                    @intCast(equalities.items[equalities.items.len - 1]),
                    Diff.init(
                        .delete,
                        try allocator.dupe(u8, last_equality.?),
                    ),
                );
                // Change second copy to insert.
                diffs.items[@intCast(equalities.items[equalities.items.len - 1] + 1)].operation = .insert;
                // Throw away the equality we just deleted.
                _ = equalities.pop();
                if (equalities.items.len > 0) {
                    _ = equalities.pop();
                }
                pointer = if (equalities.items.len > 0) equalities.items[equalities.items.len - 1] else -1;
                length_insertions1 = 0; // Reset the counters.
                length_deletions1 = 0;
                length_insertions2 = 0;
                length_deletions2 = 0;
                last_equality = null;
                changes = true;
            }
        }
        pointer += 1;
    }

    // Normalize the diff.
    if (changes) {
        try diffCleanupMerge(allocator, diffs);
    }
    try diffCleanupSemanticLossless(allocator, diffs);

    // Find any overlaps between deletions and insertions.
    // e.g: <del>abcxxx</del><ins>xxxdef</ins>
    //   -> <del>abc</del>xxx<ins>def</ins>
    // e.g: <del>xxxabc</del><ins>defxxx</ins>
    //   -> <ins>def</ins>xxx<del>abc</del>
    // Only extract an overlap if it is as big as the edit ahead or behind it.
    pointer = 1;
    while (pointer < diffs.items.len) {
        if (diffs.items[@intCast(pointer - 1)].operation == .delete and
            diffs.items[@intCast(pointer)].operation == .insert)
        {
            const deletion = diffs.items[@intCast(pointer - 1)].text;
            const insertion = diffs.items[@intCast(pointer)].text;
            const overlap_length1: usize = diffCommonOverlap(deletion, insertion);
            const overlap_length2: usize = diffCommonOverlap(insertion, deletion);
            if (overlap_length1 >= overlap_length2) {
                if (@as(f32, @floatFromInt(overlap_length1)) >= @as(f32, @floatFromInt(deletion.len)) / 2.0 or
                    @as(f32, @floatFromInt(overlap_length1)) >= @as(f32, @floatFromInt(insertion.len)) / 2.0)
                {
                    // Overlap found.
                    // Insert an equality and trim the surrounding edits.
                    try diffs.ensureUnusedCapacity(allocator, 1);
                    diffs.insertAssumeCapacity(
                        @intCast(pointer),
                        Diff.init(
                            .equal,
                            try allocator.dupe(u8, insertion[0..overlap_length1]),
                        ),
                    );
                    diffs.items[@intCast(pointer - 1)].text =
                        try allocator.dupe(u8, deletion[0 .. deletion.len - overlap_length1]);
                    allocator.free(deletion);
                    diffs.items[@intCast(pointer + 1)].text =
                        try allocator.dupe(u8, insertion[overlap_length1..]);
                    allocator.free(insertion);
                    pointer += 1;
                }
            } else {
                if (@as(f32, @floatFromInt(overlap_length2)) >= @as(f32, @floatFromInt(deletion.len)) / 2.0 or
                    @as(f32, @floatFromInt(overlap_length2)) >= @as(f32, @floatFromInt(insertion.len)) / 2.0)
                {
                    // Reverse overlap found.
                    // Insert an equality and swap and trim the surrounding edits.
                    try diffs.ensureUnusedCapacity(allocator, 1);
                    diffs.insertAssumeCapacity(
                        @intCast(pointer),
                        Diff.init(
                            .equal,
                            try allocator.dupe(u8, deletion[0..overlap_length2]),
                        ),
                    );
                    const new_minus = try allocator.dupe(u8, insertion[0 .. insertion.len - overlap_length2]);
                    errdefer allocator.free(new_minus); // necessary due to swap
                    const new_plus = try allocator.dupe(u8, deletion[overlap_length2..]);
                    allocator.free(deletion);
                    allocator.free(insertion);
                    diffs.items[@intCast(pointer - 1)].operation = .insert;
                    diffs.items[@intCast(pointer - 1)].text = new_minus;
                    diffs.items[@intCast(pointer + 1)].operation = .delete;
                    diffs.items[@intCast(pointer + 1)].text = new_plus;
                    pointer += 1;
                }
            }
            pointer += 1;
        }
        pointer += 1;
    }
}

/// Look for single edits surrounded on both sides by equalities
/// which can be shifted sideways to align the edit to a word boundary.
/// e.g: The c<ins>at c</ins>ame. -> The <ins>cat </ins>came.
pub fn diffCleanupSemanticLossless(
    allocator: std.mem.Allocator,
    diffs: *DiffList,
) DiffError!void {
    var pointer: usize = 1;
    // Intentionally ignore the first and last element (don't need checking).
    while (pointer < @as(isize, @intCast(diffs.items.len)) - 1) {
        if (diffs.items[pointer - 1].operation == .equal and
            diffs.items[pointer + 1].operation == .equal)
        {
            // This is a single edit surrounded by equalities.
            var equality_1 = std.ArrayListUnmanaged(u8){};
            defer equality_1.deinit(allocator);
            try equality_1.appendSlice(allocator, diffs.items[pointer - 1].text);

            var edit = std.ArrayListUnmanaged(u8){};
            defer edit.deinit(allocator);
            try edit.appendSlice(allocator, diffs.items[pointer].text);

            var equality_2 = std.ArrayListUnmanaged(u8){};
            defer equality_2.deinit(allocator);
            try equality_2.appendSlice(allocator, diffs.items[pointer + 1].text);

            // First, shift the edit as far left as possible.
            const common_offset = diffCommonSuffix(equality_1.items, edit.items);
            if (common_offset > 0) {
                // TODO: Use buffer
                const common_string = try allocator.dupe(u8, edit.items[edit.items.len - common_offset ..]);
                defer allocator.free(common_string);

                equality_1.items.len = equality_1.items.len - common_offset;

                // edit.items.len = edit.items.len - common_offset;
                const not_common = try allocator.dupe(u8, edit.items[0 .. edit.items.len - common_offset]);
                defer allocator.free(not_common);

                edit.clearRetainingCapacity();
                try edit.appendSlice(allocator, common_string);
                try edit.appendSlice(allocator, not_common);

                try equality_2.insertSlice(allocator, 0, common_string);
            }

            // Second, step character by character right,
            // looking for the best fit.
            var best_equality_1 = ArrayListUnmanaged(u8){};
            defer best_equality_1.deinit(allocator);
            try best_equality_1.appendSlice(allocator, equality_1.items);

            var best_edit = ArrayListUnmanaged(u8){};
            defer best_edit.deinit(allocator);
            try best_edit.appendSlice(allocator, edit.items);

            var best_equality_2 = ArrayListUnmanaged(u8){};
            defer best_equality_2.deinit(allocator);
            try best_equality_2.appendSlice(allocator, equality_2.items);

            var best_score = diffCleanupSemanticScore(equality_1.items, edit.items) +
                diffCleanupSemanticScore(edit.items, equality_2.items);

            while (edit.items.len != 0 and equality_2.items.len != 0 and edit.items[0] == equality_2.items[0]) {
                try equality_1.append(allocator, edit.items[0]);

                _ = edit.orderedRemove(0);
                try edit.append(allocator, equality_2.items[0]);

                _ = equality_2.orderedRemove(0);

                const score = diffCleanupSemanticScore(equality_1.items, edit.items) +
                    diffCleanupSemanticScore(edit.items, equality_2.items);
                // The >= encourages trailing rather than leading whitespace on
                // edits.
                if (score >= best_score) {
                    best_score = score;

                    best_equality_1.items.len = 0;
                    try best_equality_1.appendSlice(allocator, equality_1.items);

                    best_edit.items.len = 0;
                    try best_edit.appendSlice(allocator, edit.items);

                    best_equality_2.items.len = 0;
                    try best_equality_2.appendSlice(allocator, equality_2.items);
                }
            }

            if (!std.mem.eql(u8, diffs.items[pointer - 1].text, best_equality_1.items)) {
                // We have an improvement, save it back to the diff.
                if (best_equality_1.items.len != 0) {
                    const old_text = diffs.items[pointer - 1].text;
                    diffs.items[pointer - 1].text = try allocator.dupe(u8, best_equality_1.items);
                    allocator.free(old_text);
                } else {
                    const old_diff = diffs.orderedRemove(pointer - 1);
                    allocator.free(old_diff.text);
                    pointer -= 1;
                }
                const old_text1 = diffs.items[pointer].text;
                diffs.items[pointer].text = try allocator.dupe(u8, best_edit.items);
                defer allocator.free(old_text1);
                if (best_equality_2.items.len != 0) {
                    const old_text2 = diffs.items[pointer + 1].text;
                    diffs.items[pointer + 1].text = try allocator.dupe(u8, best_equality_2.items);
                    allocator.free(old_text2);
                } else {
                    const old_diff = diffs.orderedRemove(pointer + 1);
                    allocator.free(old_diff.text);
                    pointer -= 1;
                }
            }
        }
        pointer += 1;
    }
}

/// Given two strings, compute a score representing whether the internal
/// boundary falls on logical boundaries.
/// Scores range from 6 (best) to 0 (worst).
/// @param one First string.
/// @param two Second string.
/// @return The score.
fn diffCleanupSemanticScore(one: []const u8, two: []const u8) usize {
    if (one.len == 0 or two.len == 0) {
        // Edges are the best.
        return 6;
    }

    // Each port of this function behaves slightly differently due to
    // subtle differences in each language's definition of things like
    // 'whitespace'.  Since this function's purpose is largely cosmetic,
    // the choice has been made to use each language's native features
    // rather than force total conformity.
    const char1 = one[one.len - 1];
    const char2 = two[0];
    const nonAlphaNumeric1 = !std.ascii.isAlphanumeric(char1);
    const nonAlphaNumeric2 = !std.ascii.isAlphanumeric(char2);
    const whitespace1 = nonAlphaNumeric1 and std.ascii.isWhitespace(char1);
    const whitespace2 = nonAlphaNumeric2 and std.ascii.isWhitespace(char2);
    const lineBreak1 = whitespace1 and std.ascii.isControl(char1);
    const lineBreak2 = whitespace2 and std.ascii.isControl(char2);
    const blankLine1 = lineBreak1 and
        // BLANKLINEEND.IsMatch(one);
        (std.mem.endsWith(u8, one, "\n\n") or std.mem.endsWith(u8, one, "\n\r\n"));
    const blankLine2 = lineBreak2 and
        // BLANKLINESTART.IsMatch(two);
        (std.mem.startsWith(u8, two, "\n\n") or
        std.mem.startsWith(u8, two, "\r\n\n") or
        std.mem.startsWith(u8, two, "\n\r\n") or
        std.mem.startsWith(u8, two, "\r\n\r\n"));

    if (blankLine1 or blankLine2) {
        // Five points for blank lines.
        return 5;
    } else if (lineBreak1 or lineBreak2) {
        // Four points for line breaks.
        return 4;
    } else if (nonAlphaNumeric1 and !whitespace1 and whitespace2) {
        // Three points for end of sentences.
        return 3;
    } else if (whitespace1 or whitespace2) {
        // Two points for whitespace.
        return 2;
    } else if (nonAlphaNumeric1 or nonAlphaNumeric2) {
        // One point for non-alphanumeric.
        return 1;
    }
    return 0;
}

inline fn boolInt(b: bool) usize {
    if (b)
        return 1
    else
        return 0;
}

/// Reduce the number of edits by eliminating operationally trivial
/// equalities.  TODO this needs tests
pub fn diffCleanupEfficiency(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    diffs: *DiffList,
) DiffError!void {
    var changes = false;
    // Stack of indices where equalities are found.
    var equalities = std.ArrayList(usize).init(allocator);
    defer equalities.deinit();
    // Always equal to equalities[equalitiesLength-1][1]
    var last_equality: []const u8 = "";
    var ipointer: isize = 0; // Index of current position.
    // Is there an insertion operation before the last equality.
    var pre_ins = false;
    // Is there a deletion operation before the last equality.
    var pre_del = false;
    // Is there an insertion operation after the last equality.
    var post_ins = false;
    // Is there a deletion operation after the last equality.
    var post_del = false;
    while (ipointer < diffs.items.len) {
        const pointer: usize = @intCast(ipointer);
        if (diffs.items[pointer].operation == .equal) { // Equality found.
            if (diffs.items[pointer].text.len < dmp.diff_edit_cost and (post_ins or post_del)) {
                // Candidate found.
                try equalities.append(pointer);
                pre_ins = post_ins;
                pre_del = post_del;
                last_equality = diffs.items[pointer].text;
            } else {
                // Not a candidate, and can never become one.
                equalities.items.len = 0;
                last_equality = "";
            }
            post_ins = false;
            post_del = false;
        } else { // An insertion or deletion.
            if (diffs.items[pointer].operation == .delete) {
                post_del = true;
            } else {
                post_ins = true;
            }
            // Five types to be split:
            // <ins>A</ins><del>B</del>XY<ins>C</ins><del>D</del>
            // <ins>A</ins>X<ins>C</ins><del>D</del>
            // <ins>A</ins><del>B</del>X<ins>C</ins>
            // <ins>A</del>X<ins>C</ins><del>D</del>
            // <ins>A</ins><del>B</del>X<del>C</del>
            if ((last_equality.len != 0) and
                ((pre_ins and pre_del and post_ins and post_del) or
                ((last_equality.len < dmp.diff_edit_cost / 2) and
                (boolInt(pre_ins) + boolInt(pre_del) + boolInt(post_ins) + boolInt(post_del) == 3))))
            {
                // Duplicate record.
                try diffs.ensureUnusedCapacity(allocator, 1);
                diffs.insertAssumeCapacity(
                    equalities.items[equalities.items.len - 1],
                    Diff.init(
                        .delete,
                        try allocator.dupe(u8, last_equality),
                    ),
                );
                // Change second copy to insert.
                diffs.items[equalities.items[equalities.items.len - 1] + 1].operation = .insert;
                _ = equalities.pop(); // Throw away the equality we just deleted.
                last_equality = "";
                if (pre_ins and pre_del) {
                    // No changes made which could affect previous entry, keep going.
                    post_ins = true;
                    post_del = true;
                    equalities.items.len = 0;
                } else {
                    if (equalities.items.len > 0) {
                        _ = equalities.pop();
                    }

                    ipointer = if (equalities.items.len > 0) @intCast(equalities.items[equalities.items.len - 1]) else -1;
                    post_ins = false;
                    post_del = false;
                }
                changes = true;
            }
        }
        ipointer += 1;
    }

    if (changes) {
        try diffCleanupMerge(allocator, diffs);
    }
}

/// Determine if the suffix of one string is the prefix of another.
/// @param text1 First string.
/// @param text2 Second string.
/// @return The number of characters common to the end of the first
///     string and the start of the second string.
fn diffCommonOverlap(text1_in: []const u8, text2_in: []const u8) usize {
    var text1 = text1_in;
    var text2 = text2_in;

    // Cache the text lengths to prevent multiple calls.
    const text1_length = text1.len;
    const text2_length = text2.len;
    // Eliminate the null case.
    if (text1_length == 0 or text2_length == 0) {
        return 0;
    }
    // Truncate the longer string.
    if (text1_length > text2_length) {
        text1 = text1[text1_length - text2_length ..];
    } else if (text1_length < text2_length) {
        text2 = text2[0..text1_length];
    }
    const text_length = @min(text1_length, text2_length);
    // Quick check for the worst case.
    if (std.mem.eql(u8, text1, text2)) {
        return text_length;
    }

    // Start by looking for a single character match
    // and increase length until no match is found.
    // Performance analysis: https://neil.fraser.name/news/2010/11/04/
    var best: usize = 0;
    var length: usize = 1;
    const best_idx = idx: while (true) {
        const pattern = text1[text_length - length ..];
        const found = std.mem.indexOf(u8, text2, pattern) orelse
            break :idx best;

        length += found;

        if (found == 0 or std.mem.eql(u8, text1[text_length - length ..], text2[0..length])) {
            best = length;
            length += 1;
        }
    };
    if (best_idx == 0) return best_idx;
    // This would mean a truncation: lead or follow, followed by a follow
    // which differs (or it would be included in our overlap)
    if (text2[best_idx] >= 0x80 and is_follow(text2[best_idx + 1])) {
        // back out
        assert(best_idx == best);
        if (!is_follow(text2[best])) {
            // It's a lead, one back is fine
            return best - 1;
        }
        best -= 1;
        if (best == 0) return 0;
        // It's ok to get no overlap, so we ignore malformation:
        // a bunch of follows could walk back to zero, and that's
        // fine with us
        while (is_follow(text2[best])) {
            best -= 1;
            if (best == 0) return 0;
        }
        // should be a lead, but ASCII is fine, so
        if (text2[best] < 0x80) {
            return best;
        } else {
            return best - 1;
        }
    }
    return best_idx;
}

/// loc is a location in text1, compute and return the equivalent location in
/// text2.
/// e.g. "The cat" vs "The big cat", 1->1, 5->8
/// @param diffs List of Diff objects.
/// @param loc Location within text1.
/// @return Location within text2.
///
pub fn diffIndex(diffs: DiffList, u_loc: usize) usize {
    var chars1: isize = 0;
    var chars2: isize = 0;
    var last_chars1: isize = 0;
    var last_chars2: isize = 0;
    const loc: isize = @intCast(u_loc);
    //  Dummy diff
    var last_diff: Diff = Diff{ .operation = .equal, .text = "" };
    for (diffs.items) |a_diff| {
        if (a_diff.operation != .insert) {
            // Equality or deletion.
            chars1 += @intCast(a_diff.text.len);
        }
        if (a_diff.operation != .delete) {
            // Equality or insertion.
            chars2 += @intCast(a_diff.text.len);
        }
        if (chars1 > loc) {
            // Overshot the location.
            last_diff = a_diff;
            break;
        }
    }
    last_chars1 = chars1;
    last_chars2 = chars2;

    if (last_diff.text.len != 0 and last_diff.operation == .delete) {
        // The location was deleted.
        return @intCast(last_chars2);
    }
    // Add the remaining character length.
    return @intCast(last_chars2 + (loc - last_chars1));
}

/// A struct holding bookends for `diffPrittyFormat(diffs)`.
///
/// May include a function taking an allocator and the Diff,
/// which shall return the text of the Diff, appropriately munged.
/// This allows for tasks like proper HTML escaping.  Note that if
/// the function is provided, all text returned will be freed, so
/// it should always return a copy whether or not edits are needed.
pub const DiffDecorations = struct {
    delete_start: []const u8 = "",
    delete_end: []const u8 = "",
    insert_start: []const u8 = "",
    insert_end: []const u8 = "",
    equals_start: []const u8 = "",
    equals_end: []const u8 = "",
    pre_process: ?fn (Allocator, Diff) error{OutOfMemory}![]const u8 = null,
};

/// Decorations for classic Xterm printing: red for delete and
/// green for insert.
pub const xterm_classic = DiffDecorations{
    .delete_start = "\x1b[91m",
    .delete_end = "\x1b[m",
    .insert_start = "\x1b[92m",
    .insert_end = "\x1b[m",
};

/// Return text representing a pretty-formatted `DiffList`.
/// See `DiffDecorations` for how to customize this output.
pub fn diffPrettyFormat(
    allocator: Allocator,
    diffs: DiffList,
    deco: DiffDecorations,
) ![]const u8 {
    var out = std.ArrayList(u8).init(allocator);
    defer out.deinit();
    const writer = out.writer();
    _ = try writeDiffPrettyFormat(allocator, writer, diffs, deco);
    return out.toOwnedSlice();
}

/// Pretty-print a diff for output to a terminal.
pub fn diffPrettyFormatXTerm(allocator: Allocator, diffs: DiffList) ![]const u8 {
    return try diffPrettyFormat(allocator, diffs, xterm_classic);
}

/// Write a pretty-formatted `DiffList` to `writer`.  The `Allocator`
/// is only used if a custom text formatter is defined for
/// `DiffDecorations`.  Returns number of bytes written.
pub fn writeDiffPrettyFormat(
    allocator: Allocator,
    writer: anytype,
    diffs: DiffList,
    deco: DiffDecorations,
) !usize {
    var written: usize = 0;
    for (diffs.items) |d| {
        const text = if (deco.pre_process) |lambda|
            try lambda(allocator, d)
        else
            d.text;
        defer {
            if (deco.pre_process) |_|
                allocator.free(text);
        }
        switch (d.operation) {
            .delete => {
                //
                written += try writer.write(deco.delete_start);
                written += try writer.write(text);
                written += try writer.write(deco.delete_end);
            },
            .insert => {
                written += try writer.write(deco.insert_start);
                written += try writer.write(text);
                written += try writer.write(deco.insert_end);
            },
            .equal => {
                written += try writer.write(deco.equals_start);
                written += try writer.write(text);
                written += try writer.write(deco.equals_end);
            },
        }
    }
    return written;
}

///
/// Compute and return the source text (all equalities and deletions).
/// @param diffs List of `Diff` objects.
/// @return Source text.
///
pub fn diffBeforeText(allocator: Allocator, diffs: DiffList) ![]const u8 {
    var chars = ArrayListUnmanaged(u8){};
    defer chars.deinit(allocator);
    for (diffs.items) |d| {
        if (d.operation != .insert) {
            try chars.appendSlice(allocator, d.text);
        }
    }
    return chars.toOwnedSlice(allocator);
}

///
/// Compute and return the destination text (all equalities and insertions).
/// @param diffs List of `Diff` objects.
/// @return Destination text.
///
pub fn diffAfterText(allocator: Allocator, diffs: DiffList) ![]const u8 {
    var chars = ArrayListUnmanaged(u8){};
    defer chars.deinit(allocator);
    for (diffs.items) |d| {
        if (d.operation != .delete) {
            try chars.appendSlice(allocator, d.text);
        }
    }
    return chars.toOwnedSlice(allocator);
}

///
/// Compute the Levenshtein distance; the number of inserted,
/// deleted or substituted characters.
///
/// @param diffs List of Diff objects.
/// @return Number of changes.
///
pub fn diffLevenshtein(diffs: DiffList) usize {
    var inserts: usize = 0;
    var deletes: usize = 0;
    var levenshtein: usize = 0;
    for (diffs.items) |a_diff| {
        switch (a_diff.operation) {
            .insert => {
                inserts += a_diff.text.len;
            },
            .delete => {
                deletes += a_diff.text.len;
            },
            .equal => {
                // A deletion and an insertion is one substitution.
                levenshtein = @max(inserts, deletes);
                inserts = 0;
                deletes = 0;
            },
        }
    }

    return levenshtein + @max(inserts, deletes);
}

test diffLevenshtein {
    const allocator = testing.allocator;
    // These diffs don't get text freed
    {
        var diffs = DiffList{};
        defer diffs.deinit(allocator);
        try diffs.appendSlice(allocator, &.{
            Diff.init(.delete, "abc"),
            Diff.init(.insert, "1234"),
            Diff.init(.equal, "xyz"),
        });
        try testing.expectEqual(4, diffLevenshtein(diffs));
    }
    {
        var diffs = DiffList{};
        defer diffs.deinit(allocator);
        try diffs.appendSlice(allocator, &.{
            Diff.init(.equal, "xyz"),
            Diff.init(.delete, "abc"),
            Diff.init(.insert, "1234"),
        });
        try testing.expectEqual(4, diffLevenshtein(diffs));
    }
    {
        var diffs = DiffList{};
        defer diffs.deinit(allocator);
        try diffs.appendSlice(allocator, &.{
            Diff.init(.delete, "abc"),
            Diff.init(.equal, "xyz"),
            Diff.init(.insert, "1234"),
        });
        try testing.expectEqual(7, diffLevenshtein(diffs));
    }
}

//| MATCH FUNCTIONS

/// Locate the best instance of 'pattern' in 'text' near 'loc'.
/// Returns -1 if no match found.
/// @param text The text to search.
/// @param pattern The pattern to search for.
/// @param loc The location to search around.
/// @return Best match index or -1.
pub fn matchMain(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    text: []const u8,
    pattern: []const u8,
    passed_loc: usize,
) !?usize {
    // Clamp the loc to fit within text.
    const loc = @min(passed_loc, text.len);
    if (std.mem.eql(u8, text, pattern)) {
        // Shortcut (potentially not guaranteed by the algorithm)
        // TODO would be good to know what the above means...
        return 0;
    } else if (text.len == 0) {
        // Nothing to match.
        return null;
    } else if (loc + pattern.len <= text.len and std.mem.eql(u8, text[loc .. loc + pattern.len], pattern)) {
        // Perfect match at the perfect spot!  (Includes case of null pattern)
        return loc;
    } else {
        // Do a fuzzy compare.
        return dmp.matchBitap(allocator, text, pattern, loc);
    }
}

// TODO doubling the bits to fit in usize is nice and all, but there's no
// reason to be limited to that, we have bitsets which can be as large as
// we'd like.  This could be passed a comptime power-of-two size, and use
// that to make an ArrayBitSet specialized for several sizes, up to, IDK,
// 2k?  Then split very large patches only.  64, 256, 512, 1024, 2028, is
// a nice balance between code size and versatility.
// Something like this:
fn matchBitapImproved(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    text: []const u8,
    pattern: []const u8,
    loc: usize,
    UIntType: type,
) ?usize {
    assert(pattern.len < @bitSizeOf(UIntType));
    const ShiftWidth = ShiftSizeForType(UIntType);
    // Initialise the alphabet.
    var map = try matchAlphabet(allocator, pattern);
    defer map.deinit();
    // Highest score beyond which we give up.
    var threshold = dmp.threshold;
    // Is there a nearby exact match? (speedup)
    var best_loc = std.mem.indexOfPos(u8, text, pattern);
    if (best_loc) |best| {
        threshold = @min(dmp.matchBitapScore(0, best, loc, pattern), threshold);
    }
    // What about in the other direction? (speedup)
    const trunc_text = text[0..@min(loc + pattern.len, text.len)];
    best_loc = std.mem.lastIndexOf(u8, trunc_text, pattern);
    if (best_loc) |best| {
        threshold = @min(dmp.matchBitapScore(0, best, loc, pattern), threshold);
    }
    // Initialise the bit arrays.
    const shift: ShiftWidth = @intCast(pattern.len - 1);
    // 0 for a match for faster bit twiddles
    const matchmask = ~(1 << shift);
    best_loc = null;
    var bin_min: usize = undefined;
    var bin_mid: usize = undefined;
    var bin_max = pattern.len + text.len;
    // null last_rd to simplying freeing memory
    var last_rd = try allocator.alloc(UIntType, 0);
    for (0..pattern.len) |d| {
        // Scan for the best match; each iteration allows for one more error.
        // Run a binary search to determine how far from 'loc' we can stray at
        // this error level.
        bin_min = 0;
        bin_mid = bin_max;
        while (bin_min < bin_mid) {
            if (dmp.matchBitapScore(d, loc + bin_mid, loc, pattern) <= threshold) {
                bin_min = bin_mid;
            } else {
                bin_max = bin_mid;
            }
            bin_mid = (bin_max - bin_min) / 2 + bin_min;
        }
        // Use the result from this iteration as the maximum for the next.
        bin_max = bin_mid;
        var start = @max(1, loc - bin_mid + 1);
        const finish = @min(loc + bin_mid, text.len) + pattern.len;
        var rd = try allocator.alloc(UIntType, finish + 2);
        const dshift: ShiftWidth = @intCast(d);
        rd[finish + 1] = 1 << dshift;
        var j = finish;
        while (j >= start) : (j -= 1) {
            const char_match: usize = if (text.len <= j - 1 or !map.contains(text[j - 1]))
                // Out of range.
                0
            else
                map.get(text[j - 1]);
            if (d == 0) {
                // First pass: exact match.
                rd[j] = ((rd[j + 1] << 1)) & char_match;
            } else {
                // Subsequent passes: fuzzy match.
                rd[j] = ((rd[j + 1] << 1)) & char_match & (((last_rd[j + 1] & last_rd[j]) << 1)) & last_rd[j + 1];
            }
            if ((rd[j] & matchmask) != 0) {
                const score = dmp.matchBitapScore(d, j - 1, loc, pattern);
                // This match will almost certainly be better than any existing
                // match.  But check anyway.
                if (score <= threshold) {
                    // Told you so.
                    threshold = score;
                    best_loc = j - 1;
                    if (best_loc > loc) {
                        // When passing loc, don't exceed our current distance from loc.
                        start = @max(1, 2 * loc - best_loc);
                    } else {
                        // Already passed loc, downhill from here on in.
                        break;
                    }
                }
            }
        }
        if (dmp.matchBitapScore(d + 1, loc, loc, pattern) > threshold) {
            // No hope for a (better) match at greater error levels.
            break;
        }
        allocator.free(last_rd);
        last_rd = rd;
    }
    allocator.free(last_rd);
    return best_loc;
}

fn ShiftSizeForType(T: type) type {
    return switch (@typeInfo(T.Int.bits)) {
        64 => u6,
        256 => u8,
        1024 => u9,
        2048 => u10,
        else => unreachable,
    };
}

const sh_one: u64 = 1;

/// Locate the best instance of `pattern` in `text` near `loc` using the
/// Bitap algorithm.  Returns -1 if no match found.
///
/// @param text The text to search.
/// @param pattern The pattern to search for.
/// @param loc The location to search around.
/// @return Best match index or -1.
fn matchBitap(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    text: []const u8,
    pattern: []const u8,
    loc: usize,
) !?usize {
    // TODO decide what to do here:
    // assert (Match_MaxBits == 0 || pattern.Length <= Match_MaxBits)
    //    : "Pattern too long for this application.";
    assert(text.len != 0 and pattern.len != 0);

    // Initialise the alphabet.
    var map = try matchAlphabet(allocator, pattern);
    defer map.deinit();
    // Highest score beyond which we give up.
    var score_threshold = dmp.match_threshold;
    // Is there a nearby exact match? (speedup)
    // TODO obviously if we want a speedup here, we do this:
    // if (threshold == 0.0) return best_loc;  #proof in comments
    // We don't have to unwrap best_loc because the retval is ?usize already
    // #proof axiom: threshold is between 0.0 and 1.0 (doc comment)
    var best_loc = std.mem.indexOfPos(u8, text, loc, pattern);
    if (best_loc) |best| { // #proof this returns 0.0 for exact match (see comments in function)
        score_threshold = @min(dmp.matchBitapScore(0, best, loc, pattern), score_threshold);
    }
    // What about in the other direction? (speedup)
    const trunc_text = text[0..@min(loc + pattern.len, text.len)];
    best_loc = std.mem.lastIndexOf(u8, trunc_text, pattern);
    if (best_loc) |best| { // #proof same here obviously
        score_threshold = @min(dmp.matchBitapScore(0, best, loc, pattern), score_threshold);
    }
    // Initialise the bit arrays.
    const shift: u6 = @intCast(pattern.len - 1);
    const matchmask = sh_one << shift;
    best_loc = null;
    // Zig is very insistent about integer width and signedness.
    const i_textlen: isize = @intCast(text.len);
    const i_patlen: isize = @intCast(pattern.len);

    const i_loc: isize = @intCast(loc);
    var bin_min: isize = undefined;
    var bin_mid: isize = undefined;
    var bin_max: isize = i_patlen + i_textlen;
    // null last_rd to simplify freeing memory
    var last_rd: []usize = try allocator.alloc(usize, 0);
    errdefer allocator.free(last_rd);
    for (0..pattern.len) |d| {
        // Scan for the best match; each iteration allows for one more error.
        // Run a binary search to determine how far from 'loc' we can stray at
        // this error level.
        bin_min = 0;
        bin_mid = bin_max;
        while (bin_min < bin_mid) {
            // #proof lemma: if threshold == 0.0, this never happens
            if (dmp.matchBitapScore(d, @intCast(i_loc + bin_mid), loc, pattern) <= score_threshold) {
                bin_min = bin_mid;
            } else {
                bin_max = bin_mid;
            }
            bin_mid = @divTrunc(bin_max - bin_min, 2) + bin_min;
        }
        // Use the result from this iteration as the maximum for the next.
        bin_max = bin_mid;
        var start: usize = @intCast(@max(1, i_loc - bin_mid + 1));
        const finish: usize = @intCast(@min(i_loc + bin_mid, i_textlen) + i_patlen);
        // No errors below this point, so no errdefer either:
        var rd: []usize = try allocator.alloc(usize, finish + 2);
        errdefer allocator.free(rd);
        const dshift: u6 = @intCast(d);
        rd[finish + 1] = (sh_one << dshift) - 1;
        var j = finish;
        while (j >= start) : (j -= 1) {
            const char_match: usize = if (text.len <= j - 1 or !map.contains(text[j - 1]))
                // Out of range.
                0
            else
                map.get(text[j - 1]).?;
            if (d == 0) {
                // First pass: exact match.
                rd[j] = ((rd[j + 1] << 1) | 1) & char_match;
            } else {
                // Subsequent passes: fuzzy match.
                rd[j] = ((rd[j + 1] << 1) | 1) & char_match | (((last_rd[j + 1] | last_rd[j]) << 1) | 1) | last_rd[j + 1];
            }
            if ((rd[j] & matchmask) != 0) {
                const score = dmp.matchBitapScore(d, j - 1, loc, pattern);
                // This match will almost certainly be better than any existing
                // match.  But check anyway.
                // #proof: the smoking gun. This can only be equal not less.
                if (score <= score_threshold) {
                    // Told you so.
                    score_threshold = score;
                    best_loc = j - 1;
                    if (best_loc.? > loc) {
                        // When passing loc, don't exceed our current distance from loc.
                        const i_best_loc: isize = @intCast(best_loc.?);
                        start = @max(1, 2 * i_loc - i_best_loc);
                    } else {
                        // Already passed loc, downhill from here on in.
                        break;
                    }
                }
            }
        } // #proof Anything else will do this.
        // #proof d + 1 starts at 1, so (see function) this will always break.
        if (dmp.matchBitapScore(d + 1, loc, loc, pattern) > score_threshold) {
            // No hope for a (better) match at greater error levels.
            allocator.free(rd);
            break;
        }
        allocator.free(last_rd);
        last_rd = rd;
    }
    allocator.free(last_rd);
    return best_loc;
}

/// Compute and return the score for a match with e errors and x location.
/// @param e Number of errors in match.
/// @param x Location of match.
/// @param loc Expected location of match.
/// @param pattern Pattern being sought.
/// @return Overall score for match (0.0 = good, 1.0 = bad).
fn matchBitapScore(
    dmp: DiffMatchPatch,
    e: usize,
    x: usize,
    loc: usize,
    pattern: []const u8,
) f64 {
    // shortcut? TODO, proof in comments
    // if (e == 0 and x == loc) return 0.0;
    const e_float: f64 = @floatFromInt(e);
    const len_float: f64 = @floatFromInt(pattern.len);
    // if e == 0, accuracy == 0: 0/x = 0
    const accuracy = e_float / len_float;
    // if loc == x, proximity == 0
    const proximity = if (loc >= x) loc - x else x - loc;
    if (dmp.match_distance == 0) {
        // Dodge divide by zero
        if (proximity == 0) // therefore this returns 0
            return accuracy
        else
            return 1.0;
    }
    const float_match: f64 = @floatFromInt(dmp.match_distance);
    const float_proximity: f64 = @floatFromInt(proximity);
    // or this is 0 + 0/f_m aka 0
    return accuracy + (float_proximity / float_match);
}

/// Initialise the alphabet for the Bitap algorithm.
/// @param pattern The text to encode.
/// @return Hash of character locations.
fn matchAlphabet(allocator: Allocator, pattern: []const u8) !std.AutoHashMap(u8, usize) {
    var map = std.AutoHashMap(u8, usize).init(allocator);
    errdefer map.deinit();
    for (pattern) |c| {
        if (!map.contains(c)) {
            try map.put(c, 0);
        }
    }
    for (pattern, 0..) |c, i| {
        const shift: u6 = @intCast(pattern.len - i - 1);
        const value: usize = map.get(c).? | (@as(usize, 1) << shift);
        try map.put(c, value);
    }
    return map;
}

/// Initialise the alphabet for the Bitap algorithm.
/// @param pattern The text to encode.
/// @return Hash of character locations.
fn matchAlphabetImproved(allocator: Allocator, pattern: []const u8, UIntSize: type) !std.HashMap(u8, usize) {
    const ShiftType = ShiftSizeForType(UIntSize);
    var map = std.HashMap(u8, usize).init(allocator);
    errdefer map.deinit();
    for (pattern) |c| {
        if (!map.contains(c)) {
            try map.put(c, 0);
        }
    }
    for (pattern, 0..) |c, i| {
        const shift: ShiftType = @intCast(pattern.len - i - 1);
        // TODO I think we want c_mask & ~ 1 << shift here:
        const value: UIntSize = map.get(c) | (1 << shift);
        try map.put(c, value);
    }
    return map;
}

//|  PATCH FUNCTIONS

/// Increase the context until it is unique, but don't let the pattern
/// expand beyond DiffMatchPatch.match_max_bits.
///
/// @param patch The patch to grow.
/// @param text Source text.
fn patchAddContext(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    patch: *Patch,
    text: []const u8,
) !void {
    if (text.len == 0) return;
    // TODO the fixup logic here might make patterns too large?
    // It should be ok, because big patches get broken up.  Hmm.
    // Also, the SimpleNote maintained branch does it this way.
    var padding: usize = 0;
    { // Grow the pattern around the patch until unique, to set padding amount.
        var pattern = text[patch.start2 .. patch.start2 + patch.length1];
        const max_width: usize = dmp.match_max_bits - (2 * dmp.patch_margin);
        while (std.mem.indexOf(u8, text, pattern) != std.mem.lastIndexOf(u8, text, pattern) and pattern.len < max_width) {
            padding += dmp.patch_margin;
            const pat_start = if (padding > patch.start2) 0 else patch.start2 - padding;
            const pat_end = pat_start + @min(text.len, patch.start2 + patch.length1 + padding);
            pattern = text[pat_start..pat_end];
        }
    }
    // Add one chunk for good luck.
    padding += dmp.patch_margin;
    // Add the prefix.
    const prefix = pre: {
        var pre_start = if (padding > patch.start2) 0 else patch.start2 - padding;
        // Make sure we're not breaking a codepoint.
        while (is_follow(text[pre_start]) and pre_start > 0) {
            pre_start -= 1;
        } // Assuming we did everything else right, pre_end should be
        // properly placed.
        break :pre text[pre_start..patch.start2];
    };
    if (prefix.len != 0) {
        try patch.diffs.ensureUnusedCapacity(allocator, 1);
        patch.diffs.insertAssumeCapacity(0, Diff.init(
            .equal,
            try allocator.dupe(u8, prefix),
        ));
    }
    // Add the suffix.
    const suffix = post: {
        const post_start = patch.start2 + patch.length1;
        // In case we messed up somewhere:
        assert(!is_follow(text[post_start]));
        var post_end = @min(text.len, patch.start2 + patch.length1 + padding);
        // Prevent broken codepoints here as well: Lead bytes, or follow with another follow
        while (post_end + 1 < text.len and !std.ascii.isASCII(text[post_end]) and is_follow(text[post_end + 1])) {
            post_end += 1;
            // Special case: penultimate with another follow at end
            if (post_end + 2 == text.len and is_follow(text[post_end + 1])) {
                post_end += 1;
                break; // Not actually necessary, but polite.
            }
        }
        post_end = @min(post_end, text.len);
        break :post text[post_start..post_end];
    };
    if (suffix.len != 0) {
        try patch.diffs.ensureUnusedCapacity(allocator, 1);
        patch.diffs.appendAssumeCapacity(
            Diff.init(
                .equal,
                try allocator.dupe(u8, suffix),
            ),
        );
    }
    // Roll back the start points.
    patch.start1 -= prefix.len;
    patch.start2 -= prefix.len;
    // Extend the lengths.
    patch.length1 += prefix.len + suffix.len;
    patch.length2 += prefix.len + suffix.len;
}

/// Determines how to handle Diffs in a patch.  Functions which create
/// the diffs internally can use `.own`: the Diffs will be copied to
/// the patch list, new ones allocated, and old ones freed.  Then call
/// `deinit` on the DiffList, but not `deinitDiffList`.  This *must not*
/// be used if the DiffList is not immediately freed, because some of
/// the diffs will contain spuriously empty text.
///
/// Functions which operate on an existing DiffList should use `.copy`:
/// as the name indicates, copies of the Diffs will be made, and the
/// original memory must be freed separately.
const DiffHandling = enum {
    copy,
    own,
};

pub fn makePatchFromTexts(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    text1: []const u8,
    text2: []const u8,
) !PatchList {
    var diffs = try dmp.diff(allocator, text1, text2, true);
    // TODO try this with transfering patches once leaks are sorted out
    defer deinitDiffList(allocator, &diffs);
    if (diffs.items.len > 2) {
        try diffCleanupSemantic(allocator, &diffs);
        try dmp.diffCleanupEfficiency(allocator, &diffs);
    } // XXX TODO this should use .own once memory issues are sorted out
    return try dmp.makePatchInternal(allocator, text1, diffs, .copy);
}

/// @return List of Patch objects.
fn makePatchInternal(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    text: []const u8,
    diffs: DiffList,
    diff_act: DiffHandling,
) !PatchList {
    var patches = PatchList{};
    errdefer deinitPatchList(allocator, &patches);
    if (diffs.items.len == 0) {
        return patches; // Empty diff means empty patchlist
    }
    _ = diff_act; // XXX figure this out later

    var patch = Patch{};
    errdefer patch.deinit(allocator);
    var char_count1: usize = 0;
    var char_count2: usize = 0;
    // This avoids freeing the original copy of the text:
    var first_patch = true;
    var prepatch_text = text;
    defer {
        if (!first_patch)
            allocator.free(prepatch_text);
    }
    var postpatch = try std.ArrayList(u8).initCapacity(allocator, text.len);
    defer postpatch.deinit();
    try postpatch.appendSlice(text);
    for (diffs.items, 0..) |a_diff, i| {
        _ = i; // XXX
        if (patch.diffs.items.len == 0 and a_diff.operation != .equal) {
            patch.start1 = char_count1;
            patch.start2 = char_count2;
        }
        switch (a_diff.operation) {
            .insert => {
                try patch.diffs.ensureUnusedCapacity(allocator, 1);
                const d = try a_diff.clone(allocator); // if (diff_act == .copy) try a_diff.clone(allocator) else a_diff;
                patch.diffs.appendAssumeCapacity(d);
                patch.length2 += a_diff.text.len;
                try postpatch.insertSlice(char_count2, a_diff.text);
            },
            .delete => {
                try patch.diffs.ensureUnusedCapacity(allocator, 1);
                const d = try a_diff.clone(allocator); // if (diff_act == .copy) try a_diff.clone(allocator) else a_diff;
                patch.diffs.appendAssumeCapacity(d);
                patch.length1 += a_diff.text.len;
                try postpatch.replaceRange(char_count2, a_diff.text.len, "");
            },
            .equal => {
                //
                if (a_diff.text.len <= 2 * dmp.patch_margin and patch.diffs.items.len != 0 and !a_diff.eql(diffs.getLast())) {
                    // Small equality inside a patch.
                    try patch.diffs.ensureUnusedCapacity(allocator, 1);
                    const d = try a_diff.clone(allocator); // if (diff_act == .copy) try a_diff.clone(allocator) else a_diff;
                    patch.diffs.appendAssumeCapacity(d);
                    patch.length1 += a_diff.text.len;
                    patch.length2 += a_diff.text.len;
                }
                if (a_diff.text.len >= 2 * dmp.patch_margin) {
                    // Time for a new patch.
                    if (patch.diffs.items.len != 0) {
                        // Free the Diff if we own it.
                        //  if (diff_act == .own) {
                        //      allocator.free(a_diff.text);
                        //      // Replace with null patch to prevent double-free on error
                        //      diffs.items[i] = Diff{ .operation = .equal, .text = "" };
                        //  }
                        try dmp.patchAddContext(allocator, &patch, prepatch_text);
                        try patches.ensureUnusedCapacity(allocator, 1);
                        patches.appendAssumeCapacity(patch);
                        patch = Patch{};
                        // Unlike Unidiff, our patch lists have a rolling context.
                        // https://github.com/google/diff-match-patch/wiki/Unidiff
                        // Update prepatch text & pos to reflect the application of the
                        // just completed patch.
                        if (first_patch) {
                            // no free on first, we don't own the original text
                            first_patch = false;
                        } else {
                            allocator.free(prepatch_text);
                        }
                        prepatch_text = try allocator.dupe(u8, postpatch.items);
                        char_count1 = char_count2;
                    }
                }
            },
        }
        // Update the current character count.
        if (a_diff.operation != .insert) {
            char_count1 += a_diff.text.len;
        }
        if (a_diff.operation != .delete) {
            char_count2 += a_diff.text.len;
        }
    } // end for loop

    // Pick up the leftover patch if not empty.
    if (patch.diffs.items.len != 0) {
        try dmp.patchAddContext(allocator, &patch, prepatch_text);
        try patches.ensureUnusedCapacity(allocator, 1);
        patches.appendAssumeCapacity(patch);
    }
    return patches;
}

/// Compute a list of patches to turn text1 into text2.
/// text2 is not provided, diffs are the delta between text1 and text2.
///
/// @param text1 Old text.
/// @param diffs Array of Diff objects for text1 to text2.
pub fn makePatch(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    text: []const u8,
    diffs: DiffList,
) !PatchList {
    try dmp.makePatchInternal(allocator, text, diffs, .copy);
}

pub fn makePatchFromDiffs(dmp: DiffMatchPatch, allocator: Allocator, diffs: DiffList) !PatchList {
    const text1 = try diffBeforeText(allocator, diffs);
    defer allocator.free(text1);
    return try dmp.makePatch(allocator, text1, diffs, .copy);
}

/// Merge a set of patches onto the text.  Returns a tuple: the first of which
/// is the patched text, the second of which is...
///
/// TODO I'm just going to return a boolean saying whether all patches
/// were successful.  Rethink this at some point.
///
/// @param patches Array of Patch objects
/// @param text Old text.
/// @return Two element Object array, containing the new text and an array of
///      bool values.
pub fn patchApply(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    og_patches: PatchList,
    og_text: []const u8,
) !struct { []const u8, bool } {
    if (og_patches.items.len == 0) {
        // As silly as this is, we dupe the text, because something
        // passing an empty patchset isn't going to check, and will
        // end up double-freeing if we don't.  Going with 'true' as
        // the null patchset was successfully 'applied' here.
        return .{ try allocator.dupe(u8, og_text), true };
    }
    // So we can report if all patches were applied:
    var all_applied = true;
    // Deep copy the patches so that no changes are made to originals.
    const patches = try patchListClone(allocator, og_patches);
    defer patches.deinit(allocator);
    const null_padding = try patchAddPadding(allocator, patches);
    var text_array = try std.ArrayList(u8).initCapacity(og_text.len);
    defer text_array.deinit();
    try text_array.appendSlice(null_padding);
    try text_array.appendSlice(og_text);
    try text_array.appendSlice(null_padding);
    try patchSplitMax(allocator, patches);
    // delta keeps track of the offset between the expected and actual
    // location of the previous patch.  If there are patches expected at
    // positions 10 and 20, but the first patch was found at 12, delta is 2
    // and the second patch has an effective expected position of 22.
    var delta: usize = 0;
    for (patches) |a_patch| {
        const expected_loc = a_patch.start2 + delta;
        const text1 = try diffBeforeText(allocator, a_patch.diffs);
        defer allocator.free(text1);
        var maybe_start: ?usize = null;
        var maybe_end: ?usize = null;
        const m_max_b = dmp.match_max_bits;
        if (text1.len > m_max_b) {
            // patchSplitMax will only provide an oversized pattern
            // in the case of a monster delete.
            maybe_start = dmp.matchMain(
                allocator,
                text_array.items,
                text1[0..m_max_b],
                expected_loc,
            );
            if (maybe_start) |start| {
                const e_start = text1.len - m_max_b;
                maybe_end = dmp.matchMain(
                    allocator,
                    text_array.items,
                    text1[e_start..],
                    e_start + expected_loc,
                );
                // No match if a) no end_loc or b) the matches cross each other.
                if (maybe_end) |end| {
                    if (start >= end) {
                        maybe_start = null;
                    }
                } else {
                    maybe_start = null;
                }
            }
        } else {
            maybe_start = dmp.matchMain(allocator, og_text, text1, expected_loc);
        }
        if (maybe_start) |start| {
            // Found a match.  :)
            delta = start - expected_loc;
            // results[x] = true;
            const text2 = t2: {
                if (maybe_end) |end| {
                    break :t2 og_text[start..@min(end + m_max_b, og_text.len)];
                } else {
                    break :t2 og_text[start..@min(start + text1.len, og_text.len)];
                }
            };
            if (std.mem.eql(u8, text1, text2)) {
                // Perfect match, just shove the replacement text in.
                const diff_text = try diffAfterText(allocator, a_patch.diffs);
                defer allocator.free(diff_text);
                try text_array.replaceRange(start, text1.len, diff_text);
            } else {
                // Imperfect match.  Run a diff to get a framework of equivalent
                // indices.
                const diffs = try dmp.diff(
                    allocator,
                    text1,
                    text2,
                    false,
                );
                const t1_l_float: f64 = @floatFromInt(text1.len);
                const bad_match = diffLevenshtein(diffs) / t1_l_float > @This().patch_delete_threshold;
                if (text1.len > m_max_b and bad_match) {
                    // The end points match, but the content is unacceptably bad.
                    // results[x] = false;
                    all_applied = false;
                } else {
                    diffCleanupSemanticLossless(allocator, diffs);
                    var index1: usize = 0;
                    for (diffs) |a_diff| {
                        if (a_diff.operation != .equal) {
                            const index2 = diffIndex(diffs, index1);
                            if (a_diff.operation == .insert) {
                                // Insertion
                                try text_array.insertSlice(start + index2, a_diff.text);
                            } else if (a_diff.operation == .delete) {
                                // Deletion
                                try text_array.replaceRange(
                                    start + index2,
                                    diffIndex(diffs, index1 + a_diff.text.len),
                                    .{},
                                );
                            }
                            if (a_diff.operation != .delete) {
                                index1 += a_diff.text.len;
                            }
                        }
                    }
                }
            }
        } else {
            // No match found.  :(
            all_applied = false;
            // Subtract the delta for this failed patch from subsequent patches.
            delta -= a_patch.length2 - a_patch.length1;
        }
    }
    // strip padding
    try text_array.replaceRange(0, null_padding.len, .{});
    text_array.items.len -= null_padding.len;
    return .{ text_array.toOwnedSlice(), all_applied };
}

// Look through the patches and break up any which are longer than the
// maximum limit of the match algorithm.
// Intended to be called only from within patchApply.
// @param patches List of Patch objects.
fn patchSplitMax(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    patches: PatchList,
) !PatchList {
    const patch_size = dmp.match_max_bits;
    const patch_margin = dmp.patch_margin;
    const max_patch_len = patch_size - patch_size - patch_margin;
    // Mutating an array while iterating it? Sure, lets!
    var x = 0;
    while (x < patches.len) : (x += 1) {
        if (patches[x].length1 <= patch_size) continue;
        // We have a big ol' patch.
        const bigpatch = patches.orderedRemove(x);
        defer bigpatch.deinit(allocator);
        // Prevent incrementing past the next patch:
        x -= 1;
        var start1 = bigpatch.start1;
        var start2 = bigpatch.start2;
        // start with an empty precontext so that we can deinit consistently
        var precontext = try allocator.alloc(u8, 0);
        while (bigpatch.diffs.items.len != 0) {
            // Create one of several smaller patches.
            var patch = Patch{};
            var empty = true;
            patch.start1 = start1 - precontext.items.len;
            patch.start2 = start2 - precontext.items.len;
            if (precontext.len != 0) {
                patch.length2 = precontext.length;
                patch.length1 = patch.length2;
                try patch.diffs.ensureUnusedCapacity(allocator, 1);
                patch.diffs.appendAssumeCapacity(
                    allocator,
                    Diff{
                        .operation = .equal,
                        .text = try precontext.toOwnedSlice(),
                    },
                );
            }
            while (bigpatch.diffs.count != 0 and patch.length1 < max_patch_len) {
                const diff_type = bigpatch.diffs[0].operation;
                const diff_text = bigpatch.diffs[0].text;
                if (diff_type == .insert) {
                    // Insertions are harmless.
                    patch.length2 += diff_text.len;
                    start2 += diff_text.len;
                    // Move the patch (transfers ownership)
                    const diff1 = bigpatch.diffs.orderedRemove(0);
                    patch.diffs.append(diff1);
                    empty = false;
                } else if (cond: {
                    // zig fmt simply will not line break if clauses :/
                    const a = diff_type == .delete;
                    const b = patch.diffs.items.len == 1;
                    const c = patch.diffs[0].operation == .equal;
                    const d = diff_text.len > 2 * patch_size;
                    break :cond a and b and c and d;
                }) {
                    // This is a large deletion.  Let it pass in one chunk.
                    patch.length1 += diff_text.len;
                    start1 += diff_text.len;
                    empty = false;
                    // Transfer to patch:
                    try patch.diffs.ensureUnusedCapacity(allocator, 1);
                    const diff1 = bigpatch.diffs.orderedRemove(0);
                    patch.diffs.appendAssumeCapacity(diff1);
                } else {
                    // Deletion or equality.  Only take as much as we can stomach.
                    const text_end = @min(diff_text.len, patch_size - patch.length1 - patch_margin);
                    const new_diff_text = diff_text[0..text_end];
                    patch.length += new_diff_text.len;
                    start1 += new_diff_text.len;
                    if (diff_type == .equal) {
                        patch.length2 += diff_text.len;
                        start2 += diff_text.len;
                    } else {
                        empty = false;
                    }
                    // Now check if we did anything.
                    if (new_diff_text.len == diff_text.len) {
                        // We can reuse the diff.
                        try patch.diffs.ensureUnusedCapacity(allocator, 1);
                        const diff1 = bigpatch.diffs.orderedRemove(0);
                        patch.diffs.append(diff1);
                    } else {
                        // Free and dupe
                        const old_diff = bigpatch.diffs[0];
                        errdefer old_diff.deinit(allocator);
                        bigpatch.diffs[0] = Diff{
                            .operation = diff_type,
                            .text = try allocator.dupe(u8, new_diff_text),
                        };
                        old_diff.deinit(allocator);
                    }
                }
            }
            // Compute the head context for the next patch.
            const context_len: isize = precontext.len - patch_margin;
            allocator.free(precontext);
            if (context_len > 0) {
                const after_text = try diffAfterText(allocator, patch.diffs);
                defer allocator.free(after_text);
                precontext = try allocator.dupe(u8, after_text[context_len..]);
            } else {
                precontext = try allocator.alloc(u8, 0);
            }
            // Append the end context for this patch.
            const post_text = try diffBeforeText(bigpatch.diffs);
            const postcontext = post: {
                if (post_text.len > patch_margin) {
                    defer allocator.free(post_text);
                    break :post post_text[0..patch_margin];
                } else {
                    break :post post_text;
                }
            };
            if (postcontext.len != 0) {
                patch.length1 += postcontext.len;
                patch.length2 += postcontext.len;
                const maybe_last_diff = patch.diffs.getLastOrNull();
                if (maybe_last_diff) |last_diff| {
                    if (last_diff.operation == .equal) {
                        // free this diff and swap in a new one
                        defer last_diff.deinit(allocator);
                        patch.diffs.items.len -= 1;
                        const new_diff_text = try std.mem.concat(
                            allocator,
                            last_diff.text,
                            postcontext,
                        );
                        try patch.diffs.append(
                            allocator,
                            Diff{ .operation = .equal, .text = new_diff_text },
                        );
                    }
                } else {
                    // New diff from postcontext.
                    try patch.diffs.append(
                        allocator,
                        Diff{ .operation = .equal, .text = postcontext },
                    );
                }
            }
            if (!empty) {
                // Insert the next patch
                // Goes after x, and we need increment to skip:
                x += 1;
                try patches.insert(allocator, x, patch);
            }
        }
        // Free final precontext.
        allocator.free(precontext);
    }
}

/// Add some padding on text start and end so that edges can match something.
/// Intended to be called only from within patchApply.
/// @param patches Array of Patch objects.
/// @return The padding string added to each side.
fn patchAddPadding(
    dmp: DiffMatchPatch,
    allocator: Allocator,
    patches: PatchList,
) ![]const u8 {
    assert(patches.items.len != 0);
    const pad_len = dmp.patch_margin;
    var paddingcodes = try std.ArrayList(u8).initCapacity(allocator, pad_len);
    defer paddingcodes.deinit();
    {
        var control_code: u8 = 1;
        while (control_code <= pad_len) : (control_code += 1) {
            try paddingcodes.append(control_code);
        }
    }
    // Bump all the patches forward.
    for (patches) |a_patch| {
        a_patch.start1 += pad_len;
        a_patch.start2 += pad_len;
    }
    // Add some padding on start of first diff.
    var patch = patches.items[0];
    var diffs = patch.diffs;
    if (diffs.items.len == 0 or diffs.items[0].operation != .equal) {
        // Add nullPadding equality.
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.insert(
            0,
            Diff{
                .operation = .equal,
                .text = try allocator.dupe(u8, paddingcodes.items),
            },
        );
        patch.start1 -= pad_len;
        // OG code says "Should be 0" but this statement is not justified...
        assert(patch.start1 == 0);
        patch.start2 -= pad_len;
        assert(patch.start2 == 0);
        patch.length1 += pad_len;
        patch.lenght2 += pad_len;
    } else if (pad_len > diffs.items[0].text.len) {
        // Grow first equality.
        var diff1 = diffs.items[0];
        defer allocator.free(diff1.text);
        const extra_len = pad_len - diff1.text.len;
        diff1.text = try std.mem.concat(
            allocator,
            paddingcodes.items[diff1.text.len..],
            diff1.text,
        );
        patch.start1 -= extra_len;
        patch.start2 -= extra_len;
        patch.length1 += extra_len;
        patch.length2 += extra_len;
    }
    // Add some padding on end of last diff.
    patch = patches.getLast();
    diffs = patch.diffs;
    if (diffs.items.len == 0 or diffs.getLast().opeation != .equal) {
        // Add nullPadding equality.
        try diffs.ensureUnusedCapacity(allocator, 1);
        diffs.appendAssumeCapacity(
            Diff{
                .operation = .equal,
                .text = try allocator.dupe(u8, paddingcodes.items),
            },
        );
        patch.length1 += pad_len;
        patch.length2 += pad_len;
    } else if (pad_len > diffs.getLast().text.len) {
        // Grow last equality.
        var last_diff = diffs.getLast();
        defer allocator.free(last_diff.text);
        const extra_len = pad_len - last_diff.text.len;
        last_diff.text = try std.mem.concat(
            allocator,
            last_diff.text,
            paddingcodes[0..extra_len],
        );
        patch.length1 += extra_len;
        patch.length2 += extra_len;
    }
    return paddingcodes.toOwnedSlice();
}

/// Given an array of patches, return another array that is identical.
/// @param patches Array of Patch objects.
/// @return Array of Patch objects.
fn patchListClone(allocator: Allocator, patches: PatchList) !PatchList {
    var new_patches = PatchList{};
    errdefer {
        for (new_patches) |p| {
            p.deinit(allocator);
        }
    }
    new_patches.initCapacity(allocator, patches.items.len);
    for (patches) |patch| {
        try new_patches.append(allocator, try patch.clone(allocator));
    }
    return new_patches;
}

/// Take a list of patches and return a textual representation.
/// @param patches List of Patch objects.
/// @return Text representation of patches.
pub fn patchToText(allocator: Allocator, patches: PatchList) ![]const u8 {
    var text_array = std.ArrayList(u8).init(allocator);
    defer text_array.deinit();
    const writer = text_array.writer();
    try writePatch(writer, patches);
    return text_array.toOwnedSlice();
}

/// Stream a `PatchList` to the provided Writer.
pub fn writePatch(writer: anytype, patches: PatchList) !void {
    for (patches.items) |a_patch| {
        try a_patch.writeText(writer);
    }
}

/// Parse a textual representation of patches and return a List of Patch
/// objects.
/// @param textline Text representation of patches.
/// @return List of Patch objects.
/// @throws ArgumentException If invalid input.
pub fn patchFromText(allocator: Allocator, text: []const u8) !PatchList {
    if (text.len == 0) return PatchList{};
    var patches = PatchList{};
    errdefer patches.deinit(allocator);
    var cursor: usize = 0;
    while (cursor < text.len) {
        // TODO catch BadPatchString here and print diagnostic
        try patches.ensureUnusedCapacity(allocator, 1);
        const cursor_delta, const patch = try patchFromHeader(allocator, text[cursor..]);
        cursor += cursor_delta;
        patches.appendAssumeCapacity(patch);
    }
    return patches;
}

fn countDigits(text: []const u8) usize {
    var idx: usize = 0;
    while (std.ascii.isDigit(text[idx])) : (idx += 1) {}
    return idx;
}

fn patchFromHeader(allocator: Allocator, text: []const u8) !struct { usize, Patch } {
    var patch = Patch{ .diffs = DiffList{} };
    errdefer patch.deinit(allocator);
    var cursor: usize = undefined;
    if (std.mem.eql(u8, text[0..4], PATCH_HEAD)) {
        // Parse location and length in before text
        const count = 4 + countDigits(text[4..]);
        patch.start1 = std.fmt.parseInt(
            usize,
            text[4..count],
            10,
        ) catch return error.BadPatchString;
        cursor = count;
        assert(cursor > 4);
        if (text[cursor] != ',') {
            patch.start1 -= 1;
            patch.length1 = 1;
        } else {
            cursor += 1;
            const delta = countDigits(text[cursor..]);
            patch.length1 = std.fmt.parseInt(
                usize,
                text[cursor .. cursor + delta],
                10,
            ) catch return error.BadPatchString;
            assert(delta > 0);
            cursor += delta;
            if (patch.length1 != 0) {
                patch.start1 -= 1;
            }
        }
    } else return error.BadPatchString;
    // Parse location and length in after text.
    if (text[cursor] == ' ' and text[cursor + 1] == '+') {
        cursor += 2;
        const delta1 = countDigits(text[cursor..]);
        assert(delta1 > 0);
        patch.start2 = std.fmt.parseInt(
            usize,
            text[cursor .. cursor + delta1],
            10,
        ) catch return error.BadPatchString;
        cursor += delta1;
        if (text[cursor] != ',') {
            patch.start2 -= 1;
            patch.length2 = 1;
        } else {
            cursor += 1;
            const delta2 = countDigits(text[cursor..]);
            assert(delta2 > 0);
            patch.length2 = std.fmt.parseInt(
                usize,
                text[cursor .. cursor + delta2],
                10,
            ) catch return error.BadPatchString;
            cursor += delta2;
            if (patch.length2 != 0) {
                patch.start2 -= 1;
            }
        }
    } else return error.BadPatchString;
    if (std.mem.eql(u8, text[cursor .. cursor + 4], PATCH_TAIL)) {
        cursor += 4;
    } else return error.BadPatchString;
    // Eat the diffs
    var patch_lines = std.mem.splitScalar(
        u8,
        text[cursor..],
        '\n',
    );
    // `splitScalar` means blank lines, but we need that to
    // track the cursor.
    patch_loop: while (patch_lines.next()) |line| {
        cursor += line.len + 1;
        if (line.len == 0) continue;
        // Microsoft encodes spaces as +, we don't, so we don't need this:
        // line = line.Replace("+", "%2b");
        const diff_line = decodeUri(allocator, line[1..]) catch |e| {
            switch (e) {
                error.OutOfMemory => return e,
                else => return error.BadPatchString,
            }
        };
        errdefer allocator.free(diff_line);
        switch (line[0]) {
            '+' => { // Insertion
                try patch.diffs.append(
                    allocator,
                    Diff{
                        .operation = .insert,
                        .text = diff_line,
                    },
                );
            },
            '-' => { // Deletion
                try patch.diffs.append(
                    allocator,
                    Diff{
                        .operation = .delete,
                        .text = diff_line,
                    },
                );
            },
            ' ' => { // Minor equality
                try patch.diffs.append(
                    allocator,
                    Diff{
                        .operation = .equal,
                        .text = diff_line,
                    },
                );
            },
            '@' => { // Start of next patch
                // back out cursor
                cursor -= line.len - 1;
                break :patch_loop;
            },
            else => return error.BadPatchString,
        }
    } // end while
    return .{ cursor, patch };
}

/// Decode our URI-esque escaping
fn decodeUri(allocator: Allocator, line: []const u8) ![]const u8 {
    if (std.mem.indexOf(u8, line, "%")) |first| {
        // Text to decode.
        // Result will always be shorter than line:
        var new_line = try std.ArrayList(u8).initCapacity(allocator, line.len);
        defer new_line.deinit();
        try new_line.appendSlice(line[0..first]);
        var out_buf: [1]u8 = .{0};
        var codeunit = std.fmt.hexToBytes(
            &out_buf,
            line[first + 1 .. first + 3],
        ) catch return error.BadPatchString;
        try new_line.append(codeunit[0]);
        var cursor = first + 3;
        while (std.mem.indexOf(u8, line[cursor..], "%")) |next| {
            codeunit = std.fmt.hexToBytes(
                &out_buf,
                line[next + 1 .. next + 3],
            ) catch return error.BadPatchString;
            try new_line.append(codeunit[0]);
            cursor = next + 3;
        } else {
            try new_line.appendSlice(line[cursor..]);
        }
        return new_line.toOwnedSlice();
    } else {
        return allocator.dupe(u8, line);
    }
}

///
/// Borrowed from https://github.com/elerch/aws-sdk-for-zig/blob/master/src/aws_http.zig
/// under the MIT license. Thanks!
///
/// Modified to implement Unidiff escaping, documented here:
/// https://github.com/google/diff-match-patch/wiki/Unidiff
///
/// The documentation reads:
///
/// > Special characters are encoded using %xx notation. The set of
/// > characters which are encoded matches JavaScript's `encodeURI()`
/// > function, with the exception of spaces which are not encoded.
///
/// So we encode everything but the characters defined by Moz:
/// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/encodeURI
///
/// These:  !#$&'()*+,-./:;=?@_~  (and alphanumeric ASCII)
///
/// There is a nice contiguous run of 10 symbols between `&` and `/`, which we
/// can test in two comparisons, leaving these assorted:
///
///     !#$:;=?@_~
///
/// Each URI encoded byte is formed by a '%' and the two-digit
/// hexadecimal value of the byte.
///
/// Letters in the hexadecimal value must be uppercase, for example "%1A".
///
fn writeUriEncoded(writer: anytype, text: []const u8) !usize {
    const remaining_characters = "!#$:;=?@_~";
    var written: usize = 0;
    for (text) |c| {
        const should_encode = should: {
            if (c == ' ' or std.ascii.isAlphanumeric(c)) {
                break :should false;
            }
            if ('&' <= c and c <= '/') {
                break :should false;
            }
            for (remaining_characters) |r| {
                if (r == c) {
                    break :should false;
                }
            }
            break :should true;
        };

        if (!should_encode) {
            try writer.writeByte(c);
            written += 1;
            continue;
        }
        // Whatever remains, encode it
        try writer.writeByte('%');
        written += 1;
        const hexen = std.fmt.bytesToHex(&[_]u8{c}, .upper);
        written += try writer.write(&hexen);
    }
    return written;
}

fn encodeUri(allocator: std.mem.Allocator, text: []const u8) ![]u8 {
    var charlist = try std.ArrayList(u8).initCapacity(allocator, text.len);
    defer charlist.deinit();
    const writer = charlist.writer();
    _ = try writeUriEncoded(writer, text);
    return charlist.toOwnedSlice();
}

test encodeUri {
    const allocator = std.testing.allocator;
    const special_chars = "!#$&'()*+,-./:;=?@_~";
    const special_encoded = try encodeUri(allocator, special_chars);
    defer allocator.free(special_encoded);
    try testing.expectEqualStrings(special_chars, special_encoded);
    const alphaspace = " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    const alpha_encoded = try encodeUri(allocator, alphaspace);
    defer allocator.free(alpha_encoded);
    try testing.expectEqualStrings(alphaspace, alpha_encoded);
    const to_encode = "\"%<>[\\]^`{|}δ";
    const encodes = try encodeUri(allocator, to_encode);
    defer allocator.free(encodes);
    try testing.expectEqualStrings("%22%25%3C%3E%5B%5C%5D%5E%60%7B%7C%7D%CE%B4", encodes);
}

test diffCommonPrefix {
    // Detect any common suffix.
    try testing.expectEqual(@as(usize, 0), diffCommonPrefix("abc", "xyz")); // Null case
    try testing.expectEqual(@as(usize, 4), diffCommonPrefix("1234abcdef", "1234xyz")); // Non-null case
    try testing.expectEqual(@as(usize, 4), diffCommonPrefix("1234", "1234xyz")); // Whole case
}

test diffCommonSuffix {
    // Detect any common suffix.
    try testing.expectEqual(@as(usize, 0), diffCommonSuffix("abc", "xyz")); // Null case
    try testing.expectEqual(@as(usize, 4), diffCommonSuffix("abcdef1234", "xyz1234")); // Non-null case
    try testing.expectEqual(@as(usize, 4), diffCommonSuffix("1234", "xyz1234")); // Whole case
}

test diffCommonOverlap {
    // Detect any suffix/prefix overlap.
    try testing.expectEqual(@as(usize, 0), diffCommonOverlap("", "abcd")); // Null case
    try testing.expectEqual(@as(usize, 3), diffCommonOverlap("abc", "abcd")); // Whole case
    try testing.expectEqual(@as(usize, 0), diffCommonOverlap("123456", "abcd")); // No overlap
    try testing.expectEqual(@as(usize, 3), diffCommonOverlap("123456xxx", "xxxabcd")); // Overlap

    // Some overly clever languages (C#) may treat ligatures as equal to their
    // component letters.  E.g. U+FB01 == 'fi'
    try testing.expectEqual(@as(usize, 0), diffCommonOverlap("fi", "\u{fb01}")); // Unicode
}

fn testDiffHalfMatch(
    allocator: std.mem.Allocator,
    params: struct {
        dmp: DiffMatchPatch,
        before: []const u8,
        after: []const u8,
        expected: ?HalfMatchResult,
    },
) !void {
    const maybe_result = try params.dmp.diffHalfMatch(allocator, params.before, params.after);
    defer if (maybe_result) |result| result.deinit(allocator);
    try testing.expectEqualDeep(params.expected, maybe_result);
}

fn testDiffHalfMatchLeak(allocator: Allocator) !void {
    const dmp = DiffMatchPatch{};
    const text1 = "The quick brown fox jumps over the lazy dog.";
    const text2 = "That quick brown fox jumped over a lazy dog.";
    var diffs = try dmp.diff(allocator, text2, text1, true);
    deinitDiffList(allocator, &diffs);
}

test "diffHalfMatch leak regression test" {
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatchLeak, .{});
}

test diffHalfMatch {
    const one_timeout: DiffMatchPatch = .{ .diff_timeout = 1 };

    // No match #1
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "1234567890",
        .after = "abcdef",
        .expected = null,
    }});

    // No match #2
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "12345",
        .after = "23",
        .expected = null,
    }});

    // Single matches
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "1234567890",
        .after = "a345678z",
        .expected = .{
            .prefix_before = "12",
            .suffix_before = "90",
            .prefix_after = "a",
            .suffix_after = "z",
            .common_middle = "345678",
        },
    }});

    // Single Match #2
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "a345678z",
        .after = "1234567890",
        .expected = .{
            .prefix_before = "a",
            .suffix_before = "z",
            .prefix_after = "12",
            .suffix_after = "90",
            .common_middle = "345678",
        },
    }});

    // Single Match #3
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "abc56789z",
        .after = "1234567890",
        .expected = .{
            .prefix_before = "abc",
            .suffix_before = "z",
            .prefix_after = "1234",
            .suffix_after = "0",
            .common_middle = "56789",
        },
    }});

    // Single Match #4
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "a23456xyz",
        .after = "1234567890",
        .expected = .{
            .prefix_before = "a",
            .suffix_before = "xyz",
            .prefix_after = "1",
            .suffix_after = "7890",
            .common_middle = "23456",
        },
    }});

    // Multiple matches #1
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "121231234123451234123121",
        .after = "a1234123451234z",
        .expected = .{
            .prefix_before = "12123",
            .suffix_before = "123121",
            .prefix_after = "a",
            .suffix_after = "z",
            .common_middle = "1234123451234",
        },
    }});

    // Multiple Matches #2
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "x-=-=-=-=-=-=-=-=-=-=-=-=",
        .after = "xx-=-=-=-=-=-=-=",
        .expected = .{
            .prefix_before = "",
            .suffix_before = "-=-=-=-=-=",
            .prefix_after = "x",
            .suffix_after = "",
            .common_middle = "x-=-=-=-=-=-=-=",
        },
    }});

    // Multiple Matches #3
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "-=-=-=-=-=-=-=-=-=-=-=-=y",
        .after = "-=-=-=-=-=-=-=yy",
        .expected = .{
            .prefix_before = "-=-=-=-=-=",
            .suffix_before = "",
            .prefix_after = "",
            .suffix_after = "y",
            .common_middle = "-=-=-=-=-=-=-=y",
        },
    }});

    // Other cases

    // Optimal diff would be -q+x=H-i+e=lloHe+Hu=llo-Hew+y not -qHillo+x=HelloHe-w+Hulloy
    // Non-optimal halfmatch
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = one_timeout,
        .before = "qHilloHelloHew",
        .after = "xHelloHeHulloy",
        .expected = .{
            .prefix_before = "qHillo",
            .suffix_before = "w",
            .prefix_after = "x",
            .suffix_after = "Hulloy",
            .common_middle = "HelloHe",
        },
    }});

    // Non-optimal halfmatch
    try testing.checkAllAllocationFailures(testing.allocator, testDiffHalfMatch, .{.{
        .dmp = .{ .diff_timeout = 0 },
        .before = "qHilloHelloHew",
        .after = "xHelloHeHulloy",
        .expected = null,
    }});
}

test diffLinesToChars {
    const allocator = testing.allocator;
    // Convert lines down to characters.
    var tmp_array_list = std.ArrayList([]const u8).init(allocator);
    defer tmp_array_list.deinit();
    try tmp_array_list.append("");
    try tmp_array_list.append("alpha\n");
    try tmp_array_list.append("beta\n");

    var result = try diffLinesToChars(allocator, "alpha\nbeta\nalpha\n", "beta\nalpha\nbeta\n");
    try testing.expectEqualStrings("\u{0001}\u{0002}\u{0001}", result.chars_1); // Shared lines #1
    try testing.expectEqualStrings("\u{0002}\u{0001}\u{0002}", result.chars_2); // Shared lines #2
    try testing.expectEqualDeep(tmp_array_list.items, result.line_array.items); // Shared lines #3

    tmp_array_list.items.len = 0;
    try tmp_array_list.append("");
    try tmp_array_list.append("alpha\r\n");
    try tmp_array_list.append("beta\r\n");
    try tmp_array_list.append("\r\n");
    result.deinit(allocator);

    result = try diffLinesToChars(allocator, "", "alpha\r\nbeta\r\n\r\n\r\n");
    try testing.expectEqualStrings("", result.chars_1); // Empty string and blank lines #1
    try testing.expectEqualStrings("\u{0001}\u{0002}\u{0003}\u{0003}", result.chars_2); // Empty string and blank lines #2
    try testing.expectEqualDeep(tmp_array_list.items, result.line_array.items); // Empty string and blank lines #3

    tmp_array_list.items.len = 0;
    try tmp_array_list.append("");
    try tmp_array_list.append("a");
    try tmp_array_list.append("b");
    result.deinit(allocator);

    result = try diffLinesToChars(allocator, "a", "b");
    try testing.expectEqualStrings("\u{0001}", result.chars_1); // No linebreaks #1.
    try testing.expectEqualStrings("\u{0002}", result.chars_2); // No linebreaks #2.
    try testing.expectEqualDeep(tmp_array_list.items, result.line_array.items); // No linebreaks #3.
    result.deinit(allocator);

    // TODO: More than 256 to reveal any 8-bit limitations but this requires
    // some unicode logic that I don't want to deal with
    //
    // Casting to Unicode is straightforward and should sort correctly, I'm
    // more concerned about the weird behavior when the 'char' is equal to a
    // newline.  Uncomment the EqualSlices below to see what I mean.
    // I think there's some cleanup logic in the actual linediff that should
    // take care of the problem, but I don't like it.

    const n: u8 = 255;
    tmp_array_list.items.len = 0;

    var line_list = std.ArrayList(u8).init(allocator);
    defer line_list.deinit();
    var char_list = std.ArrayList(u8).init(allocator);
    defer char_list.deinit();

    var i: u8 = 1;
    while (i < n) : (i += 1) {
        try tmp_array_list.append(&.{ i, '\n' });
        try line_list.appendSlice(&.{ i, '\n' });
        try char_list.append(i);
    }
    try testing.expectEqual(@as(usize, n - 1), tmp_array_list.items.len); // Test initialization fail #1
    try testing.expectEqual(@as(usize, n - 1), char_list.items.len); // Test initialization fail #2
    try tmp_array_list.insert(0, "");
    result = try diffLinesToChars(allocator, line_list.items, "");
    defer result.deinit(allocator);
    // TODO: This isn't equal, should it be?
    // try testing.expectEqualSlices(u8, char_list.items, result.chars_1);
    try testing.expectEqualStrings("", result.chars_2);
    // TODO this is wrong because of the max_value I think?
    // try testing.expectEqualDeep(tmp_array_list.items, result.line_array.items);
}

fn testDiffCharsToLines(
    allocator: std.mem.Allocator,
    params: struct {
        diffs: []const Diff,
        line_array: []const []const u8,
        expected: []const Diff,
    },
) !void {
    var diffs = try DiffList.initCapacity(allocator, params.diffs.len);
    defer deinitDiffList(allocator, &diffs);

    for (params.diffs) |item| {
        diffs.appendAssumeCapacity(.{ .operation = item.operation, .text = try allocator.dupe(u8, item.text) });
    }

    try diffCharsToLines(allocator, diffs.items, params.line_array);

    try testing.expectEqualDeep(params.expected, diffs.items);
}

test diffCharsToLines {
    // Convert chars up to lines.
    var diff_list = DiffList{};
    defer deinitDiffList(testing.allocator, &diff_list);
    try diff_list.ensureTotalCapacity(testing.allocator, 2);
    diff_list.appendSliceAssumeCapacity(&.{
        Diff.init(.equal, try testing.allocator.dupe(u8, "\u{0001}\u{0002}\u{0001}")),
        Diff.init(.insert, try testing.allocator.dupe(u8, "\u{0002}\u{0001}\u{0002}")),
    });
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCharsToLines, .{.{
        .diffs = diff_list.items,
        .line_array = &[_][]const u8{
            "",
            "alpha\n",
            "beta\n",
        },
        .expected = &.{
            .{ .operation = .equal, .text = "alpha\nbeta\nalpha\n" },
            .{ .operation = .insert, .text = "beta\nalpha\nbeta\n" },
        },
    }});

    // TODO: Implement exhaustive tests
}

fn testDiffCleanupMerge(allocator: std.mem.Allocator, params: struct {
    input: []const Diff,
    expected: []const Diff,
}) !void {
    var diffs = try DiffList.initCapacity(allocator, params.input.len);
    defer deinitDiffList(allocator, &diffs);

    for (params.input) |item| {
        diffs.appendAssumeCapacity(.{ .operation = item.operation, .text = try allocator.dupe(u8, item.text) });
    }

    try diffCleanupMerge(allocator, &diffs);

    try testing.expectEqualDeep(params.expected, diffs.items);
}

test diffCleanupMerge {
    // Cleanup a messy diff.

    // No change case
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "b" },
            .{ .operation = .insert, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "b" },
            .{ .operation = .insert, .text = "c" },
        },
    }});

    // Merge equalities
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .equal, .text = "b" },
            .{ .operation = .equal, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "abc" },
        },
    }});

    // Merge deletions
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .delete, .text = "b" },
            .{ .operation = .delete, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abc" },
        },
    }});

    // Merge insertions
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .insert, .text = "a" },
            .{ .operation = .insert, .text = "b" },
            .{ .operation = .insert, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .insert, .text = "abc" },
        },
    }});

    // Merge interweave
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .insert, .text = "b" },
            .{ .operation = .delete, .text = "c" },
            .{ .operation = .insert, .text = "d" },
            .{ .operation = .equal, .text = "e" },
            .{ .operation = .equal, .text = "f" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "ac" },
            .{ .operation = .insert, .text = "bd" },
            .{ .operation = .equal, .text = "ef" },
        },
    }});

    // Prefix and suffix detection
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .insert, .text = "abc" },
            .{ .operation = .delete, .text = "dc" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "d" },
            .{ .operation = .insert, .text = "b" },
            .{ .operation = .equal, .text = "c" },
        },
    }});

    // Prefix and suffix detection with equalities
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "x" },
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .insert, .text = "abc" },
            .{ .operation = .delete, .text = "dc" },
            .{ .operation = .equal, .text = "y" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "xa" },
            .{ .operation = .delete, .text = "d" },
            .{ .operation = .insert, .text = "b" },
            .{ .operation = .equal, .text = "cy" },
        },
    }});

    // Slide edit left
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .insert, .text = "ba" },
            .{ .operation = .equal, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .insert, .text = "ab" },
            .{ .operation = .equal, .text = "ac" },
        },
    }});

    // Slide edit right
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "c" },
            .{ .operation = .insert, .text = "ab" },
            .{ .operation = .equal, .text = "a" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "ca" },
            .{ .operation = .insert, .text = "ba" },
        },
    }});

    // Slide edit left recursive
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "b" },
            .{ .operation = .equal, .text = "c" },
            .{ .operation = .delete, .text = "ac" },
            .{ .operation = .equal, .text = "x" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abc" },
            .{ .operation = .equal, .text = "acx" },
        },
    }});

    // Slide edit right recursive
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "x" },
            .{ .operation = .delete, .text = "ca" },
            .{ .operation = .equal, .text = "c" },
            .{ .operation = .delete, .text = "b" },
            .{ .operation = .equal, .text = "a" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "xca" },
            .{ .operation = .delete, .text = "cba" },
        },
    }});

    // Empty merge
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "b" },
            .{ .operation = .insert, .text = "ab" },
            .{ .operation = .equal, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .insert, .text = "a" },
            .{ .operation = .equal, .text = "bc" },
        },
    }});

    // Empty equality
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupMerge, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "" },
            .{ .operation = .insert, .text = "a" },
            .{ .operation = .equal, .text = "b" },
        },
        .expected = &.{
            .{ .operation = .insert, .text = "a" },
            .{ .operation = .equal, .text = "b" },
        },
    }});
}

fn testDiffCleanupSemanticLossless(
    allocator: std.mem.Allocator,
    params: struct {
        input: []const Diff,
        expected: []const Diff,
    },
) !void {
    var diffs = try DiffList.initCapacity(allocator, params.input.len);
    defer deinitDiffList(allocator, &diffs);

    for (params.input) |item| {
        diffs.appendAssumeCapacity(.{ .operation = item.operation, .text = try allocator.dupe(u8, item.text) });
    }

    try diffCleanupSemanticLossless(allocator, &diffs);

    try testing.expectEqualDeep(params.expected, diffs.items);
}

fn sliceToDiffList(allocator: Allocator, diff_slice: []const Diff) !DiffList {
    var diff_list = DiffList{};
    errdefer deinitDiffList(allocator, &diff_list);
    try diff_list.ensureTotalCapacity(allocator, diff_slice.len);
    for (diff_slice) |d| {
        diff_list.appendAssumeCapacity(Diff.init(
            d.operation,
            try allocator.dupe(u8, d.text),
        ));
    }
    return diff_list;
}

test diffCleanupSemanticLossless {
    // Null case
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &[_]Diff{},
        .expected = &[_]Diff{},
    }});

    //defer deinitDiffList(allocator, &diffs);
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "AAA\r\n\r\nBBB" },
            .{ .operation = .insert, .text = "\r\nDDD\r\n\r\nBBB" },
            .{ .operation = .equal, .text = "\r\nEEE" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "AAA\r\n\r\n" },
            .{ .operation = .insert, .text = "BBB\r\nDDD\r\n\r\n" },
            .{ .operation = .equal, .text = "BBB\r\nEEE" },
        },
    }});

    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "AAA\r\nBBB" },
            .{ .operation = .insert, .text = " DDD\r\nBBB" },
            .{ .operation = .equal, .text = " EEE" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "AAA\r\n" },
            .{ .operation = .insert, .text = "BBB DDD\r\n" },
            .{ .operation = .equal, .text = "BBB EEE" },
        },
    }});

    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "The c" },
            .{ .operation = .insert, .text = "ow and the c" },
            .{ .operation = .equal, .text = "at." },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "The " },
            .{ .operation = .insert, .text = "cow and the " },
            .{ .operation = .equal, .text = "cat." },
        },
    }});

    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "The-c" },
            .{ .operation = .insert, .text = "ow-and-the-c" },
            .{ .operation = .equal, .text = "at." },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "The-" },
            .{ .operation = .insert, .text = "cow-and-the-" },
            .{ .operation = .equal, .text = "cat." },
        },
    }});

    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .equal, .text = "ax" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .equal, .text = "aax" },
        },
    }});

    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "xa" },
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .equal, .text = "a" },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "xaa" },
            .{ .operation = .delete, .text = "a" },
        },
    }});

    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemanticLossless, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "The xxx. The " },
            .{ .operation = .insert, .text = "zzz. The " },
            .{ .operation = .equal, .text = "yyy." },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "The xxx." },
            .{ .operation = .insert, .text = " The zzz." },
            .{ .operation = .equal, .text = " The yyy." },
        },
    }});
}

fn rebuildtexts(allocator: std.mem.Allocator, diffs: DiffList) ![2][]const u8 {
    var text = [2]std.ArrayList(u8){
        std.ArrayList(u8).init(allocator),
        std.ArrayList(u8).init(allocator),
    };
    errdefer {
        text[0].deinit();
        text[1].deinit();
    }

    for (diffs.items) |myDiff| {
        if (myDiff.operation != .insert) {
            try text[0].appendSlice(myDiff.text);
        }
        if (myDiff.operation != .delete) {
            try text[1].appendSlice(myDiff.text);
        }
    }
    return .{
        try text[0].toOwnedSlice(),
        try text[1].toOwnedSlice(),
    };
}

fn testRebuildTexts(allocator: Allocator, diffs: DiffList, params: struct {
    before: []const u8,
    after: []const u8,
}) !void {
    const texts = try rebuildtexts(allocator, diffs);
    defer {
        allocator.free(texts[0]);
        allocator.free(texts[1]);
    }
    try testing.expectEqualStrings(params.before, texts[0]);
    try testing.expectEqualStrings(params.after, texts[1]);
}

test rebuildtexts {
    {
        var diffs = try sliceToDiffList(testing.allocator, &.{
            .{ .operation = .insert, .text = "abcabc" },
            .{ .operation = .equal, .text = "defdef" },
            .{ .operation = .delete, .text = "ghighi" },
        });
        defer deinitDiffList(testing.allocator, &diffs);
        try testing.checkAllAllocationFailures(testing.allocator, testRebuildTexts, .{
            diffs,
            .{
                .before = "defdefghighi",
                .after = "abcabcdefdef",
            },
        });
    }
    {
        var diffs = try sliceToDiffList(testing.allocator, &.{
            .{ .operation = .insert, .text = "xxx" },
            .{ .operation = .delete, .text = "yyy" },
        });
        defer deinitDiffList(testing.allocator, &diffs);
        try testing.checkAllAllocationFailures(testing.allocator, testRebuildTexts, .{
            diffs,
            .{
                .before = "yyy",
                .after = "xxx",
            },
        });
    }
    {
        var diffs = try sliceToDiffList(testing.allocator, &.{
            .{ .operation = .equal, .text = "xyz" },
            .{ .operation = .equal, .text = "pdq" },
        });
        defer deinitDiffList(testing.allocator, &diffs);
        try testing.checkAllAllocationFailures(testing.allocator, testRebuildTexts, .{
            diffs,
            .{
                .before = "xyzpdq",
                .after = "xyzpdq",
            },
        });
    }
}

fn testDiffBisect(
    allocator: std.mem.Allocator,
    params: struct {
        dmp: DiffMatchPatch,
        before: []const u8,
        after: []const u8,
        deadline: u64,
        expected: []const Diff,
    },
) !void {
    var diffs = try params.dmp.diffBisect(allocator, params.before, params.after, params.deadline);
    defer deinitDiffList(allocator, &diffs);
    try testing.expectEqualDeep(params.expected, diffs.items);
}

test diffBisect {
    const this: DiffMatchPatch = .{ .diff_timeout = 0 };

    const a = "cat";
    const b = "map";

    // Normal
    try testing.checkAllAllocationFailures(testing.allocator, testDiffBisect, .{.{
        .dmp = this,
        .before = a,
        .after = b,
        .deadline = std.math.maxInt(u64), // Travis TODO not sure if maxInt(u64) is correct for  DateTime.MaxValue
        .expected = &.{
            .{ .operation = .delete, .text = "c" },
            .{ .operation = .insert, .text = "m" },
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "t" },
            .{ .operation = .insert, .text = "p" },
        },
    }});

    // Timeout
    try testing.checkAllAllocationFailures(testing.allocator, testDiffBisect, .{.{
        .dmp = this,
        .before = a,
        .after = b,
        .deadline = 0, // Travis TODO not sure if 0 is correct for  DateTime.MinValue
        .expected = &.{
            .{ .operation = .delete, .text = "cat" },
            .{ .operation = .insert, .text = "map" },
        },
    }});
}

fn testDiff(
    allocator: std.mem.Allocator,
    params: struct {
        dmp: DiffMatchPatch,
        before: []const u8,
        after: []const u8,
        check_lines: bool,
        expected: []const Diff,
    },
) !void {
    var diffs = try params.dmp.diff(allocator, params.before, params.after, params.check_lines);
    defer deinitDiffList(allocator, &diffs);
    try testing.expectEqualDeep(params.expected, diffs.items);
}

test diff {
    const this: DiffMatchPatch = .{ .diff_timeout = 0 };

    //  Null case.
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "",
        .after = "",
        .check_lines = false,
        .expected = &[_]Diff{},
    }});

    //  Equality.
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "abc",
        .after = "abc",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .equal, .text = "abc" },
        },
    }});

    // Simple insertion.
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "abc",
        .after = "ab123c",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .equal, .text = "ab" },
            .{ .operation = .insert, .text = "123" },
            .{ .operation = .equal, .text = "c" },
        },
    }});

    // Simple deletion.
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "a123bc",
        .after = "abc",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "123" },
            .{ .operation = .equal, .text = "bc" },
        },
    }});

    // Two insertions.
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "abc",
        .after = "a123b456c",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .insert, .text = "123" },
            .{ .operation = .equal, .text = "b" },
            .{ .operation = .insert, .text = "456" },
            .{ .operation = .equal, .text = "c" },
        },
    }});

    // Two deletions.
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "a123b456c",
        .after = "abc",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "123" },
            .{ .operation = .equal, .text = "b" },
            .{ .operation = .delete, .text = "456" },
            .{ .operation = .equal, .text = "c" },
        },
    }});

    // Simple case #1
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "a",
        .after = "b",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .insert, .text = "b" },
        },
    }});

    // Simple case #2
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "Apples are a fruit.",
        .after = "Bananas are also fruit.",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .delete, .text = "Apple" },
            .{ .operation = .insert, .text = "Banana" },
            .{ .operation = .equal, .text = "s are a" },
            .{ .operation = .insert, .text = "lso" },
            .{ .operation = .equal, .text = " fruit." },
        },
    }});

    // Simple case #3
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "ax\t",
        .after = "\u{0680}x\x00",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .insert, .text = "\u{0680}" },
            .{ .operation = .equal, .text = "x" },
            .{ .operation = .delete, .text = "\t" },
            .{ .operation = .insert, .text = "\x00" },
        },
    }});

    // Overlap #1
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "1ayb2",
        .after = "abxab",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .delete, .text = "1" },
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "y" },
            .{ .operation = .equal, .text = "b" },
            .{ .operation = .delete, .text = "2" },
            .{ .operation = .insert, .text = "xab" },
        },
    }});

    // Overlap #2
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "abcy",
        .after = "xaxcxabc",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .insert, .text = "xaxcx" },
            .{ .operation = .equal, .text = "abc" },
            .{ .operation = .delete, .text = "y" },
        },
    }});

    // Overlap #3
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "ABCDa=bcd=efghijklmnopqrsEFGHIJKLMNOefg",
        .after = "a-bcd-efghijklmnopqrs",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .delete, .text = "ABCD" },
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .delete, .text = "=" },
            .{ .operation = .insert, .text = "-" },
            .{ .operation = .equal, .text = "bcd" },
            .{ .operation = .delete, .text = "=" },
            .{ .operation = .insert, .text = "-" },
            .{ .operation = .equal, .text = "efghijklmnopqrs" },
            .{ .operation = .delete, .text = "EFGHIJKLMNOefg" },
        },
    }});

    // Large equality
    try testing.checkAllAllocationFailures(testing.allocator, testDiff, .{.{
        .dmp = this,
        .before = "a [[Pennsylvania]] and [[New",
        .after = " and [[Pennsylvania]]",
        .check_lines = false,
        .expected = &.{
            .{ .operation = .insert, .text = " " },
            .{ .operation = .equal, .text = "a" },
            .{ .operation = .insert, .text = "nd" },
            .{ .operation = .equal, .text = " [[Pennsylvania]]" },
            .{ .operation = .delete, .text = " and [[New" },
        },
    }});

    const allocator = testing.allocator;
    // TODO these tests should be checked for allocation failure

    // Increase the text lengths by 1024 times to ensure a timeout.
    {
        const a = "`Twas brillig, and the slithy toves\nDid gyre and gimble in the wabe:\nAll mimsy were the borogoves,\nAnd the mome raths outgrabe.\n" ** 1024;
        const b = "I am the very model of a modern major general,\nI've information vegetable, animal, and mineral,\nI know the kings of England, and I quote the fights historical,\nFrom Marathon to Waterloo, in order categorical.\n" ** 1024;

        const with_timout: DiffMatchPatch = .{
            .diff_timeout = 100, // 100ms
        };

        const start_time = std.time.milliTimestamp();
        {
            var time_diff = try with_timout.diff(allocator, a, b, false);
            defer deinitDiffList(allocator, &time_diff);
        }
        const end_time = std.time.milliTimestamp();

        // Test that we took at least the timeout period.
        try testing.expect(with_timout.diff_timeout <= end_time - start_time); // diff: Timeout min.
        // Test that we didn't take forever (be forgiving).
        // Theoretically this test could fail very occasionally if the
        // OS task swaps or locks up for a second at the wrong moment.
        try testing.expect((with_timout.diff_timeout) * 10000 * 2 > end_time - start_time); // diff: Timeout max.
    }

    {
        // Test the linemode speedup.
        // Must be long to pass the 100 char cutoff.
        const a = "1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n";
        const b = "abcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\nabcdefghij\n";

        var diff_checked = try this.diff(allocator, a, b, true);
        defer deinitDiffList(allocator, &diff_checked);

        var diff_unchecked = try this.diff(allocator, a, b, false);
        defer deinitDiffList(allocator, &diff_unchecked);

        try testing.expectEqualDeep(diff_checked.items, diff_unchecked.items); // diff: Simple line-mode.
    }

    {
        const a = "1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890";
        const b = "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghij";

        var diff_checked = try this.diff(allocator, a, b, true);
        defer deinitDiffList(allocator, &diff_checked);

        var diff_unchecked = try this.diff(allocator, a, b, false);
        defer deinitDiffList(allocator, &diff_unchecked);

        try testing.expectEqualDeep(diff_checked.items, diff_unchecked.items); // diff: Single line-mode.
    }

    {
        // diff: Overlap line-mode.
        const a = "1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n1234567890\n";
        const b = "abcdefghij\n1234567890\n1234567890\n1234567890\nabcdefghij\n1234567890\n1234567890\n1234567890\nabcdefghij\n1234567890\n1234567890\n1234567890\nabcdefghij\n";

        var diffs_linemode = try this.diff(allocator, a, b, true);
        defer deinitDiffList(allocator, &diffs_linemode);

        const texts_linemode = try rebuildtexts(allocator, diffs_linemode);
        defer {
            allocator.free(texts_linemode[0]);
            allocator.free(texts_linemode[1]);
        }

        var diffs_textmode = try this.diff(allocator, a, b, false);
        defer deinitDiffList(allocator, &diffs_textmode);

        const texts_textmode = try rebuildtexts(allocator, diffs_textmode);
        defer {
            allocator.free(texts_textmode[0]);
            allocator.free(texts_textmode[1]);
        }

        try testing.expectEqualStrings(texts_textmode[0], texts_linemode[0]);
        try testing.expectEqualStrings(texts_textmode[1], texts_linemode[1]);
    }
}

fn testDiffCleanupSemantic(
    allocator: std.mem.Allocator,
    params: struct {
        input: []const Diff,
        expected: []const Diff,
    },
) !void {
    var diffs = try DiffList.initCapacity(allocator, params.input.len);
    defer deinitDiffList(allocator, &diffs);

    for (params.input) |item| {
        diffs.appendAssumeCapacity(.{ .operation = item.operation, .text = try allocator.dupe(u8, item.text) });
    }

    try diffCleanupSemantic(allocator, &diffs);

    try testing.expectEqualDeep(params.expected, diffs.items);
}

test diffCleanupSemantic {
    // Null case.
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &[_]Diff{},
        .expected = &[_]Diff{},
    }});

    // No elimination #1
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .insert, .text = "cd" },
            .{ .operation = .equal, .text = "12" },
            .{ .operation = .delete, .text = "e" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .insert, .text = "cd" },
            .{ .operation = .equal, .text = "12" },
            .{ .operation = .delete, .text = "e" },
        },
    }});

    // No elimination #2
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "abc" },
            .{ .operation = .insert, .text = "ABC" },
            .{ .operation = .equal, .text = "1234" },
            .{ .operation = .delete, .text = "wxyz" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abc" },
            .{ .operation = .insert, .text = "ABC" },
            .{ .operation = .equal, .text = "1234" },
            .{ .operation = .delete, .text = "wxyz" },
        },
    }});

    // Simple elimination
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "a" },
            .{ .operation = .equal, .text = "b" },
            .{ .operation = .delete, .text = "c" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abc" },
            .{ .operation = .insert, .text = "b" },
        },
    }});

    // Backpass elimination
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .equal, .text = "cd" },
            .{ .operation = .delete, .text = "e" },
            .{ .operation = .equal, .text = "f" },
            .{ .operation = .insert, .text = "g" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abcdef" },
            .{ .operation = .insert, .text = "cdfg" },
        },
    }});

    // Multiple elimination
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .insert, .text = "1" },
            .{ .operation = .equal, .text = "A" },
            .{ .operation = .delete, .text = "B" },
            .{ .operation = .insert, .text = "2" },
            .{ .operation = .equal, .text = "_" },
            .{ .operation = .insert, .text = "1" },
            .{ .operation = .equal, .text = "A" },
            .{ .operation = .delete, .text = "B" },
            .{ .operation = .insert, .text = "2" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "AB_AB" },
            .{ .operation = .insert, .text = "1A2_1A2" },
        },
    }});

    // Word boundaries
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .equal, .text = "The c" },
            .{ .operation = .delete, .text = "ow and the c" },
            .{ .operation = .equal, .text = "at." },
        },
        .expected = &.{
            .{ .operation = .equal, .text = "The " },
            .{ .operation = .delete, .text = "cow and the " },
            .{ .operation = .equal, .text = "cat." },
        },
    }});

    // No overlap elimination
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "abcxx" },
            .{ .operation = .insert, .text = "xxdef" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abcxx" },
            .{ .operation = .insert, .text = "xxdef" },
        },
    }});

    // Overlap elimination
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "abcxxx" },
            .{ .operation = .insert, .text = "xxxdef" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abc" },
            .{ .operation = .equal, .text = "xxx" },
            .{ .operation = .insert, .text = "def" },
        },
    }});

    // Reverse overlap elimination
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "xxxabc" },
            .{ .operation = .insert, .text = "defxxx" },
        },
        .expected = &.{
            .{ .operation = .insert, .text = "def" },
            .{ .operation = .equal, .text = "xxx" },
            .{ .operation = .delete, .text = "abc" },
        },
    }});

    // Two overlap eliminations
    try testing.checkAllAllocationFailures(testing.allocator, testDiffCleanupSemantic, .{.{
        .input = &.{
            .{ .operation = .delete, .text = "abcd1212" },
            .{ .operation = .insert, .text = "1212efghi" },
            .{ .operation = .equal, .text = "----" },
            .{ .operation = .delete, .text = "A3" },
            .{ .operation = .insert, .text = "3BC" },
        },
        .expected = &.{
            .{ .operation = .delete, .text = "abcd" },
            .{ .operation = .equal, .text = "1212" },
            .{ .operation = .insert, .text = "efghi" },
            .{ .operation = .equal, .text = "----" },
            .{ .operation = .delete, .text = "A" },
            .{ .operation = .equal, .text = "3" },
            .{ .operation = .insert, .text = "BC" },
        },
    }});
}

fn testDiffCleanupEfficiency(
    allocator: Allocator,
    dmp: DiffMatchPatch,
    params: struct {
        input: []const Diff,
        expected: []const Diff,
    },
) !void {
    var diffs = try DiffList.initCapacity(allocator, params.input.len);
    defer deinitDiffList(allocator, &diffs);
    for (params.input) |item| {
        diffs.appendAssumeCapacity(.{ .operation = item.operation, .text = try allocator.dupe(u8, item.text) });
    }
    try dmp.diffCleanupEfficiency(allocator, &diffs);

    try testing.expectEqualDeep(params.expected, diffs.items);
}

test "diffCleanupEfficiency" {
    const allocator = testing.allocator;
    var dmp = DiffMatchPatch{};
    dmp.diff_edit_cost = 4;
    { // Null case.
        var diffs = DiffList{};
        try dmp.diffCleanupEfficiency(allocator, &diffs);
        try testing.expectEqualDeep(DiffList{}, diffs);
    }
    { // No elimination.
        const dslice: []const Diff = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .insert, .text = "12" },
            .{ .operation = .equal, .text = "wxyz" },
            .{ .operation = .delete, .text = "cd" },
            .{ .operation = .insert, .text = "34" },
        };
        try testing.checkAllAllocationFailures(
            testing.allocator,
            testDiffCleanupEfficiency,
            .{
                dmp,
                .{ .input = dslice, .expected = dslice },
            },
        );
    }
    { // Four-edit elimination.
        const dslice: []const Diff = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .insert, .text = "12" },
            .{ .operation = .equal, .text = "xyz" },
            .{ .operation = .delete, .text = "cd" },
            .{ .operation = .insert, .text = "34" },
        };
        const d_after: []const Diff = &.{
            .{ .operation = .delete, .text = "abxyzcd" },
            .{ .operation = .insert, .text = "12xyz34" },
        };
        try testing.checkAllAllocationFailures(
            testing.allocator,
            testDiffCleanupEfficiency,
            .{
                dmp,
                .{ .input = dslice, .expected = d_after },
            },
        );
    }
    { // Three-edit elimination.
        const dslice: []const Diff = &.{
            .{ .operation = .insert, .text = "12" },
            .{ .operation = .equal, .text = "x" },
            .{ .operation = .delete, .text = "cd" },
            .{ .operation = .insert, .text = "34" },
        };
        const d_after: []const Diff = &.{
            .{ .operation = .delete, .text = "xcd" },
            .{ .operation = .insert, .text = "12x34" },
        };
        try testing.checkAllAllocationFailures(
            testing.allocator,
            testDiffCleanupEfficiency,
            .{
                dmp,
                .{ .input = dslice, .expected = d_after },
            },
        );
    }
    { // Backpass elimination.
        const dslice: []const Diff = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .insert, .text = "12" },
            .{ .operation = .equal, .text = "xy" },
            .{ .operation = .insert, .text = "34" },
            .{ .operation = .equal, .text = "z" },
            .{ .operation = .delete, .text = "cd" },
            .{ .operation = .insert, .text = "56" },
        };
        const d_after: []const Diff = &.{
            .{ .operation = .delete, .text = "abxyzcd" },
            .{ .operation = .insert, .text = "12xy34z56" },
        };
        try testing.checkAllAllocationFailures(
            testing.allocator,
            testDiffCleanupEfficiency,
            .{
                dmp,
                .{ .input = dslice, .expected = d_after },
            },
        );
    }
    { // High cost elimination.
        dmp.diff_edit_cost = 5;
        const dslice: []const Diff = &.{
            .{ .operation = .delete, .text = "ab" },
            .{ .operation = .insert, .text = "12" },
            .{ .operation = .equal, .text = "wxyz" },
            .{ .operation = .delete, .text = "cd" },
            .{ .operation = .insert, .text = "34" },
        };
        const d_after: []const Diff = &.{
            .{ .operation = .delete, .text = "abwxyzcd" },
            .{ .operation = .insert, .text = "12wxyz34" },
        };
        try testing.checkAllAllocationFailures(
            testing.allocator,
            testDiffCleanupEfficiency,
            .{
                dmp,
                .{ .input = dslice, .expected = d_after },
            },
        );
        dmp.diff_edit_cost = 4;
    }
}

test "diff before and after text" {
    const dmp = DiffMatchPatch{};
    const allocator = testing.allocator;
    const before = "The cat in the hat.";
    const after = "The bat in the belfry.";
    var diffs = try dmp.diff(allocator, before, after, false);
    defer deinitDiffList(allocator, &diffs);
    const before1 = try diffBeforeText(allocator, diffs);
    defer allocator.free(before1);
    const after1 = try diffAfterText(allocator, diffs);
    defer allocator.free(after1);
    try testing.expectEqualStrings(before, before1);
    try testing.expectEqualStrings(after, after1);
}

test diffIndex {
    const dmp = DiffMatchPatch{};
    {
        var diffs = try dmp.diff(
            testing.allocator,
            "The midnight train",
            "The blue midnight train",
            false,
        );
        defer deinitDiffList(testing.allocator, &diffs);
        try testing.expectEqual(0, diffIndex(diffs, 0));
        try testing.expectEqual(9, diffIndex(diffs, 4));
    }
    {
        var diffs = try dmp.diff(
            testing.allocator,
            "Better still to live and learn",
            "Better yet to learn and live",
            false,
        );
        defer deinitDiffList(testing.allocator, &diffs);
        try testing.expectEqual(11, diffIndex(diffs, 13));
        try testing.expectEqual(20, diffIndex(diffs, 21));
    }
}

test "diffPrettyFormat" {
    const test_deco = DiffDecorations{
        .delete_start = "<+>",
        .delete_end = "</+>",
        .insert_start = "<->",
        .insert_end = "</->",
        .equals_start = "<=>",
        .equals_end = "</=>",
    };
    const dmp = DiffMatchPatch{};
    const allocator = std.testing.allocator;
    var diffs = try dmp.diff(
        allocator,
        "A thing of beauty is a joy forever",
        "Singular beauty is enjoyed forever",
        false,
    );
    defer deinitDiffList(allocator, &diffs);
    try diffCleanupSemantic(allocator, &diffs);
    const out_text = try diffPrettyFormat(allocator, diffs, test_deco);
    defer allocator.free(out_text);
    try testing.expectEqualStrings(
        "<+>A thing of</+><->Singular</-><=> beauty is </=><+>a </+><->en</-><=>joy</=><->ed</-><=> forever</=>",
        out_text,
    );
}

fn testMapSubsetEquality(left: anytype, right: anytype) !void {
    var map_iter = left.iterator();
    while (map_iter.next()) |entry| {
        const key = entry.key_ptr.*;
        const value = entry.value_ptr.*;
        try testing.expectEqual(value, right.get(key));
    }
}
test "matchAlphabet" {
    var map = std.AutoHashMap(u8, usize).init(testing.allocator);
    defer map.deinit();
    try map.put('a', 4);
    try map.put('b', 2);
    try map.put('c', 1);
    var bitap_map = try matchAlphabet(testing.allocator, "abc");
    defer bitap_map.deinit();
    try testMapSubsetEquality(map, bitap_map);
    map.clearRetainingCapacity();
    try map.put('a', 37);
    try map.put('b', 18);
    try map.put('c', 8);
    var bitap_map2 = try matchAlphabet(testing.allocator, "abcaba");
    defer bitap_map2.deinit();
    try testMapSubsetEquality(map, bitap_map2);
}

fn testMatchBitap(
    allocator: Allocator,
    dmp: DiffMatchPatch,
    params: struct {
        text: []const u8,
        pattern: []const u8,
        loc: usize,
        expect: ?usize,
    },
) !void {
    const best_loc = try dmp.matchBitap(
        allocator,
        params.text,
        params.pattern,
        params.loc,
    );
    try testing.expectEqual(params.expect, best_loc);
}

test matchBitap {
    var dmp = DiffMatchPatch{};
    dmp.match_distance = 500;
    dmp.match_threshold = 0.5;
    // Exact match #1.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "fgh",
                .loc = 5,
                .expect = 5,
            },
        },
    );
    // Exact match #2.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "fgh",
                .loc = 0,
                .expect = 5,
            },
        },
    );
    // Fuzzy match #1
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "efxhi",
                .loc = 0,
                .expect = 4,
            },
        },
    );
    // Fuzzy match #2.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "cdefxyhijk",
                .loc = 5,
                .expect = 2,
            },
        },
    );
    // Fuzzy match #3.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "bxy",
                .loc = 1,
                .expect = null,
            },
        },
    );
    // Overflow.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "123456789xx0",
                .pattern = "3456789x0",
                .loc = 2,
                .expect = 2,
            },
        },
    );
    //Before start match.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdef",
                .pattern = "xxabc",
                .loc = 4,
                .expect = 0,
            },
        },
    );
    //
    // Beyond end match.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdef",
                .pattern = "defyy",
                .loc = 4,
                .expect = 3,
            },
        },
    );
    //  Oversized pattern.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdef",
                .pattern = "xabcdefy",
                .loc = 0,
                .expect = 0,
            },
        },
    );
    dmp.match_threshold = 0.4;
    // Threshold #1.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "efxyhi",
                .loc = 1,
                .expect = 4,
            },
        },
    );
    dmp.match_threshold = 0.3;
    //  Threshold #2.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "efxyhi",
                .loc = 1,
                .expect = null,
            },
        },
    );
    dmp.match_threshold = 0.0;
    //  Threshold #3.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijk",
                .pattern = "bcdef",
                .loc = 1,
                .expect = 1,
            },
        },
    );
    dmp.match_threshold = 0.5;
    //  Multiple select #1.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdexyzabcde",
                .pattern = "abccde",
                .loc = 5,
                .expect = 8,
            },
        },
    );
    dmp.match_distance = 10; // Strict location.
    //  Distance test #1.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijklmnopqrstuvwxyz",
                .pattern = "abcdefg",
                .loc = 1,
                .expect = 0,
            },
        },
    );
    // Distance test #2.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijklmnopqrstuvwxyz",
                .pattern = "abcdxxefg",
                .loc = 1,
                .expect = 0,
            },
        },
    );
    dmp.match_distance = 1000; // Loose location.
    //  Distance test #3.
    try testing.checkAllAllocationFailures(
        testing.allocator,
        testMatchBitap,
        .{
            dmp,
            .{
                .text = "abcdefghijklmnopqrstuvwxyz",
                .pattern = "abcdefg",
                .loc = 24,
                .expect = 0,
            },
        },
    );
}

test matchMain {
    var dmp = DiffMatchPatch{};
    dmp.match_threshold = 0.5;
    dmp.match_distance = 100;
    const allocator = testing.allocator;
    // Equality.
    try testing.expectEqual(0, dmp.matchMain(
        allocator,
        "abcdefg",
        "abcdefg",
        1000,
    ));
    // Null text
    try testing.expectEqual(null, dmp.matchMain(
        allocator,
        "",
        "abcdefg",
        1,
    ));
    // Null pattern.
    try testing.expectEqual(3, dmp.matchMain(
        allocator,
        "abcdefg",
        "",
        3,
    ));
    // Exact match.
    try testing.expectEqual(3, dmp.matchMain(
        allocator,
        "abcdefg",
        "de",
        3,
    ));
    // Beyond end match.
    try testing.expectEqual(3, dmp.matchMain(
        allocator,
        "abcdef",
        "defy",
        4,
    ));

    // Oversized pattern.
    try testing.expectEqual(0, dmp.matchMain(
        allocator,
        "abcdef",
        "abcdefy",
        0,
    ));
    dmp.match_threshold = 0.7;
    // Complex match.
    try testing.expectEqual(4, dmp.matchMain(
        allocator,
        "I am the very model of a modern major general.",
        " that berry ",
        5,
    ));
    dmp.match_threshold = 0.5;
}

fn testPatchToText(allocator: Allocator) !void {
    //
    var p: Patch = Patch{
        .start1 = 20,
        .start2 = 21,
        .length1 = 18,
        .length2 = 17,
        .diffs = try sliceToDiffList(allocator, &.{
            .{ .operation = .equal, .text = "jump" },
            .{ .operation = .delete, .text = "s" },
            .{ .operation = .insert, .text = "ed" },
            .{ .operation = .equal, .text = " over " },
            .{ .operation = .delete, .text = "the" },
            .{ .operation = .insert, .text = "a" },
            .{ .operation = .equal, .text = "\nlaz" },
        }),
    };
    defer p.deinit(allocator);
    const strp = "@@ -21,18 +22,17 @@\n jump\n-s\n+ed\n  over \n-the\n+a\n %0Alaz\n";
    const patch_str = try p.asText(allocator);
    defer allocator.free(patch_str);
    try testing.expectEqualStrings(strp, patch_str);
}

test "patch to text" {
    try std.testing.checkAllAllocationFailures(
        testing.allocator,
        testPatchToText,
        .{},
    );
}

fn testPatchRoundTrip(allocator: Allocator, patch_in: []const u8) !void {
    var patches = try patchFromText(allocator, patch_in);
    defer deinitPatchList(allocator, &patches);
    const patch_out = try patchToText(allocator, patches);
    defer allocator.free(patch_out);
    try testing.expectEqualStrings(patch_in, patch_out);
}

test "patch from text" {
    const allocator = testing.allocator;
    var p0 = try patchFromText(allocator, "");
    defer deinitPatchList(allocator, &p0);
    try testing.expectEqual(0, p0.items.len);
    try std.testing.checkAllAllocationFailures(
        testing.allocator,
        testPatchRoundTrip,
        .{"@@ -21,18 +22,17 @@\n jump\n-s\n+ed\n  over \n-the\n+a\n %0Alaz\n"},
    );
    try std.testing.checkAllAllocationFailures(
        allocator,
        testPatchRoundTrip,
        .{"@@ -1 +1 @@\n-a\n+b\n"},
    );
    try std.testing.checkAllAllocationFailures(
        testing.allocator,
        testPatchRoundTrip,
        .{"@@ -1,3 +0,0 @@\n-abc\n"},
    );
    try std.testing.checkAllAllocationFailures(
        testing.allocator,
        testPatchRoundTrip,
        .{"@@ -0,0 +1,3 @@\n+abc\n"},
    );
    try testing.expectError(error.BadPatchString, patchFromText(allocator, "Bad\nPatch\nString\n"));
}

fn testPatchAddContext(
    allocator: Allocator,
    dmp: DiffMatchPatch,
    patch_text: []const u8,
    text: []const u8,
    expect: []const u8,
) !void {
    _, var patch = try patchFromHeader(allocator, patch_text);
    defer patch.deinit(allocator);
    const patch_og = try patch.asText(allocator);
    defer allocator.free(patch_og);
    try testing.expectEqualStrings(patch_text, patch_og);
    try dmp.patchAddContext(allocator, &patch, text);
    const patch_out = try patch.asText(allocator);
    defer allocator.free(patch_out);
    try testing.expectEqualStrings(expect, patch_out);
}

test "testPatchAddContext" {
    const allocator = testing.allocator;
    var dmp = DiffMatchPatch{};
    dmp.patch_margin = 4;
    // Simple case.
    try std.testing.checkAllAllocationFailures(
        allocator,
        testPatchAddContext,
        .{
            dmp,
            "@@ -21,4 +21,10 @@\n-jump\n+somersault\n",
            "The quick brown fox jumps over the lazy dog.",
            "@@ -17,12 +17,18 @@\n fox \n-jump\n+somersault\n s ov\n",
        },
    );
    // Not enough trailing context.
    try std.testing.checkAllAllocationFailures(
        allocator,
        testPatchAddContext,
        .{
            dmp,
            "@@ -21,4 +21,10 @@\n-jump\n+somersault\n",
            "The quick brown fox jumps.",
            "@@ -17,10 +17,16 @@\n fox \n-jump\n+somersault\n s.\n",
        },
    );
    // Not enough leading context.
    try std.testing.checkAllAllocationFailures(
        allocator,
        testPatchAddContext,
        .{
            dmp,
            "@@ -3 +3,2 @@\n-e\n+at\n",
            "The quick brown fox jumps.",
            "@@ -1,7 +1,8 @@\n Th\n-e\n+at\n  qui\n",
        },
    );
    // Ambiguity.
    try std.testing.checkAllAllocationFailures(
        allocator,
        testPatchAddContext,
        .{
            dmp,
            "@@ -3 +3,2 @@\n-e\n+at\n",
            "The quick brown fox jumps.  The quick brown fox crashes.",
            "@@ -1,27 +1,28 @@\n Th\n-e\n+at\n  quick brown fox jumps. \n",
        },
    );
    // TODO: This will need some patches which check the Unicode handling.
}

fn testMakePatch(allocator: Allocator) !void {
    var dmp = DiffMatchPatch{};
    dmp.match_max_bits = 32; // Need this for compat with translated tests
    var null_patch = try dmp.makePatchFromTexts(allocator, "", "");
    defer deinitPatchList(allocator, &null_patch);
    const null_patch_text = try patchToText(allocator, null_patch);
    defer allocator.free(null_patch_text);
    try testing.expectEqualStrings("", null_patch_text);
    const text1 = "The quick brown fox jumps over the lazy dog.";
    const text2 = "That quick brown fox jumped over a lazy dog.";
    const expectedPatch = "@@ -1,8 +1,7 @@\n Th\n-at\n+e\n  qui\n@@ -21,17 +21,18 @@\n jump\n-ed\n+s\n  over \n-a\n+the\n  laz\n";
    // The second patch must be "-21,17 +21,18", not "-22,17 +21,18" due to rolling context.
    var diffs = try dmp.diff(allocator, text2, text1, true);
    deinitDiffList(allocator, &diffs);
    if (false) {
        var patches = try dmp.makePatchFromTexts(allocator, text2, text1);
        defer deinitPatchList(allocator, &patches);
        const patch_text = try patchToText(allocator, patches);
        defer allocator.free(patch_text);
        try testing.expectEqualStrings(expectedPatch, patch_text);
    }
}

test "testMakePatch" {
    if (true) {
        try testing.checkAllAllocationFailures(
            testing.allocator,
            testMakePatch,
            .{},
        );
    }
}
