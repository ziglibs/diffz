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
match_threshold: f32 = 0.05,

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
        if (d.text.len > 0) {
            allocator.free(d.text);
        }
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
    diffs: DiffList, // TODO This should be a Diff
    /// Start of patch in before text
    start1: usize = 0,
    length1: usize = 0,
    /// Start of patch in after text
    start2: usize = 0,
    length2: usize = 0,

    pub fn toString(patch: Patch) ![]const u8 {
        // TODO
        _ = patch;
    }

    pub fn writeTo(writer: anytype) !usize {
        // TODO
        _ = writer;
    }

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
        deinitDiffList(allocator, patch.diffs);
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
        try patch.writeText(writer, patch);
        return text_array.toOwnedSlice();
    }

    const format = std.fmt.format;

    /// Stream textual patch representation to Writer.  See `asText`
    /// for more information.
    pub fn writeText(writer: anytype, patch: Patch) !void {
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
        _ = writer.write(PATCH_TAIL);
        // Escape the body of the patch with %xx notation.
        for (patch.diffs) |a_diff| {
            switch (a_diff.operation) {
                .insert => try writer.writeByte('+'),
                .delete => try writer.writeByte('-'),
                .equal => try writer.writeByte('='),
            }
            _ = try writeUriEncoded(writer, diff.text);
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
        if (before.len != 0) {
            try diffs.append(allocator, Diff.init(.equal, try allocator.dupe(u8, before)));
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

    // Restore the prefix and suffix.
    if (common_prefix.len != 0) {
        try diffs.insert(allocator, 0, Diff.init(.equal, try allocator.dupe(u8, common_prefix)));
    }
    if (common_suffix.len != 0) {
        try diffs.append(allocator, Diff.init(.equal, try allocator.dupe(u8, common_suffix)));
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
    var diffs = DiffList{};

    if (before.len == 0) {
        // Just add some text (speedup).
        try diffs.append(allocator, Diff.init(.insert, try allocator.dupe(u8, after)));
        return diffs;
    }

    if (after.len == 0) {
        // Just delete some text (speedup).
        try diffs.append(allocator, Diff.init(.delete, try allocator.dupe(u8, before)));
        return diffs;
    }

    const long_text = if (before.len > after.len) before else after;
    const short_text = if (before.len > after.len) after else before;

    if (std.mem.indexOf(u8, long_text, short_text)) |index| {
        // Shorter text is inside the longer text (speedup).
        const op: Diff.Operation = if (before.len > after.len)
            .delete
        else
            .insert;
        // No need to adjust this index since any split is already valid
        try diffs.append(allocator, Diff.init(op, try allocator.dupe(u8, long_text[0..index])));
        try diffs.append(allocator, Diff.init(.equal, try allocator.dupe(u8, short_text)));
        try diffs.append(allocator, Diff.init(op, try allocator.dupe(u8, long_text[index + short_text.len ..])));
        return diffs;
    }

    if (short_text.len == 1) {
        // Single character string.
        // After the previous speedup, the character can't be an equality.
        try diffs.append(allocator, Diff.init(.delete, try allocator.dupe(u8, before)));
        try diffs.append(allocator, Diff.init(.insert, try allocator.dupe(u8, after)));
        return diffs;
    }

    // Check to see if the problem can be split in two.
    if (try dmp.diffHalfMatch(allocator, before, after)) |half_match| {
        // A half-match was found, sort out the return data.
        defer half_match.deinit(allocator);
        // Send both pairs off for separate processing.
        const diffs_a = try dmp.diffInternal(
            allocator,
            half_match.prefix_before,
            half_match.prefix_after,
            check_lines,
            deadline,
        );
        var diffs_b = try dmp.diffInternal(
            allocator,
            half_match.suffix_before,
            half_match.suffix_after,
            check_lines,
            deadline,
        );
        defer diffs_b.deinit(allocator);

        // Merge the results.
        diffs = diffs_a;
        try diffs.append(allocator, Diff.init(.equal, try allocator.dupe(u8, half_match.common_middle)));
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
    // Check again based on the third quarter.
    const half_match_2 = try dmp.diffHalfMatchInternal(allocator, long_text, short_text, (long_text.len + 1) / 2);

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
        return .{
            .prefix_before = try allocator.dupe(u8, best_long_text_a),
            .suffix_before = try allocator.dupe(u8, best_long_text_b),
            .prefix_after = try allocator.dupe(u8, best_short_text_a),
            .suffix_after = try allocator.dupe(u8, best_short_text_b),
            .common_middle = try best_common.toOwnedSlice(allocator),
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
    try diffs.append(allocator, Diff.init(.delete, try allocator.dupe(u8, before)));
    try diffs.append(allocator, Diff.init(.insert, try allocator.dupe(u8, after)));
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
    var diffsb = try dmp.diffInternal(allocator, text1b, text2b, false, deadline);
    defer diffsb.deinit(allocator);

    try diffs.appendSlice(allocator, diffsb.items);
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
        text.items.len = 0;
        var j: usize = 0;
        while (j < d.text.len) : (j += 1) {
            try text.appendSlice(allocator, line_array[d.text[j]]);
        }
        allocator.free(d.text);
        d.text = try allocator.dupe(u8, text.items);
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
                                defer allocator.free(ot);
                                @memcpy(nt[0..ot.len], ot);
                                @memcpy(nt[ot.len..], text_insert.items[0..common_length]);
                                diffs.items[ii].text = nt;
                            } else { // The prefix is at the start of the diffs
                                const text = try allocator.dupe(u8, text_insert.items[0..common_length]);
                                try diffs.insert(allocator, 0, Diff.init(.equal, text));
                                pointer += 1; // Keep pointer pointed at current diff
                            }
                            // Remove merged prefixes
                            try text_insert.replaceRange(allocator, 0, common_length, &.{});
                            try text_delete.replaceRange(allocator, 0, common_length, &.{});
                        }
                        // Factor out any common suffixies.
                        common_length = diffCommonSuffix(text_insert.items, text_delete.items);
                        if (common_length != 0) {
                            // Move the common part to the equal diff
                            const old_text = diffs.items[pointer].text;
                            defer allocator.free(old_text);
                            diffs.items[pointer].text = try std.mem.concat(allocator, u8, &.{
                                text_insert.items[text_insert.items.len - common_length ..],
                                old_text,
                            }); // Remove it from the ends of the insert/delete pair
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
                        try diffs.insert(allocator, pointer, Diff.init(
                            .delete,
                            try allocator.dupe(u8, text_delete.items),
                        ));
                        pointer += 1;
                    }
                    if (text_insert.items.len != 0) {
                        try diffs.insert(allocator, pointer, Diff.init(
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
                const pt = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer - 1].text,
                    diffs.items[pointer].text[0 .. diffs.items[pointer].text.len -
                        diffs.items[pointer - 1].text.len],
                });
                const p1t = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer - 1].text,
                    diffs.items[pointer + 1].text,
                });
                const old_pt = diffs.items[pointer].text;
                defer allocator.free(old_pt);
                const old_pt1t = diffs.items[pointer + 1].text;
                defer allocator.free(old_pt1t);
                diffs.items[pointer].text = pt;
                diffs.items[pointer + 1].text = p1t;
                freeRangeDiffList(allocator, diffs, pointer - 1, 1);
                try diffs.replaceRange(allocator, pointer - 1, 1, &.{});
                changes = true;
            } else if (std.mem.startsWith(u8, diffs.items[pointer].text, diffs.items[pointer + 1].text)) {
                const pm1t = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer - 1].text,
                    diffs.items[pointer + 1].text,
                });
                const pt = try std.mem.concat(allocator, u8, &.{
                    diffs.items[pointer].text[diffs.items[pointer + 1].text.len..],
                    diffs.items[pointer + 1].text,
                });
                const old_ptm1 = diffs.items[pointer - 1].text;
                defer allocator.free(old_ptm1);
                const old_pt = diffs.items[pointer].text;
                defer allocator.free(old_pt);
                diffs.items[pointer - 1].text = pm1t;
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
                try diffs.insert(
                    allocator,
                    @intCast(equalities.items[equalities.items.len - 1]),
                    Diff.init(.delete, try allocator.dupe(u8, last_equality.?)),
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
                    defer allocator.free(deletion);
                    defer allocator.free(insertion);
                    try diffs.insert(
                        allocator,
                        @intCast(pointer),
                        Diff.init(.equal, try allocator.dupe(u8, insertion[0..overlap_length1])),
                    );
                    diffs.items[@intCast(pointer - 1)].text =
                        try allocator.dupe(u8, deletion[0 .. deletion.len - overlap_length1]);
                    diffs.items[@intCast(pointer + 1)].text =
                        try allocator.dupe(u8, insertion[overlap_length1..]);
                    pointer += 1;
                }
            } else {
                if (@as(f32, @floatFromInt(overlap_length2)) >= @as(f32, @floatFromInt(deletion.len)) / 2.0 or
                    @as(f32, @floatFromInt(overlap_length2)) >= @as(f32, @floatFromInt(insertion.len)) / 2.0)
                {
                    // Reverse overlap found.
                    // Insert an equality and swap and trim the surrounding edits.
                    defer allocator.free(deletion);
                    defer allocator.free(insertion);
                    try diffs.insert(
                        allocator,
                        @intCast(pointer),
                        Diff.init(.equal, try allocator.dupe(u8, deletion[0..overlap_length2])),
                    );
                    diffs.items[@intCast(pointer - 1)].operation = .insert;
                    const new_minus = try allocator.dupe(u8, insertion[0 .. insertion.len - overlap_length2]);
                    diffs.items[@intCast(pointer - 1)].text = new_minus;
                    diffs.items[@intCast(pointer + 1)].operation = .delete;
                    const new_plus = try allocator.dupe(u8, deletion[overlap_length2..]);
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
                    allocator.free(diffs.items[pointer - 1].text);
                    diffs.items[pointer - 1].text = try allocator.dupe(u8, best_equality_1.items);
                } else {
                    const old_diff = diffs.orderedRemove(pointer - 1);
                    allocator.free(old_diff.text);
                    pointer -= 1;
                }
                allocator.free(diffs.items[pointer].text);
                diffs.items[pointer].text = try allocator.dupe(u8, best_edit.items);
                if (best_equality_2.items.len != 0) {
                    allocator.free(diffs.items[pointer + 1].text);
                    diffs.items[pointer + 1].text = try allocator.dupe(u8, best_equality_2.items);
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

/// Reduce the number of edits by eliminating operationally trivial
/// equalities.  TODO this needs tests
pub fn diffCleanupEfficiency(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    diffs: *DiffList,
) DiffError!void {
    var changes = false;
    // Stack of indices where equalities are found.
    var equalities = DiffList{};
    defer deinitDiffList(allocator, equalities);
    // Always equal to equalities[equalitiesLength-1][1]
    var last_equality = "";
    var pointer: isize = 0; // Index of current position.
    // Is there an insertion operation before the last equality.
    var pre_ins = false;
    // Is there a deletion operation before the last equality.
    var pre_del = false;
    // Is there an insertion operation after the last equality.
    var post_ins = false;
    // Is there a deletion operation after the last equality.
    var post_del = false;
    while (pointer < diffs.len) {
        if (diffs.items[pointer].operation == .equal) { // Equality found.
            if (diffs.items[pointer].text.len < dmp.diff_edit_cost and (post_ins or post_del)) {
                // Candidate found.
                equalities.Push(pointer);
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
            if ((last_equality.Length != 0) and
                ((pre_ins and pre_del and post_ins and post_del) or
                ((last_equality.Length < dmp.diff_edit_cost / 2) and
                ((if (pre_ins) 1 else 0) + (if (pre_del) 1 else 0) + (if (post_ins) 1 else 0) + (if (post_del) 1 else 0)) == 3)))
            {
                // Duplicate record.
                try diffs.insert(
                    allocator,
                    equalities.items[equalities.items.len - 1],
                    Diff.init(.delete, try allocator.dupe(u8, last_equality)),
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

                    pointer = if (equalities.items.len > 0) equalities.items[equalities.items.len - 1] else -1;
                    post_ins = false;
                    post_del = false;
                }
                changes = true;
            }
        }
        pointer += 1;
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
pub fn diffIndex(diffs: DiffList, loc: usize) usize {
    //      int chars1 = 0;
    //      int chars2 = 0;
    //      int last_chars1 = 0;
    //      int last_chars2 = 0;
    var chars1: usize = 0;
    var chars2: usize = 0;
    var last_chars1: usize = 0;
    var last_chars2: usize = 0;
    //  Dummy diff
    var last_diff: Diff = Diff{ .operation = .equal, .text = "" };
    for (diffs) |a_diff| {
        if (a_diff.operation != .insert) {
            // Equality or deletion.
            chars1 += a_diff.text.len;
        }
        if (a_diff.operation != .delete) {
            // Equality or insertion.
            chars2 += a_diff.text.len;
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
        return last_chars2;
    }
    // Add the remaining character length.
    return last_chars2 + (loc - last_chars1);
}

/// A struct holding bookends for `diffPrittyFormat(diffs)`.
///
/// May include a function taking an allocator and the diff,
/// which shall return the text of the diff, appropriately munged.
/// Note that if the function is provided, all text returned will
/// be freed, so it should always return a copy whether or not
/// edits are needed.
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
    var out = ArrayListUnmanaged(u8){};
    defer out.deinit(allocator);
    const writer = out.writer();
    _ = try writeDiffPrettyFormat(allocator, writer, diffs, deco);
    return out.toOwnedSlice(allocator);
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
    for (diffs) |d| {
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
            .equals => {
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
    for (diffs) |d| {
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
    for (diffs) |d| {
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
    for (diffs) |d| {
        switch (d.operation) {
            .insert => {
                inserts += d.text.len;
            },
            .delete => {
                deletes += d.text.len;
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

//| MATCH FUNCTIONS

/// Locate the best instance of 'pattern' in 'text' near 'loc'.
/// Returns -1 if no match found.
/// @param text The text to search.
/// @param pattern The pattern to search for.
/// @param loc The location to search around.
/// @return Best match index or -1.
pub fn matchMain(allocator: Allocator, text: []const u8, pattern: []const u8, passed_loc: usize) ?usize {
    // Clamp the loc to fit within text.
    const loc = @min(passed_loc, text.len);
    if (std.mem.eql(u8, text, pattern)) {
        // Shortcut (potentially not guaranteed by the algorithm)
        // TODO would be good to know what the above means...
        return 0;
    } else if (text.len == 0) {
        // Nothing to match.
        return null;
    } else if (loc + pattern.len <= text.len and std.mem.eql(u8, text[loc..pattern.length], pattern)) {
        // Perfect match at the perfect spot!  (Includes case of null pattern)
        return loc;
    } else {
        // Do a fuzzy compare.
        return matchBitap(allocator, text, pattern, loc);
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
    var threshold = @This().threshold;
    // Is there a nearby exact match? (speedup)
    var best_loc = std.mem.indexOfPos(u8, text, pattern);
    if (best_loc) |best| {
        threshold = @min(matchBitapScore(0, best, loc, pattern), threshold);
    }
    // What about in the other direction? (speedup)
    const trunc_text = text[0..@min(loc + pattern.len, text.len)];
    best_loc = std.mem.lastIndexOf(u8, trunc_text, pattern);
    if (best_loc) |best| {
        threshold = @min(matchBitapScore(0, best, loc, pattern), threshold);
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
            if (matchBitapScore(d, loc + bin_mid, loc, pattern) <= threshold) {
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
                const score = matchBitapScore(d, j - 1, loc, pattern);
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
        if (matchBitapScore(d + 1, loc, loc, pattern) > threshold) {
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

/// Locate the best instance of `pattern` in `text` near `loc` using the
/// Bitap algorithm.  Returns -1 if no match found.
///
/// @param text The text to search.
/// @param pattern The pattern to search for.
/// @param loc The location to search around.
/// @return Best match index or -1.
fn matchBitap(
    allocator: Allocator,
    text: []const u8,
    pattern: []const u8,
    loc: usize,
) !?usize {
    // TODO decide what to do here:
    // assert (Match_MaxBits == 0 || pattern.Length <= Match_MaxBits)
    //    : "Pattern too long for this application.";

    // Initialise the alphabet.
    var map = try matchAlphabet(allocator, pattern);
    defer map.deinit();
    // Highest score beyond which we give up.
    var threshold = @This().threshold;
    // Is there a nearby exact match? (speedup)
    var best_loc = std.mem.indexOfPos(u8, text, pattern);
    if (best_loc) |best| {
        threshold = @min(matchBitapScore(0, best, loc, pattern), threshold);
    }
    // TODO obviously if we want a speedup here, we do this:
    // if (threshold == 0.0) return best_loc;
    // We don't have to unwrap best_loc because the retval is ?usize already
    // What about in the other direction? (speedup)
    const trunc_text = text[0..@min(loc + pattern.len, text.len)];
    best_loc = std.mem.lastIndexOf(u8, trunc_text, pattern);
    if (best_loc) |best| {
        threshold = @min(matchBitapScore(0, best, loc, pattern), threshold);
    }
    // Initialise the bit arrays.
    const shift: u6 = @intCast(pattern.len - 1);
    const matchmask = 1 << shift;
    best_loc = null;
    var bin_min: usize = undefined;
    var bin_mid: usize = undefined;
    var bin_max = pattern.len + text.len;
    // null last_rd to simplying freeing memory
    var last_rd: []usize = try allocator.alloct(usize, 0);
    for (0..pattern.len) |d| {
        // Scan for the best match; each iteration allows for one more error.
        // Run a binary search to determine how far from 'loc' we can stray at
        // this error level.
        bin_min = 0;
        bin_mid = bin_max;
        while (bin_min < bin_mid) {
            if (matchBitapScore(d, loc + bin_mid, loc, pattern) <= threshold) {
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
        var rd: []usize = allocator.alloc(usize, finish + 2);
        const dshift: u6 = @intCast(d);
        rd[finish + 1] = (1 << dshift) - 1;
        var j = finish;
        while (j >= start) : (j -= 1) {
            const char_match: usize = if (text.len <= j - 1 or !map.contains(text[j - 1]))
                // Out of range.
                0
            else
                map.get(text[j - 1]);
            if (d == 0) {
                // First pass: exact match.
                rd[j] = ((rd[j + 1] << 1) | 1) & char_match;
            } else {
                // Subsequent passes: fuzzy match.
                rd[j] = ((rd[j + 1] << 1) | 1) & char_match | (((last_rd[j + 1] | last_rd[j]) << 1) | 1) | last_rd[j + 1];
            }
            if ((rd[j] & matchmask) != 0) {
                const score = matchBitapScore(d, j - 1, loc, pattern);
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
        if (matchBitapScore(d + 1, loc, loc, pattern) > threshold) {
            // No hope for a (better) match at greater error levels.
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
fn matchBitapScore(e: usize, x: usize, loc: usize, pattern: []const u8) f64 {
    // shortcut? TODO, proof in comments
    // if (e == 0 and x == loc) return 0.0;
    const e_float: f32 = @floatFromInt(e);
    const len_float: f32 = @floatFromInt(pattern.len);
    // if e == 0, accuracy == 0: 0/x = 0
    const accuracy = e_float / len_float;
    // if loc == x, proximity == 0
    const proximity = if (loc >= x) loc - x else x - loc;
    if (@This().match_distance == 0) {
        // Dodge divide by zero
        if (proximity == 0) // therefore this returns 0
            return accuracy
        else
            return 1.0;
    }
    const float_match: f64 = @floatFromInt(@This().match_distance);
    // or this is 0 + 0/f_m aka 0
    return accuracy + (proximity / float_match);
}

/// Initialise the alphabet for the Bitap algorithm.
/// @param pattern The text to encode.
/// @return Hash of character locations.
fn matchAlphabet(allocator: Allocator, pattern: []const u8) !std.HashMap(u8, usize) {
    var map = std.HashMap(u8, usize).init(allocator);
    errdefer map.deinit();
    for (pattern) |c| {
        if (!map.contains(c)) {
            try map.put(c, 0);
        }
    }
    for (pattern, 0..) |c, i| {
        const shift: u6 = @intCast(pattern.len - i - 1);
        const value: usize = map.get(c) | (1 << shift);
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

///
/// Increase the context until it is unique, but don't let the pattern
/// expand beyond DiffMatchPatch.match_max_bits.
///
/// @param patch The patch to grow.
/// @param text Source text.
fn patchAddContext(allocator: Allocator, patch: *Patch, text: []const u8) !void {
    if (text.len == 0) return;
    // TODO the fixup logic here might make patterns too large?
    // It should be ok, because big patches get broken up.  Hmm.
    var padding = 0;
    { // Grow the pattern around the patch until unique, to set padding amount.
        var pattern = text[patch.start2 .. patch.start2 + patch.length1];
        const max_width: usize = @This().match_max_bits - (2 * @This().patch_margin);
        while (std.mem.indexOf(u8, text, pattern) != std.mem.lastIndexOf(u8, text, pattern) and pattern.len < max_width) {
            padding += @This().patch_margin;
            const pat_start = @max(0, patch.start2 - padding);
            const pat_end = pat_start + @min(text.len, patch.start2 + patch.length1 + padding);
            pattern = text[pat_start..pat_end];
        }
    }
    // Add one chunk for good luck.
    padding += @This().patch_margin;
    // Add the prefix.
    const prefix = pre: {
        var pre_start = @max(0, patch.start2 - padding);
        // Make sure we're not breaking a codepoint.
        while (is_follow(text[pre_start]) and pre_start > 0) {
            pre_start -= 1;
        } // Assuming we did everything else right, pre_end should be
        // properly placed.
        const pre_end = pre_start + patch.start2;
        break :pre text[pre_start..pre_end];
    };
    if (prefix.len != 0) {
        try patch.diffs.append(
            allocator,
            Diff{
                .operation = .equal,
                .text = try allocator.dupe(u8, prefix),
            },
        );
    }
    // Add the suffix.
    const suffix = post: {
        const post_start = patch.start2 + patch.length1;
        // In case we messed up somewhere:
        assert(!is_follow(text[post_start]));
        var post_end = post_start + @min(text.len, patch.start2 + patch.length1 + padding);
        // Prevent broken codepoints here as well: Lead bytes, or follow with another follow
        while (!std.ascii.isASCII(text[post_end]) and post_end + 1 < text.len and is_follow(text[post_end + 1])) {
            post_end += 1;
            // Special case: penultimate with another follow at end
            if (post_end + 2 == text.len and is_follow(text[post_end + 1])) {
                post_end += 1;
                break; // Not actually necessary, but polite.
            }
        }
        break :post text[post_start..post_end];
    };
    if (suffix.len != 0) {
        try patch.diffs.append(
            allocator,
            Diff{
                .operation = .equal,
                .text = try allocator.dupe(u8, suffix),
            },
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

/// @return List of Patch objects.
fn makePatchInternal(
    allocator: Allocator,
    text: []const u8,
    diffs: DiffList,
    diff_act: DiffHandling,
) !PatchList {
    const patches = PatchList{};
    if (diffs.items.len == 0) {
        return patches; // Empty diff means empty patchlist
    }

    var patch = Patch{};
    var char_count1 = 0;
    var char_count2 = 0;
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
    for (diffs) |a_diff| {
        if (patch.diffs.items.len == 0 and a_diff.operation != .equal) {
            patch.start1 = char_count1;
            patch.start2 = char_count2;
        }
        switch (a_diff.operation) {
            .insert => {
                const d = if (diff_act == .copy) a_diff.clone(allocator) else a_diff;
                try patch.diffs.append(allocator, d);
                patch.length2 += a_diff.text.len;
                try postpatch.insertSlice(char_count2, a_diff.text);
            },
            .delete => {
                //
                const d = if (diff_act == .copy) a_diff.clone(allocator) else a_diff;
                try patch.diffs.append(allocator, d);
                patch.length1 += a_diff.text.len;
                try postpatch.replaceRange(char_count2, a_diff.text.len, .{});
            },
            .equal => {
                //
                if (a_diff.text.len <= 2 * @This().patch_margin and patch.diffs.items.len != 0 and a_diff != diffs.items[diffs.items.len]) {
                    // Small equality inside a patch.
                    const d = if (diff_act == .copy) a_diff.clone(allocator) else a_diff;
                    try patch.diffs.append(allocator, d);
                    patch.length1 += a_diff.text.len;
                    patch.length2 += a_diff.text.len;
                }
                if (a_diff.text.len >= 2 * @This().patch_margin) {
                    // Time for a new patch.
                    if (patch.diffs.items.len != 0) {
                        // free the Diff if we own it
                        if (diff_act == .own) {
                            allocator.free(a_diff.text);
                        }
                        try patchAddContext(allocator, patch, prepatch_text);
                        try patches.append(allocator, patch);
                        patch = Patch{};
                        // Unlike Unidiff, our patch lists have a rolling context.
                        // https://github.com/google/diff-match-patch/wiki/Unidiff
                        // Update prepatch text & pos to reflect the application of the
                        // just completed patch.
                        if (first_patch) {
                            // no free on first
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
        if (a_diff.operation != .remove) {
            char_count2 += a_diff.text.len;
        }
    } // end for loop

    // Pick up the leftover patch if not empty.
    if (patch.diffs.items.len != 0) {
        try patchAddContext(allocator, patch, prepatch_text);
        try patches.append(allocator, patch);
    }
}

/// Compute a list of patches to turn text1 into text2.
/// text2 is not provided, diffs are the delta between text1 and text2.
///
/// @param text1 Old text.
/// @param diffs Array of Diff objects for text1 to text2.
pub fn makePatch(allocator: Allocator, text: []const u8, diffs: DiffList) !PatchList {
    try makePatchInternal(allocator, text, diffs, .copy);
}

pub fn makePatchFromTexts(allocator: Allocator, text1: []const u8, text2: []const u8) !PatchList {
    const diffs = try diff(@This(), allocator, text1, text2, true);
    if (diffs.items.len > 2) {
        try diffCleanupSemantic(diffs);
        try diffCleanupEfficiency(diffs);
    }
    return try makePatchInternal(allocator, text1, diffs, .own);
}

pub fn makePatchFromDiffs(allocator: Allocator, diffs: DiffList) !PatchList {
    const text1 = try diffBeforeText(allocator, diffs);
    return try makePatch(allocator, text1, diffs, .copy);
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
pub fn patchApply(allocator: Allocator, og_patches: PatchList, og_text: []const u8) !struct { []const u8, bool } {
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
    const null_padding = try patchAddPadding(patches);
    var text_array = try std.ArrayList(u8).initCapacity(og_text.len);
    defer text_array.deinit();
    text_array.appendSlice(null_padding);
    text_array.appendSlice(og_text);
    text_array.appendSlice(null_padding);
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
        const m_max_b = @This().match_max_bits;
        if (text1.len > m_max_b) {
            // patchSplitMax will only provide an oversized pattern
            // in the case of a monster delete.
            maybe_start = matchMain(
                allocator,
                text_array.items,
                text1[0..m_max_b],
                expected_loc,
            );
            if (maybe_start) |start| {
                const e_start = text1.len - m_max_b;
                maybe_end = matchMain(
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
            maybe_start = matchMain(allocator, og_text, text1, expected_loc);
        }
        if (maybe_start) |start| {
            // Found a match.  :)
            delta = start - expected_loc;
            //          results[x] = true;
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
                const diffs = try diff(
                    @This(),
                    allocator,
                    text1,
                    text2,
                    false,
                );
                const t1_l_float: f64 = @floatFromInt(text1.len);
                const bad_match = diffLevenshtein(diffs) / t1_l_float > @This().patch_delete_threshold;
                if (text1.len > m_max_b and bad_match) {
                    // The end points match, but the content is unacceptably bad.
                    //              results[x] = false;
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
fn patchSplitMax(allocator: Allocator, patches: PatchList) !PatchList {
    const patch_size = @This().match_max_bits;
    const patch_margin = @This().patch_margin;
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
                try patch.diffs.append(
                    allocator,
                    Diff{
                        .operation = .equal,
                        .text = precontext.toOwnedSlice(),
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
                    const diff1 = bigpatch.diffs.orderedRemove(0);
                    try patch.diffs.append(allocator, diff1);
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
                        const diff1 = bigpatch.diffs.orderedRemove(0);
                        try patch.diffs.append(allocator, diff1);
                    } else {
                        // Free and dupe
                        const old_diff = bigpatch.diffs[0];
                        defer old_diff.deinit(allocator);
                        bigpatch.diffs[0] = Diff{
                            .operation = diff_type,
                            .text = try allocator.dupe(u8, new_diff_text),
                        };
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
            } else {
                // We didn't allocate memory, but it's polite to free it (?)
                allocator.free(postcontext);
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
fn patchAddPadding(allocator: Allocator, patches: PatchList) ![]const u8 {
    assert(patches.items.len != 0);
    const pad_len = @This().patch_margin;
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
        try diffs.insert(
            allocator,
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
        diffs.append(
            allocator,
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
    const text_array = try std.ArrayList(u8).init(allocator);
    defer text_array.deinit();
    const writer = text_array.writer();
    try writePatch(writer, patches);
    return text_array.toOwnedSlice();
}

/// Stream a `PatchList` to the provided Writer.
pub fn writePatch(writer: anytype, patches: PatchList) !void {
    for (patches) |a_patch| {
        try a_patch.writePatch(writer);
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
    var cursor = 0;
    while (cursor < text.len) {
        // TODO catch BadPatchString here and print diagnostic
        const cursor_delta, const patch = try patchFromHeader(allocator, text[cursor..]);
        cursor += cursor_delta;
        try patches.append(allocator, patch);
    }
}

fn countDigits(text: []const u8) usize {
    var idx = 0;
    while (std.ascii.isDigit(text[idx])) : (idx += 1) {}
    return idx;
}

fn patchFromHeader(allocator: Allocator, text: []const u8) !struct { usize, Patch } {
    var patch = Patch{};
    var cursor: usize = undefined;
    if (std.mem.eql(u8, text[0..4], PATCH_HEAD)) {
        // Parse location and length in before text
        patch.start1 = std.fmt.parseInt(
            usize,
            text[4..],
            10,
        ) catch return error.BadPatchString;
        cursor = 4 + countDigits(text[4..]);
        assert(cursor > 4);
        if (text[cursor] != ',') {
            cursor += 1;
            patch.start1 -= 1;
            patch.length1 = 1;
        } else {
            cursor += 1;
            patch.length1 = std.fmt.parseInt(
                usize,
                text[cursor..],
                10,
            ) catch return error.BadPatchString;
            const delta = countDigits(text[cursor..]);
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
        patch.start2 = std.fmt.parseInt(
            usize,
            text[cursor..],
            10,
        ) catch return error.BadPatchString;
        const delta1 = 4 + countDigits(text[4..]);
        assert(delta1 > 0);
        cursor += delta1;
        if (text[cursor] != ',') {
            cursor += 1;
            patch.start2 -= 1;
            patch.length2 = 1;
        } else {
            cursor += 1;
            patch.length2 = std.fmt.parseInt(
                usize,
                text[cursor..],
                10,
            ) catch return error.BadPatchString;
            const delta2 = countDigits(text[cursor..]);
            assert(delta2 > 1);
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
    const patch_lines = std.mem.splitScalar(
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
        const diff_line = try decodeUri(allocator, line) catch return error.BadPatchString;
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
    if (std.mem.indexOf(u8, line, '%')) |first| {
        // Text to decode.
        // Result will always be shorter than line:
        var new_line = try std.ArrayList(u8).initCapacity(allocator, line.len);
        defer new_line.init;
        try new_line.appendSlice(line[0..first]);
        var out_buf: [1]u8 = .{0};
        var codeunit = try std.fmt.hexToBytes(&out_buf, line[first + 1 .. first + 3]);
        try new_line.append(codeunit[0]);
        var cursor = first + 3;
        while (std.mem.indexOf(u8, line[cursor..], '%')) |next| {
            codeunit = try std.fmt.hexToBytes(&out_buf, line[next + 1 .. next + 3]);
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

// DONE [✅]: Allocate all text in diffs to
// not cause segfault while freeing; not a problem
// at the moment because we don't free anything :(
// (or was it??)

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

test diffHalfMatch {
    const allocator = testing.allocator;

    var one_timeout = DiffMatchPatch{};
    one_timeout.diff_timeout = 1;
    const dh1 = try one_timeout.diffHalfMatch(allocator, "1234567890", "abcdef");
    try testing.expectEqual(
        @as(?HalfMatchResult, null),
        dh1,
    ); // No match #1
    const dh2 = try one_timeout.diffHalfMatch(allocator, "12345", "23");
    try testing.expectEqual(
        @as(?HalfMatchResult, null),
        dh2,
    ); // No match #2

    // Single matches
    var dh3 = (try one_timeout.diffHalfMatch(allocator, "1234567890", "a345678z")).?;
    defer dh3.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "12",
        .suffix_before = "90",
        .prefix_after = "a",
        .suffix_after = "z",
        .common_middle = "345678",
    }, dh3); // Single Match #1

    var dh4 = (try one_timeout.diffHalfMatch(allocator, "a345678z", "1234567890")).?;
    defer dh4.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "a",
        .suffix_before = "z",
        .prefix_after = "12",
        .suffix_after = "90",
        .common_middle = "345678",
    }, dh4); // Single Match #2

    var dh5 = (try one_timeout.diffHalfMatch(allocator, "abc56789z", "1234567890")).?;
    defer dh5.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "abc",
        .suffix_before = "z",
        .prefix_after = "1234",
        .suffix_after = "0",
        .common_middle = "56789",
    }, dh5); // Single Match #3

    var dh6 = (try one_timeout.diffHalfMatch(allocator, "a23456xyz", "1234567890")).?;
    defer dh6.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "a",
        .suffix_before = "xyz",
        .prefix_after = "1",
        .suffix_after = "7890",
        .common_middle = "23456",
    }, dh6); // Single Match #4

    // Multiple matches
    var dh7 = (try one_timeout.diffHalfMatch(allocator, "121231234123451234123121", "a1234123451234z")).?;
    defer dh7.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "12123",
        .suffix_before = "123121",
        .prefix_after = "a",
        .suffix_after = "z",
        .common_middle = "1234123451234",
    }, dh7); // Multiple Matches #1

    var dh8 = (try one_timeout.diffHalfMatch(allocator, "x-=-=-=-=-=-=-=-=-=-=-=-=", "xx-=-=-=-=-=-=-=")).?;
    defer dh8.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "",
        .suffix_before = "-=-=-=-=-=",
        .prefix_after = "x",
        .suffix_after = "",
        .common_middle = "x-=-=-=-=-=-=-=",
    }, dh8); // Multiple Matches #2

    var dh9 = (try one_timeout.diffHalfMatch(allocator, "-=-=-=-=-=-=-=-=-=-=-=-=y", "-=-=-=-=-=-=-=yy")).?;
    defer dh9.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "-=-=-=-=-=",
        .suffix_before = "",
        .prefix_after = "",
        .suffix_after = "y",
        .common_middle = "-=-=-=-=-=-=-=y",
    }, dh9); // Multiple Matches #3

    // Other cases
    // Optimal diff would be -q+x=H-i+e=lloHe+Hu=llo-Hew+y not -qHillo+x=HelloHe-w+Hulloy
    var dh10 = (try one_timeout.diffHalfMatch(allocator, "qHilloHelloHew", "xHelloHeHulloy")).?;
    defer dh10.deinit(allocator);
    try testing.expectEqualDeep(HalfMatchResult{
        .prefix_before = "qHillo",
        .suffix_before = "w",
        .prefix_after = "x",
        .suffix_after = "Hulloy",
        .common_middle = "HelloHe",
    }, dh10); // Non-optimal halfmatch

    one_timeout.diff_timeout = 0;
    try testing.expectEqualDeep(@as(?HalfMatchResult, null), try one_timeout.diffHalfMatch(allocator, "qHilloHelloHew", "xHelloHeHulloy")); // Non-optimal halfmatch
}

test diffLinesToChars {
    const allocator = std.testing.allocator;
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

test diffCharsToLines {
    const allocator = std.testing.allocator;
    const equal_a = Diff.init(.equal, try allocator.dupe(u8, "a"));
    defer allocator.free(equal_a.text);
    const insert_a = Diff.init(.insert, try allocator.dupe(u8, "a"));
    defer allocator.free(insert_a.text);
    const equal_b = Diff.init(.equal, try allocator.dupe(u8, "b"));
    defer allocator.free(equal_b.text);
    const delete_b = Diff.init(.delete, try allocator.dupe(u8, "b"));
    defer allocator.free(delete_b.text);
    try testing.expect(equal_a.eql(equal_a));
    try testing.expect(!insert_a.eql(equal_a));
    try testing.expect(!equal_a.eql(equal_b));
    try testing.expect(!equal_a.eql(delete_b));

    // Convert chars up to lines.
    var diffs = DiffList{};
    defer deinitDiffList(allocator, &diffs);
    try diffs.appendSlice(allocator, &.{
        Diff{ .operation = .equal, .text = try allocator.dupe(u8, "\u{0001}\u{0002}\u{0001}") },
        Diff{ .operation = .insert, .text = try allocator.dupe(u8, "\u{0002}\u{0001}\u{0002}") },
    });
    var tmp_vector = std.ArrayList([]const u8).init(allocator);
    defer tmp_vector.deinit();
    try tmp_vector.append("");
    try tmp_vector.append("alpha\n");
    try tmp_vector.append("beta\n");
    try diffCharsToLines(allocator, diffs.items, tmp_vector.items);

    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        Diff.init(.equal, "alpha\nbeta\nalpha\n"),
        Diff.init(.insert, "beta\nalpha\nbeta\n"),
    }), diffs.items);

    // TODO: Implement exhaustive tests
}

test diffCleanupMerge {
    const allocator = std.testing.allocator;
    // Cleanup a messy diff.
    var diffs = DiffList{};
    defer deinitDiffList(allocator, &diffs);

    try testing.expectEqualDeep(
        @as([]const Diff, &[0]Diff{}),
        diffs.items,
    ); // Null case

    try diffs.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "b"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "c"),
        },
    });
    try diffCleanupMerge(allocator, &diffs);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ .{ .operation = .equal, .text = "a" }, .{ .operation = .delete, .text = "b" }, .{ .operation = .insert, .text = "c" } }), diffs.items); // No change case

    var diffs2 = DiffList{};
    defer deinitDiffList(allocator, &diffs2);
    try diffs2.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "b"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "c"),
        },
    });
    try diffCleanupMerge(allocator, &diffs2);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .equal, .text = "abc" },
    }), diffs2.items); // Merge equalities

    var diffs3 = DiffList{};
    defer deinitDiffList(allocator, &diffs3);
    try diffs3.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "b"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "c"),
        },
    });
    try diffCleanupMerge(allocator, &diffs3);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .delete, .text = "abc" },
    }), diffs3.items); // Merge deletions

    var diffs4 = DiffList{};
    defer deinitDiffList(allocator, &diffs4);
    try diffs4.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "b"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "c"),
        },
    });
    try diffCleanupMerge(allocator, &diffs4);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .insert, .text = "abc" },
    }), diffs4.items); // Merge insertions

    var diffs5 = DiffList{};
    defer deinitDiffList(allocator, &diffs5);
    try diffs5.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "b"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "c"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "d"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "e"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "f"),
        },
    });
    try diffCleanupMerge(allocator, &diffs5);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .delete, .text = "ac" },
        .{ .operation = .insert, .text = "bd" },
        .{ .operation = .equal, .text = "ef" },
    }), diffs5.items); // Merge interweave

    var diffs6 = DiffList{};
    defer deinitDiffList(allocator, &diffs6);
    try diffs6.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "abc"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "dc"),
        },
    });
    try diffCleanupMerge(allocator, &diffs6);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .equal, .text = "a" },
        .{ .operation = .delete, .text = "d" },
        .{ .operation = .insert, .text = "b" },
        .{ .operation = .equal, .text = "c" },
    }), diffs6.items); // Prefix and suffix detection

    var diffs7 = DiffList{};
    defer deinitDiffList(allocator, &diffs7);
    try diffs7.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "x"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "abc"),
        },
        .{
            .operation = .delete,
            .text = try allocator.dupe(u8, "dc"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "y"),
        },
    });
    try diffCleanupMerge(allocator, &diffs7);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .equal, .text = "xa" },
        .{ .operation = .delete, .text = "d" },
        .{ .operation = .insert, .text = "b" },
        .{ .operation = .equal, .text = "cy" },
    }), diffs7.items); // Prefix and suffix detection with equalities

    var diffs8 = DiffList{};
    defer deinitDiffList(allocator, &diffs8);
    try diffs8.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "a"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "ba"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "c"),
        },
    });
    try diffCleanupMerge(allocator, &diffs8);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .insert, .text = "ab" },
        .{ .operation = .equal, .text = "ac" },
    }), diffs8.items); // Slide edit left

    var diffs9 = DiffList{};
    defer deinitDiffList(allocator, &diffs9);
    try diffs9.appendSlice(allocator, &[_]Diff{
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "c"),
        },
        .{
            .operation = .insert,
            .text = try allocator.dupe(u8, "ab"),
        },
        .{
            .operation = .equal,
            .text = try allocator.dupe(u8, "a"),
        },
    });
    try diffCleanupMerge(allocator, &diffs9);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        .{ .operation = .equal, .text = "ca" },
        .{ .operation = .insert, .text = "ba" },
    }), diffs9.items); // Slide edit right

    var diffs10 = DiffList{};
    defer deinitDiffList(allocator, &diffs10);
    try diffs10.appendSlice(allocator, &[_]Diff{
        Diff.init(
            .equal,
            try allocator.dupe(u8, "a"),
        ),
        Diff.init(
            .delete,
            try allocator.dupe(u8, "b"),
        ),
        Diff.init(
            .equal,
            try allocator.dupe(u8, "c"),
        ),
        Diff.init(
            .delete,
            try allocator.dupe(u8, "ac"),
        ),
        Diff.init(
            .equal,
            try allocator.dupe(u8, "x"),
        ),
    });
    try diffCleanupMerge(allocator, &diffs10);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        Diff.init(.delete, "abc"),
        Diff.init(.equal, "acx"),
    }), diffs10.items); // Slide edit left recursive

    var diffs11 = DiffList{};
    defer deinitDiffList(allocator, &diffs11);
    try diffs11.appendSlice(allocator, &[_]Diff{
        Diff.init(
            .equal,
            try allocator.dupe(u8, "x"),
        ),
        Diff.init(
            .delete,
            try allocator.dupe(u8, "ca"),
        ),
        Diff.init(
            .equal,
            try allocator.dupe(u8, "c"),
        ),
        Diff.init(
            .delete,
            try allocator.dupe(u8, "b"),
        ),
        Diff.init(
            .equal,
            try allocator.dupe(u8, "a"),
        ),
    });
    try diffCleanupMerge(allocator, &diffs11);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        Diff.init(.equal, "xca"),
        Diff.init(.delete, "cba"),
    }), diffs11.items); // Slide edit right recursive

    var diffs12 = DiffList{};
    defer deinitDiffList(allocator, &diffs12);
    try diffs12.appendSlice(allocator, &[_]Diff{
        Diff.init(
            .delete,
            try allocator.dupe(u8, "b"),
        ),
        Diff.init(
            .insert,
            try allocator.dupe(u8, "ab"),
        ),
        Diff.init(
            .equal,
            try allocator.dupe(u8, "c"),
        ),
    });
    try diffCleanupMerge(allocator, &diffs12);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        Diff.init(.insert, "a"),
        Diff.init(.equal, "bc"),
    }), diffs12.items); // Empty merge

    var diffs13 = DiffList{};
    defer deinitDiffList(allocator, &diffs13);
    try diffs13.appendSlice(allocator, &[_]Diff{
        Diff.init(.equal, ""),
        Diff.init(.insert, try allocator.dupe(u8, "a")),
        Diff.init(.equal, try allocator.dupe(u8, "b")),
    });
    try diffCleanupMerge(allocator, &diffs13);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{
        Diff.init(.insert, "a"),
        Diff.init(.equal, "b"),
    }), diffs13.items); // Empty equality
}

test diffCleanupSemanticLossless {
    const allocator = std.testing.allocator;
    var diffs = DiffList{};
    try diffCleanupSemanticLossless(allocator, &diffs);
    try testing.expectEqualDeep(@as([]const Diff, &[0]Diff{}), diffs.items); // Null case

    var diffs2 = DiffList{};
    defer deinitDiffList(allocator, &diffs2);
    try diffs2.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "AAA\r\n\r\nBBB")),
        Diff.init(.insert, try allocator.dupe(u8, "\r\nDDD\r\n\r\nBBB")),
        Diff.init(.equal, try allocator.dupe(u8, "\r\nEEE")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs2);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "AAA\r\n\r\n"),
        Diff.init(.insert, "BBB\r\nDDD\r\n\r\n"),
        Diff.init(.equal, "BBB\r\nEEE"),
    }), diffs2.items);

    var diffs3 = DiffList{};
    defer deinitDiffList(allocator, &diffs3);
    try diffs3.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "AAA\r\nBBB")),
        Diff.init(.insert, try allocator.dupe(u8, " DDD\r\nBBB")),
        Diff.init(.equal, try allocator.dupe(u8, " EEE")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs3);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "AAA\r\n"),
        Diff.init(.insert, "BBB DDD\r\n"),
        Diff.init(.equal, "BBB EEE"),
    }), diffs3.items);

    var diffs4 = DiffList{};
    defer deinitDiffList(allocator, &diffs4);
    try diffs4.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "The c")),
        Diff.init(.insert, try allocator.dupe(u8, "ow and the c")),
        Diff.init(.equal, try allocator.dupe(u8, "at.")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs4);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "The "),
        Diff.init(.insert, "cow and the "),
        Diff.init(.equal, "cat."),
    }), diffs4.items);

    var diffs5 = DiffList{};
    defer deinitDiffList(allocator, &diffs5);
    try diffs5.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "The-c")),
        Diff.init(.insert, try allocator.dupe(u8, "ow-and-the-c")),
        Diff.init(.equal, try allocator.dupe(u8, "at.")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs5);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "The-"),
        Diff.init(.insert, "cow-and-the-"),
        Diff.init(.equal, "cat."),
    }), diffs5.items);

    var diffs6 = DiffList{};
    defer deinitDiffList(allocator, &diffs6);
    try diffs6.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "a")),
        Diff.init(.delete, try allocator.dupe(u8, "a")),
        Diff.init(.equal, try allocator.dupe(u8, "ax")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs6);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.delete, "a"),
        Diff.init(.equal, "aax"),
    }), diffs6.items);

    var diffs7 = DiffList{};
    defer deinitDiffList(allocator, &diffs7);
    try diffs7.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "xa")),
        Diff.init(.delete, try allocator.dupe(u8, "a")),
        Diff.init(.equal, try allocator.dupe(u8, "a")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs7);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "xaa"),
        Diff.init(.delete, "a"),
    }), diffs7.items);

    var diffs8 = DiffList{};
    defer deinitDiffList(allocator, &diffs8);
    try diffs8.appendSlice(allocator, &.{
        Diff.init(.equal, try allocator.dupe(u8, "The xxx. The ")),
        Diff.init(.insert, try allocator.dupe(u8, "zzz. The ")),
        Diff.init(.equal, try allocator.dupe(u8, "yyy.")),
    });
    try diffCleanupSemanticLossless(allocator, &diffs8);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "The xxx."),
        Diff.init(.insert, " The zzz."),
        Diff.init(.equal, " The yyy."),
    }), diffs8.items);
}

fn rebuildtexts(allocator: std.mem.Allocator, diffs: DiffList) ![2][]const u8 {
    var text = [2]std.ArrayList(u8){
        std.ArrayList(u8).init(allocator),
        std.ArrayList(u8).init(allocator),
    };

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

test diffBisect {
    const allocator = std.testing.allocator;
    // Normal.
    const a = "cat";
    const b = "map";
    // Since the resulting diff hasn't been normalized, it would be ok if
    // the insertion and deletion pairs are swapped.
    // If the order changes, tweak this test as required.
    var diffs = DiffList{};
    defer deinitDiffList(allocator, &diffs);
    var this = default;
    try diffs.appendSlice(allocator, &.{
        Diff.init(.delete, try allocator.dupe(u8, "c")),
        Diff.init(.insert, try allocator.dupe(u8, "m")),
        Diff.init(.equal, try allocator.dupe(u8, "a")),
        Diff.init(.delete, try allocator.dupe(u8, "t")),
        Diff.init(.insert, try allocator.dupe(u8, "p")),
    });
    // Travis TODO not sure if maxInt(u64) is correct for  DateTime.MaxValue
    var diff_bisect = try this.diffBisect(
        allocator,
        a,
        b,
        std.math.maxInt(u64),
    );
    defer deinitDiffList(allocator, &diff_bisect);
    try testing.expectEqualDeep(diffs, diff_bisect); // Normal.

    // Timeout.
    var diffs2 = DiffList{};
    defer deinitDiffList(allocator, &diffs2);
    try diffs2.appendSlice(allocator, &.{
        Diff.init(.delete, try allocator.dupe(u8, "cat")),
        Diff.init(.insert, try allocator.dupe(u8, "map")),
    });
    // Travis TODO not sure if 0 is correct for  DateTime.MinValue
    var diff_bisect2 = try this.diffBisect(allocator, a, b, 0);
    defer deinitDiffList(allocator, &diff_bisect2);
    try testing.expectEqualDeep(diffs2, diff_bisect2); // Timeout.
}

const talloc = testing.allocator;

test diff {
    var arena = std.heap.ArenaAllocator.init(talloc);
    defer arena.deinit();
    const allocator = std.testing.allocator;

    // Perform a trivial diff.
    var diffs = DiffList{};
    defer diffs.deinit(arena.allocator());
    var this = DiffMatchPatch{};
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "", "", false)).items); // diff: Null case.

    // TODO This is the last set of tests using the arena.  Someone should
    // rewrite them not to do so. -Sam
    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{Diff.init(.equal, "abc")});
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "abc", "abc", false)).items); // diff: Equality.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.equal, "ab"), Diff.init(.insert, "123"), Diff.init(.equal, "c") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "abc", "ab123c", false)).items); // diff: Simple insertion.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.equal, "a"), Diff.init(.delete, "123"), Diff.init(.equal, "bc") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "a123bc", "abc", false)).items); // diff: Simple deletion.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.equal, "a"), Diff.init(.insert, "123"), Diff.init(.equal, "b"), Diff.init(.insert, "456"), Diff.init(.equal, "c") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "abc", "a123b456c", false)).items); // diff: Two insertions.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.equal, "a"), Diff.init(.delete, "123"), Diff.init(.equal, "b"), Diff.init(.delete, "456"), Diff.init(.equal, "c") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "a123b456c", "abc", false)).items); // diff: Two deletions.

    // Perform a real diff.
    // Switch off the timeout.
    this.diff_timeout = 0;
    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.delete, "a"), Diff.init(.insert, "b") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "a", "b", false)).items); // diff: Simple case #1.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.delete, "Apple"), Diff.init(.insert, "Banana"), Diff.init(.equal, "s are a"), Diff.init(.insert, "lso"), Diff.init(.equal, " fruit.") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "Apples are a fruit.", "Bananas are also fruit.", false)).items); // diff: Simple case #2.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.delete, "a"), Diff.init(.insert, "\u{0680}"), Diff.init(.equal, "x"), Diff.init(.delete, "\t"), Diff.init(.insert, "\x00") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "ax\t", "\u{0680}x\x00", false)).items); // diff: Simple case #3.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.delete, "1"), Diff.init(.equal, "a"), Diff.init(.delete, "y"), Diff.init(.equal, "b"), Diff.init(.delete, "2"), Diff.init(.insert, "xab") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "1ayb2", "abxab", false)).items); // diff: Overlap #1.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.insert, "xaxcx"), Diff.init(.equal, "abc"), Diff.init(.delete, "y") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "abcy", "xaxcxabc", false)).items); // diff: Overlap #2.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.delete, "ABCD"), Diff.init(.equal, "a"), Diff.init(.delete, "="), Diff.init(.insert, "-"), Diff.init(.equal, "bcd"), Diff.init(.delete, "="), Diff.init(.insert, "-"), Diff.init(.equal, "efghijklmnopqrs"), Diff.init(.delete, "EFGHIJKLMNOefg") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "ABCDa=bcd=efghijklmnopqrsEFGHIJKLMNOefg", "a-bcd-efghijklmnopqrs", false)).items); // diff: Overlap #3.

    diffs.items.len = 0;
    try diffs.appendSlice(arena.allocator(), &.{ Diff.init(.insert, " "), Diff.init(.equal, "a"), Diff.init(.insert, "nd"), Diff.init(.equal, " [[Pennsylvania]]"), Diff.init(.delete, " and [[New") });
    try testing.expectEqualDeep(diffs.items, (try this.diff(arena.allocator(), "a [[Pennsylvania]] and [[New", " and [[Pennsylvania]]", false)).items); // diff: Large equality.

    // end of Arena Zone

    this.diff_timeout = 100; // 100ms
    // Increase the text lengths by 1024 times to ensure a timeout.
    {
        const a = "`Twas brillig, and the slithy toves\nDid gyre and gimble in the wabe:\nAll mimsy were the borogoves,\nAnd the mome raths outgrabe.\n" ** 1024;
        const b = "I am the very model of a modern major general,\nI've information vegetable, animal, and mineral,\nI know the kings of England, and I quote the fights historical,\nFrom Marathon to Waterloo, in order categorical.\n" ** 1024;
        const start_time = std.time.milliTimestamp();
        var time_diff = try this.diff(allocator, a, b, false);
        defer deinitDiffList(allocator, &time_diff);
        const end_time = std.time.milliTimestamp();
        // Test that we took at least the timeout period.
        try testing.expect(this.diff_timeout <= end_time - start_time); // diff: Timeout min.
        // Test that we didn't take forever (be forgiving).
        // Theoretically this test could fail very occasionally if the
        // OS task swaps or locks up for a second at the wrong moment.
        try testing.expect((this.diff_timeout) * 10000 * 2 > end_time - start_time); // diff: Timeout max.
        this.diff_timeout = 0;
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
        try testing.expectEqualDeep(diff_checked, diff_unchecked); // diff: Simple line-mode.
    }
    {
        const a = "1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890";
        const b = "abcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghijabcdefghij";
        var diff_checked = try this.diff(allocator, a, b, true);
        defer deinitDiffList(allocator, &diff_checked);
        var diff_unchecked = try this.diff(allocator, a, b, false);
        defer deinitDiffList(allocator, &diff_unchecked);
        try testing.expectEqualDeep(diff_checked, diff_unchecked); // diff: Single line-mode.
    }

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
    try testing.expectEqualDeep(texts_textmode, texts_linemode); // diff: Overlap line-mode.
}

test "Unicode diffs" {
    const allocator = std.testing.allocator;
    const this = DiffMatchPatch{};
    var greek_diff = try this.diff(
        allocator,
        "αβγ",
        "αβδ",
        false,
    );
    defer deinitDiffList(allocator, &greek_diff);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "αβ"),
        Diff.init(.delete, "γ"),
        Diff.init(.insert, "δ"),
    }), greek_diff.items);
    // ө is 0xd3, 0xa9, թ is 0xd6, 0xa9
    var prefix_diff = try this.diff(
        allocator,
        "abө",
        "abթ",
        false,
    );
    defer deinitDiffList(allocator, &prefix_diff);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "ab"),
        Diff.init(.delete, "ө"),
        Diff.init(.insert, "թ"),
    }), prefix_diff.items);
    var mid_diff = try this.diff(
        allocator,
        "αөβ",
        "αթβ",
        false,
    );
    defer deinitDiffList(allocator, &mid_diff);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "α"),
        Diff.init(.delete, "ө"),
        Diff.init(.insert, "թ"),
        Diff.init(.equal, "β"),
    }), mid_diff.items);

    var mid_prefix = try this.diff(
        allocator,
        "αβλ",
        "αδλ",
        false,
    );
    defer deinitDiffList(allocator, &mid_prefix);
    try testing.expectEqualDeep(@as([]const Diff, &.{
        Diff.init(.equal, "α"),
        Diff.init(.delete, "β"),
        Diff.init(.insert, "δ"),
        Diff.init(.equal, "λ"),
    }), mid_prefix.items);
}

test diffCleanupSemantic {
    const alloc = std.testing.allocator;
    // Cleanup semantically trivial equalities.
    // Null case.
    var diffs_empty = DiffList{};
    defer deinitDiffList(alloc, &diffs_empty);
    // var this = default;
    try diffCleanupSemantic(alloc, &diffs_empty);
    try testing.expectEqual(@as(usize, 0), diffs_empty.items.len); // Null case

    var diffs = DiffList{};
    defer deinitDiffList(alloc, &diffs);
    diffs.items.len = 0;
    try diffs.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "ab")),
        Diff.init(.insert, try alloc.dupe(u8, "cd")),
        Diff.init(.equal, try alloc.dupe(u8, "12")),
        Diff.init(.delete, try alloc.dupe(u8, "e")),
    });
    try diffCleanupSemantic(alloc, &diffs);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // No elimination #1
        Diff.init(.delete, "ab"),
        Diff.init(.insert, "cd"),
        Diff.init(.equal, "12"),
        Diff.init(.delete, "e"),
    }), diffs.items);

    var diffs2 = DiffList{};
    defer deinitDiffList(alloc, &diffs2);
    diffs2.items.len = 0;
    try diffs2.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "abc")),
        Diff.init(.insert, try alloc.dupe(u8, "ABC")),
        Diff.init(.equal, try alloc.dupe(u8, "1234")),
        Diff.init(.delete, try alloc.dupe(u8, "wxyz")),
    });
    try diffCleanupSemantic(alloc, &diffs2);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // No elimination #2
        Diff.init(.delete, "abc"),
        Diff.init(.insert, "ABC"),
        Diff.init(.equal, "1234"),
        Diff.init(.delete, "wxyz"),
    }), diffs2.items);

    var diffs3 = DiffList{};
    defer deinitDiffList(alloc, &diffs3);
    try diffs3.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "a")),
        Diff.init(.equal, try alloc.dupe(u8, "b")),
        Diff.init(.delete, try alloc.dupe(u8, "c")),
    });
    try diffCleanupSemantic(alloc, &diffs3);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Simple elimination
        Diff.init(.delete, "abc"),
        Diff.init(.insert, "b"),
    }), diffs3.items);

    var diffs4 = DiffList{};
    defer deinitDiffList(alloc, &diffs4);
    try diffs4.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "ab")),
        Diff.init(.equal, try alloc.dupe(u8, "cd")),
        Diff.init(.delete, try alloc.dupe(u8, "e")),
        Diff.init(.equal, try alloc.dupe(u8, "f")),
        Diff.init(.insert, try alloc.dupe(u8, "g")),
    });
    try diffCleanupSemantic(alloc, &diffs4);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Backpass elimination
        Diff.init(.delete, "abcdef"),
        Diff.init(.insert, "cdfg"),
    }), diffs4.items);

    var diffs5 = DiffList{};
    defer deinitDiffList(alloc, &diffs5);
    try diffs5.appendSlice(alloc, &.{
        Diff.init(.insert, try alloc.dupe(u8, "1")),
        Diff.init(.equal, try alloc.dupe(u8, "A")),
        Diff.init(.delete, try alloc.dupe(u8, "B")),
        Diff.init(.insert, try alloc.dupe(u8, "2")),
        Diff.init(.equal, try alloc.dupe(u8, "_")),
        Diff.init(.insert, try alloc.dupe(u8, "1")),
        Diff.init(.equal, try alloc.dupe(u8, "A")),
        Diff.init(.delete, try alloc.dupe(u8, "B")),
        Diff.init(.insert, try alloc.dupe(u8, "2")),
    });
    try diffCleanupSemantic(alloc, &diffs5);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Multiple elimination
        Diff.init(.delete, "AB_AB"),
        Diff.init(.insert, "1A2_1A2"),
    }), diffs5.items);

    var diffs6 = DiffList{};
    defer deinitDiffList(alloc, &diffs6);
    try diffs6.appendSlice(alloc, &.{
        Diff.init(.equal, try alloc.dupe(u8, "The c")),
        Diff.init(.delete, try alloc.dupe(u8, "ow and the c")),
        Diff.init(.equal, try alloc.dupe(u8, "at.")),
    });
    try diffCleanupSemantic(alloc, &diffs6);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Word boundaries
        Diff.init(.equal, "The "),
        Diff.init(.delete, "cow and the "),
        Diff.init(.equal, "cat."),
    }), diffs6.items);

    var diffs7 = DiffList{};
    defer deinitDiffList(alloc, &diffs7);
    try diffs7.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "abcxx")),
        Diff.init(.insert, try alloc.dupe(u8, "xxdef")),
    });
    try diffCleanupSemantic(alloc, &diffs7);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // No overlap elimination
        Diff.init(.delete, "abcxx"),
        Diff.init(.insert, "xxdef"),
    }), diffs7.items);

    var diffs8 = DiffList{};
    defer deinitDiffList(alloc, &diffs8);
    try diffs8.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "abcxxx")),
        Diff.init(.insert, try alloc.dupe(u8, "xxxdef")),
    });
    try diffCleanupSemantic(alloc, &diffs8);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Overlap elimination
        Diff.init(.delete, "abc"),
        Diff.init(.equal, "xxx"),
        Diff.init(.insert, "def"),
    }), diffs8.items);

    var diffs9 = DiffList{};
    defer deinitDiffList(alloc, &diffs9);
    try diffs9.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "xxxabc")),
        Diff.init(.insert, try alloc.dupe(u8, "defxxx")),
    });
    try diffCleanupSemantic(alloc, &diffs9);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Reverse overlap elimination
        Diff.init(.insert, "def"),
        Diff.init(.equal, "xxx"),
        Diff.init(.delete, "abc"),
    }), diffs9.items);

    var diffs10 = DiffList{};
    defer deinitDiffList(alloc, &diffs10);
    try diffs10.appendSlice(alloc, &.{
        Diff.init(.delete, try alloc.dupe(u8, "abcd1212")),
        Diff.init(.insert, try alloc.dupe(u8, "1212efghi")),
        Diff.init(.equal, try alloc.dupe(u8, "----")),
        Diff.init(.delete, try alloc.dupe(u8, "A3")),
        Diff.init(.insert, try alloc.dupe(u8, "3BC")),
    });
    try diffCleanupSemantic(alloc, &diffs10);
    try testing.expectEqualDeep(@as([]const Diff, &[_]Diff{ // Two overlap eliminations
        Diff.init(.delete, "abcd"),
        Diff.init(.equal, "1212"),
        Diff.init(.insert, "efghi"),
        Diff.init(.equal, "----"),
        Diff.init(.delete, "A"),
        Diff.init(.equal, "3"),
        Diff.init(.insert, "BC"),
    }), diffs10.items);
}
