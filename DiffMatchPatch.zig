const DiffMatchPatch = @This();

const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const Allocator = std.mem.Allocator;
const ArrayListUnmanaged = std.ArrayListUnmanaged;
const DiffList = ArrayListUnmanaged(Diff);

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
};

/// Number of milliseconds to map a diff before giving up (0 for infinity).
diff_timeout: u64 = 1000,
/// Cost of an empty edit operation in terms of edit characters.
diff_edit_cost: u16 = 4,

/// At what point is no match declared (0.0 = perfection, 1.0 = very loose).
match_threshold: f32 = 0.5,
/// How far to search for a match (0 = exact location, 1000+ = broad match).
/// A match this many characters away from the expected location will add
/// 1.0 to the score (0.0 is a perfect match).
match_distance: u32 = 1000,
/// The number of bits in an int.
match_max_bits: u16 = 32,

/// When deleting a large block of text (over ~64 characters), how close
/// do the contents have to be to match the expected contents. (0.0 =
/// perfection, 1.0 = very loose).  Note that Match_Threshold controls
/// how closely the end points of a delete need to match.
patch_delete_threshold: f32 = 0.5,
/// Chunk size for context length.
patch_margin: u16 = 4,

pub const DiffError = error{OutOfMemory};

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
    if (dmp.diff_timeout <= 0) {
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
            while (x1 < before_length and
                y1 < after_length and before[@intCast(x1)] == after[@intCast(y1)])
            {
                x1 += 1;
                y1 += 1;
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
            while (x2 < before_length and y2 < after_length and
                before[@intCast(before_length - x2 - 1)] ==
                after[@intCast(after_length - y2 - 1)])
            {
                x2 += 1;
                y2 += 1;
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
/// equalities.
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
    while (pointer < diffs.Count) {
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
    while (true) {
        const pattern = text1[text_length - length ..];
        const found = std.mem.indexOf(u8, text2, pattern) orelse
            return best;

        length += found;

        if (found == 0 or std.mem.eql(u8, text1[text_length - length ..], text2[0..length])) {
            best = length;
            length += 1;
        }
    }
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
        Diff.init(.insert, "δ"),
        Diff.init(.equal, "γ"),
    }), greek_diff.items);
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
