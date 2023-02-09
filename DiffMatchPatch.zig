const DiffMatchPatch = @This();

const std = @import("std");
const ArrayListUnmanaged = std.ArrayListUnmanaged;

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
};

/// Number of microseconds to map a diff before giving up (0 for infinity).
diff_timeout: i64 = 1 * std.time.us_per_s,
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

/// It is recommended that you use an Arena for this operation.
pub fn diff(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    /// If false, then don't run a line-level diff first
    /// to identify the changed areas. If true, then run
    /// a faster slightly less optimal diff.
    check_lines: bool,
) DiffError!ArrayListUnmanaged(Diff) {
    const deadline = std.time.microTimestamp() + dmp.diff_timeout;
    return dmp.diffInternal(allocator, before, after, check_lines, deadline);
}

fn diffInternal(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    check_lines: bool,
    deadline: u64,
) DiffError!ArrayListUnmanaged(Diff) {
    // Check for equality (speedup).
    var diffs = ArrayListUnmanaged(Diff){};
    if (std.mem.eql(u8, before, after)) {
        if (before.len != 0) {
            diffs.append(allocator, Diff{ .operation = .equal, .text = before });
        }
        return diffs;
    }

    // Trim off common prefix (speedup).
    var common_length = diffCommonPrefix(before, after);
    const common_prefix = before[0..common_length];
    var trimmed_before = before[common_length..];
    var trimmed_after = after[common_length..];

    // Trim off common suffix (speedup).
    common_length = diffCommonSuffix(before, after);
    var common_suffix = before[before.len - common_length ..];
    trimmed_before = trimmed_before[0 .. before.len - common_length];
    trimmed_after = trimmed_after[0 .. after.len - common_length];

    // Compute the diff on the middle block.
    diffs = try dmp.diffCompute(allocator, before, after, checklines, deadline);

    // Restore the prefix and suffix.
    if (common_prefix.len != 0) {
        try diffs.insert(allocator, 0, Diff{ .operation = .equal, .text = common_prefix });
    }
    if (common_suffix.len != 0) {
        try diffs.append(allocator, Diff{ .operation = .equal, .text = common_suffix });
    }

    diffCleanupMerge(diffs);
    return diffs;
}

fn diffCommonPrefix(before: []const u8, after: []const u8) usize {
    const n = std.math.min(before.len, after.len);
    var i: usize = 0;

    while (i < n) : (i += 1) {
        if (before[i] != after[i]) {
            return i;
        }
    }

    return n;
}

fn diffCommonSuffix(before: []const u8, after: []const u8) usize {
    const n = std.math.min(before.len, after.len);
    var i: usize = 1;

    while (i <= n) : (i += 1) {
        if (before[before.len - i] != after[after.len - i]) {
            return i - 1;
        }
    }

    return n;
}

fn diffCompute(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    check_lines: bool,
    deadline: u64,
) DiffError!ArrayListUnmanaged(Diff) {
    var diffs = ArrayListUnmanaged(Diff){};

    if (before.len == 0) {
        // Just add some text (speedup).
        try diffs.append(allocator, Diff{ .operation = .insert, .text = after });
        return diffs;
    }

    if (after.len == 0) {
        // Just delete some text (speedup).
        try diffs.append(allocator, Diff{ .operation = .delete, .text = before });
        return diffs;
    }

    const long_text = if (before.len > after.len) before else after;
    const short_text = if (before.len > after.len) after else before;

    var short_text_in_long_text_index = std.mem.indexOf(u8, long_text, short_text);
    if (short_text_in_long_text_index) |index| {
        // Shorter text is inside the longer text (speedup).
        const op: Diff.Operation = if (before.len > after.len)
            .delete
        else
            .insert;
        try diffs.insert(allocator, Diff{ .operation = op, .text = long_text[0..index] });
        try diffs.insert(allocator, Diff{ .operation = .equal, .text = short_text });
        try diffs.insert(allocator, Diff{ .operation = op, .text = long_text[index + short_text.len ..] });
        return diffs;
    }

    if (short_text.len == 1) {
        // Single character string.
        // After the previous speedup, the character can't be an equality.
        try diffs.insert(allocator, Diff{ .operation = .delete, .text = before });
        try diffs.insert(allocator, Diff{ .operation = .insert, .text = after });
        return diffs;
    }

    // Check to see if the problem can be split in two.
    var maybe_half_match = dmp.diffHalfMatch(allocator, before, after);
    if (maybe_half_match) |half_match| {
        // A half-match was found, sort out the return data.

        // Send both pairs off for separate processing.
        var diffs_a = try dmp.diffInternal(allocator, half_match.prefix_before, half_match.prefix_after, check_lines, deadline);
        var diffs_b = try dmp.diffInternal(allocator, half_match.suffix_before, half_match.suffix_after, check_lines, deadline);
        defer diffs_b.deinit(allocator);

        // Merge the results.
        diffs = diffs_a;
        try diffs.append(allocator, Diff{ .operation = .equal, .text = half_match.common_middle });
        try diffs.appendSlice(allocator, diffs_b);
        return diffs;
    }

    if (check_lines and before.len > 100 and after.len > 100) {
        return diffLineMode(text1, text2, deadline);
    }

    return diffBisect(text1, text2, deadline);
}

const HalfMatchResult = ?struct {
    prefix_before: []const u8,
    suffix_before: []const u8,
    prefix_after: []const u8,
    suffix_after: []const u8,
    common_middle: []const u8,
};

fn diffHalfMatch(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
) DiffError!HalfMatchResult {
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
    var half_match_1 = try dmp.diffHalfMatchInternal(allocator, long_text, short_text, (long_text.len + 3) / 4);
    // Check again based on the third quarter.
    var half_match_2 = try dmp.diffHalfMatchInternal(allocator, long_text, short_text, (long_text.len + 1) / 2);

    var half_match: HalfMatchResult = undefined;
    if (half_match_1 == null and half_match_2 == null) {
        return null;
    } else if (half_match_2 == null) {
        half_match = half_match_1.?;
    } else if (half_match_1 == null) {
        half_match = half_match_2.?;
    } else {
        // Both matched. Select the longest.
        half_match = if (half_match_1.common_midle.len > half_match_2.common_midle.len) half_match_1 else half_match_2;
    }

    // A half-match was found, sort out the return data.
    if (before.len > after.len) {
        return half_match;
    } else {
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

fn diffHalfMatchInternal(
    _: DiffMatchPatch,
    allocator: std.mem.Allocator,
    long_text: []const u8,
    short_text: []const u8,
    i: usize,
) DiffError!HalfMatchResult {
    // Start with a 1/4 length Substring at position i as a seed.
    const seed = long_text[i .. long_text.len / 4];
    var j: isize = -1;

    var best_common = std.ArrayListUnmanaged(u8){};
    var best_long_text_a = "";
    var best_long_text_b = "";
    var best_short_text_a = "";
    var best_short_text_b = "";

    while (j < short_text.length and b: {
        j = (std.mem.indexOf(u8, short_text[j + 1 ..], seed, j + 1) orelse break :b false) + j + 1;
        break :b true;
    }) {
        var prefix_length = diffCommonPrefix(long_text[i..], short_text[j..]);
        var suffix_length = diffCommonSuffix(long_text[0..i], short_text[0..j]);
        if (best_common.items.len < suffix_length + prefix_length) {
            best_common.items.len = 0;
            try best_common.appendSlice(allocator, short_text[j - suffix_length .. suffix_length]);
            try best_common.appendSlice(allocator, short_text[j..prefix_length]);

            best_long_text_a = long_text[0 .. i - suffix_length];
            best_long_text_b = long_text[i + prefix_length ..];
            best_short_text_a = short_text[0 .. j - suffix_length];
            best_short_text_b = short_text[j + prefix_length ..];
        }
    }
    if (best_common.Length * 2 >= long_text.Length) {
        return .{
            .prefix_before = best_long_text_a,
            .suffix_before = best_long_text_b,
            .prefix_after = best_short_text_a,
            .suffix_after = best_short_text_b,
            .common_middle = best_common,
        };
    } else {
        return null;
    }
}
