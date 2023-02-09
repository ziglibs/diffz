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
        return .{
            .prefix_before = half_match[2],
            .suffix_before = half_match[3],
            .prefix_after = half_match[0],
            .suffix_after = half_match[1],
            .common_middle = half_match[4],
        };
    }
}

fn diffHalfMatchInternal(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
) DiffError!HalfMatchResult {
    // Start with a 1/4 length Substring at position i as a seed.
    //   string seed = longtext.Substring(i, longtext.Length / 4);
    //   int j = -1;
    //   string best_common = string.Empty;
    //   string best_long_text_a = string.Empty, best_long_text_b = string.Empty;
    //   string best_short_text_a = string.Empty, best_short_text_b = string.Empty;
    //   while (j < shorttext.Length && (j = shorttext.IndexOf(seed, j + 1,
    //       StringComparison.Ordinal)) != -1) {
    //     int prefixLength = diffCommonPrefix(longtext.Substring(i),
    //                                          shorttext.Substring(j));
    //     int suffixLength = diffCommonSuffix(longtext.Substring(0, i),
    //                                          shorttext.Substring(0, j));
    //     if (best_common.Length < suffixLength + prefixLength) {
    //       best_common = shorttext.Substring(j - suffixLength, suffixLength)
    //           + shorttext.Substring(j, prefixLength);
    //       best_long_text_a = longtext.Substring(0, i - suffixLength);
    //       best_long_text_b = longtext.Substring(i + prefixLength);
    //       best_short_text_a = shorttext.Substring(0, j - suffixLength);
    //       best_short_text_b = shorttext.Substring(j + prefixLength);
    //     }
    //   }
    //   if (best_common.Length * 2 >= longtext.Length) {
    //     return new string[]{best_long_text_a, best_long_text_b,
    //         best_short_text_a, best_short_text_b, best_common};
    //   } else {
    //     return null;
    //   }
}

//
// Do a quick line-level diff on both strings, then rediff the parts for
// greater accuracy.
// This speedup can produce non-minimal diffs.
// @param text1 Old string to be diffed.
// @param text2 New string to be diffed.
// @param deadline Time when the diff should be complete by.
// @return List of Diff objects.
//
fn diff_lineMode(
    text1: []const u8,
    text2: []const u8,
    deadline: u64,
) DiffError!ArrayListUnmanaged(Diff) {
    // Scan the text on a line-by-line basis first.
    var a = diff_linesToChars(text1, text2);
    text1 = a[0];
    text2 = a[1];
    var linearray = a[2];

    var diffs: std.ArrayListUnmanaged(Diff) =
        diff_main(text1, text2, false, deadline);

    // Convert the diff back to original text.
    diff_charsToLines(diffs, linearray);
    // Eliminate freak matches (e.g. blank lines)
    diff_cleanupSemantic(diffs);

    // Rediff any replacement blocks, this time character-by-character.
    // Add a dummy entry at the end.
    try diffs.append(allocator, Diff(.equal, ""));
    var pointer: usize = 0;
    var count_delete: usize = 0;
    var count_insert: usize = 0;
    var text_delete: ArrayListUnmanaged(u8) = .{};
    var text_insert: ArrayListUnmanaged(u8) = .{};
    defer {
        text_delete.deinit(allocator);
        text_insert.deinit(allocator);
    }
    while (pointer < diffs.len) {
        switch (diffs[pointer].operation) {
            .insert => {
                count_insert += 1;
                // text_insert += diffs[pointer].text;
                text_insert.append(allocator, diffs[pointer].text);
            },
            .delete => {
                count_delete += 1;
                // text_delete += diffs[pointer].text;
                text_delete.append(allocator, diffs[pointer].text);
            },
            .equal => {
                // Upon reaching an equality, check for prior redundancies.
                if (count_delete >= 1 and count_insert >= 1) {
                    // Delete the offending records and add the merged ones.
                    // diffs.RemoveRange(pointer - count_delete - count_insert, count_delete + count_insert);
                    diffs.replaceRange(
                        allocator,
                        pointer - count_delete - count_insert,
                        count_delete + count_insert,
                        &.{},
                    );
                    pointer = pointer - count_delete - count_insert;
                    var subDiff = this.diff_main(text_delete, text_insert, false, deadline);
                    // diffs.InsertRange(pointer, subDiff);
                    try diffs.insertSlice(allocator, pointer, subDiff);
                    pointer = pointer + subDiff.items.len;
                }
                count_insert = 0;
                count_delete = 0;
                text_delete.items.len = 0;
                text_insert.items.len = 0;
            },
        }
        pointer += 1;
    }
    // diffs.RemoveAt(diffs.Count - 1); // Remove the dummy entry at the end.
    diffs.items.len -= 1;

    return diffs;
}
