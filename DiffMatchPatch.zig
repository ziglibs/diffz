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
    const deadline = @intCast(u64, std.time.microTimestamp()) + dmp.diff_timeout;
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
    diffs = try dmp.diffCompute(allocator, before, after, check_lines, deadline);

    // Restore the prefix and suffix.
    if (common_prefix.len != 0) {
        try diffs.insert(allocator, 0, Diff{ .operation = .equal, .text = common_prefix });
    }
    if (common_suffix.len != 0) {
        try diffs.append(allocator, Diff{ .operation = .equal, .text = common_suffix });
    }

    dmp.diffCleanupMerge(allocator, diffs);
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
        return dmp.diffLineMode(allocator, before, after, deadline);
    }

    return dmp.diffBisect(allocator, before, after, deadline);
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

    while (j < short_text.len and b: {
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
    if (best_common.len * 2 >= long_text.len) {
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

fn diffBisect(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    before: []const u8,
    after: []const u8,
    deadline: u64,
) ArrayListUnmanaged(Diff) {
    const max_d = (before.len + after.len + 1) / 2;
    const v_offset = max_d;
    const v_length = 2 * max_d;

    var v1 = try ArrayListUnmanaged(isize).initCapacity(v_length);
    v1.items.len = v_length;
    var v2 = try ArrayListUnmanaged(isize).initCapacity(v_length);
    v2.items.len = v_length;

    var x: usize = 0;
    while (x < v_length) : (x += 1) {
        v1.items[x] = -1;
        v2.items[x] = -1;
    }
    v1.items[v_offset + 1] = 0;
    v2.items[v_offset + 1] = 0;
    var delta = before.len - after.len;
    // If the total number of characters is odd, then the front path will
    // collide with the reverse path.
    var front = (delta % 2 != 0);
    // Offsets for start and end of k loop.
    // Prevents mapping of space beyond the grid.
    var k1start: usize = 0;
    var k1end: usize = 0;
    var k2start: usize = 0;
    var k2end: usize = 0;

    var d: usize = 0;
    while (d < max_d) : (d += 1) {
        // Bail out if deadline is reached.
        if (@intCast(u64, std.time.microTimestamp()) > deadline) {
            break;
        }

        // Walk the front path one step.
        var k1: isize = -d + k1start;
        while (k1 <= d - k1end) : (k1 += 2) {
            var k1_offset: isize = v_offset + k1;
            var x1: isize = 0;
            if (k1 == -d or k1 != d and v1[k1_offset - 1] < v1[k1_offset + 1]) {
                x1 = v1[k1_offset + 1];
            } else {
                x1 = v1[k1_offset - 1] + 1;
            }
            var y1: isize = x1 - k1;
            while (x1 < before.len and y1 < after.len and before[x1] == after[y1]) {
                x1 += 1;
                y1 += 1;
            }
            v1[k1_offset] = x1;
            if (x1 > before.len) {
                // Ran off the right of the graph.
                k1end += 2;
            } else if (y1 > after.len) {
                // Ran off the bottom of the graph.
                k1start += 2;
            } else if (front) {
                var k2_offset: isize = v_offset + delta - k1;
                if (k2_offset >= 0 and k2_offset < v_length and v2[k2_offset] != -1) {
                    // Mirror x2 onto top-left coordinate system.
                    var x2: isize = before.len - v2[k2_offset];
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
            var k2_offset: isize = v_offset + k2;
            var x2: isize = 0;
            if (k2 == -d or k2 != d and v2[k2_offset - 1] < v2[k2_offset + 1]) {
                x2 = v2[k2_offset + 1];
            } else {
                x2 = v2[k2_offset - 1] + 1;
            }
            var y2: isize = x2 - k2;
            while (x2 < before.len and y2 < after.len and before[before.len - x2 - 1] == after[after.len - y2 - 1]) {
                x2 += 1;
                y2 += 1;
            }
            v2[k2_offset] = x2;
            if (x2 > before.len) {
                // Ran off the left of the graph.
                k2end += 2;
            } else if (y2 > after.len) {
                // Ran off the top of the graph.
                k2start += 2;
            } else if (!front) {
                var k1_offset: isize = v_offset + delta - k2;
                if (k1_offset >= 0 and k1_offset < v_length and v1[k1_offset] != -1) {
                    var x1: isize = v1[k1_offset];
                    var y1: isize = v_offset + x1 - k1_offset;
                    // Mirror x2 onto top-left coordinate system.
                    x2 = before.len - v2[k2_offset];
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
    var diffs = ArrayListUnmanaged(Diff);
    try diffs.append(allocator, Diff{ .operation = .delete, .text = before });
    try diffs.append(allocator, Diff{ .operation = .insert, .text = after });
    return diffs;
}

fn diffBisectSplit(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    text1: []const u8,
    text2: []const u8,
    x: isize,
    y: isize,
    deadline: u64,
) DiffError!ArrayListUnmanaged(Diff) {
    const text1a = text1[0..x];
    const text2a = text2[0..y];
    const text1b = text1[x..];
    const text2b = text2[y..];

    // Compute both diffs serially.
    var diffs = try dmp.diffInternal(allocator, text1a, text2a, false, deadline);
    var diffsb = try dmp.diffInternal(allocator, text1b, text2b, false, deadline);
    defer diffs.deinit(allocator);

    try diffs.appendSlice(allocator, diffsb);
    return diffs;
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
fn diffLineMode(
    dmp: DiffMatchPatch,
    allocator: std.mem.Allocator,
    text1: []const u8,
    text2: []const u8,
    deadline: u64,
) DiffError!ArrayListUnmanaged(Diff) {
    // Scan the text on a line-by-line basis first.
    var a = try diffLinesToChars(allocator, text1, text2);
    text1 = a.chars_1;
    text2 = a.chars_2;
    var line_array = a.line_array;

    var diffs: ArrayListUnmanaged(Diff) = try dmp.diffInternal(allocator, text1, text2, false, deadline);

    // Convert the diff back to original text.
    try diffCharsToLines(allocator, diffs.items, line_array.items);
    // Eliminate freak matches (e.g. blank lines)
    try diffCleanupSemantic(allocator, diffs);

    // Rediff any replacement blocks, this time character-by-character.
    // Add a dummy entry at the end.
    try diffs.append(allocator, Diff{ .operation = .equal, .text = "" });

    var pointer: usize = 0;
    var count_delete: usize = 0;
    var count_insert: usize = 0;
    var text_delete = ArrayListUnmanaged(u8){};
    var text_insert = ArrayListUnmanaged(u8){};
    defer {
        text_delete.deinit(allocator);
        text_insert.deinit(allocator);
    }

    while (pointer < diffs.len) {
        switch (diffs.items[pointer].operation) {
            .insert => {
                count_insert += 1;
                // text_insert += diffs.items[pointer].text;
                try text_insert.append(allocator, diffs.items[pointer].text);
            },
            .delete => {
                count_delete += 1;
                // text_delete += diffs.items[pointer].text;
                try text_delete.append(allocator, diffs.items[pointer].text);
            },
            .equal => {
                // Upon reaching an equality, check for prior redundancies.
                if (count_delete >= 1 and count_insert >= 1) {
                    // Delete the offending records and add the merged ones.
                    // diffs.RemoveRange(pointer - count_delete - count_insert, count_delete + count_insert);
                    try diffs.replaceRange(
                        allocator,
                        pointer - count_delete - count_insert,
                        count_delete + count_insert,
                        &.{},
                    );
                    pointer = pointer - count_delete - count_insert;
                    var sub_diff = dmp.diffInternal(allocator, text_delete, text_insert, false, deadline);
                    // diffs.InsertRange(pointer, sub_diff);
                    try diffs.insertSlice(allocator, pointer, sub_diff);
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
    // diffs.RemoveAt(diffs.Count - 1); // Remove the dummy entry at the end.
    diffs.items.len -= 1;

    return diffs;
}

const LinesToCharsResult = struct {
    chars_1: []const u8,
    chars_2: []const u8,
    line_array: ArrayListUnmanaged([]const u8),
};

fn diffLinesToChars(allocator: std.mem.Allocator, text1: []const u8, text2: []const u8) error{OutOfMemory}!LinesToCharsResult {
    var line_array = ArrayListUnmanaged([]const u8){};
    var line_hash = std.StringHashMapUnmanaged(usize){};
    // e.g. line_array[4] == "Hello\n"
    // e.g. line_hash.get("Hello\n") == 4

    // "\x00" is a valid character, but various debuggers don't like it.
    // So we'll insert a junk entry to avoid generating a null character.
    try line_array.append(allocator, "");

    // Allocate 2/3rds of the space for text1, the rest for text2.
    var chars1 = diffLinesToCharsMunge(allocator, text1, &line_array, &line_hash, 40000);
    var chars2 = diffLinesToCharsMunge(allocator, text2, &line_array, &line_hash, 65535);
    return .{ .chars_1 = chars1, .chars_2 = chars2, .line_array = line_array };
}

fn diffLinesToCharsMunge(
    allocator: std.mem.Allocator,
    text: []const u8,
    line_array: *ArrayListUnmanaged([]const u8),
    line_hash: *std.StringHashMapUnmanaged(usize),
    max_lines: usize,
) error{OutOfMemory}![]const u8 {
    var line_start: isize = 0;
    var line_end: isize = -1;
    var line = "";
    var chars = ArrayListUnmanaged(u8){};
    // Walk the text, pulling out a Substring for each line.
    // text.split('\n') would would temporarily double our memory footprint.
    // Modifying text would create many large strings to garbage collect.
    while (line_end < text.len - 1) {
        line_end = b: {
            break :b (std.mem.indexOf(u8, text[line_start..], "\n") orelse break :b text.len - 1) + line_start;
        };
        line = text[line_start .. line_end + 1 - line_start];

        if (line_hash.get(line)) |value| {
            try chars.append(allocator, value);
        } else {
            if (line_array.items.len == max_lines) {
                // Bail out at 65535 because char 65536 == char 0.
                line = text[line_start..];
                line_end = text.len;
            }
            try line_array.append(allocator, line);
            try line_hash.put(allocator, line, line_array.items.len - 1);
            try chars.append(allocator, line_array.items.len - 1);
        }
        line_start = line_end + 1;
    }
    return try chars.items.toOwnedSlice(allocator);
}

fn diffCharsToLines(allocator: std.mem.Allocator, diffs: []Diff, line_array: []const []const u8) void {
    var text = ArrayListUnmanaged(u8){};
    for (diffs) |d| {
        text.items.len = 0;
        var j: usize = 0;
        while (j < diff.text.Length) : (j += 1) {
            try text.append(allocator, line_array[d.text[j]]);
        }
        d.text = text;
    }
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

//
// Reorder and merge like edit sections.  Merge equalities.
// Any edit section can move as long as it doesn't cross an equality.
// @param diffs List of Diff objects.
//
fn diffCleanupMerge(diffs: std.ArrayListUnmanaged(Diff), allocator: mem.Allocator) !void {
    // Add a dummy entry at the end.
    try diffs.append(allocator, Diff{ .operation = .equal, .text = "" });
    var pointer: usize = 0;
    var count_delete: usize = 0;
    var count_insert: usize = 0;
    var text_delete = std.ArrayList(u8).init(allocator);
    var text_insert = std.ArrayList(u8).init(allocator);
    var commonlength: usize = undefined;
    while (pointer < diffs.items.len) {
        switch (diffs[pointer].operation) {
            .insert => {
                count_insert += 1;
                text_insert += diffs[pointer].text;
                pointer += 1;
            },
            .delete => {
                count_delete += 1;
                text_delete += diffs[pointer].text;
                pointer += 1;
            },
            .equal => {
                // Upon reaching an equality, check for prior redundancies.
                if (count_delete + count_insert > 1) {
                    if (count_delete != 0 and count_insert != 0) {
                        // Factor out any common prefixies.
                        commonlength = this.diffCommonPrefix(text_insert, text_delete);
                        if (commonlength != 0) {
                            if ((pointer - count_delete - count_insert) > 0 and
                                diffs[pointer - count_delete - count_insert - 1].operation == .equal)
                            {
                                // diffs[pointer - count_delete - count_insert - 1].text
                                //     += text_insert.Substring(0, commonlength);
                                try diffs[pointer - count_delete - count_insert - 1].text.append(allocator, text_insert.items[0..commonlength]);
                            } else {
                                // diffs.Insert(0, new Diff(Operation.EQUAL,
                                //    text_insert.Substring(0, commonlength)));
                                const text = std.ArrayListUnmanaged(u8){ .items = try allocator.dupe(u8, text_insert.items[0..commonlength]) };
                                diffs.insert(0, Diff{ .operation = .equal, .text = text });
                                pointer += 1;
                            }
                            text_insert.len = commonlength;
                            text_delete.len = commonlength;
                        }
                        // Factor out any common suffixies.
                        // @ZigPort this seems very wrong
                        commonlength = this.diffCommonSuffix(text_insert, text_delete);
                        if (commonlength != 0) {
                            diffs[pointer].text = try std.mem.concat(allocator, &.{ text_insert.items[
                                text_insert.items.len - commonlength
                            ], diffs[pointer].text });
                            text_insert.items.len -= commonlength;
                            text_delete.items.len -= commonlength;
                        }
                    }
                    // Delete the offending records and add the merged ones.
                    pointer -= count_delete + count_insert;
                    try diffs.replaceRange(allocator, pointer, count_delete + count_insert, &.{});

                    if (text_delete.items.len != 0) {
                        try diffs.replaceRange(allocator, pointer, 0, &.{Diff{ .operation = .delete, .text = text_delete }});
                        pointer += 1;
                    }
                    if (text_insert.Length != 0) {
                        try diffs.replaceRange(allocator, pointer, 0, &.{Diff{ .operation = .insert, .text = text_insert }});
                        pointer += 1;
                    }
                    pointer += 1;
                } else if (pointer != 0 and diffs[pointer - 1].operation == .equal) {
                    // Merge this equality with the previous one.
                    try diffs[pointer - 1].text.append(allocator, diffs[pointer].text.items);
                    diffs.orderedRemove(pointer);
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
    if (diffs[diffs.items.len - 1].text.items.len == 0) {
        diffs.items.len -= 1;
    }

    // Second pass: look for single edits surrounded on both sides by
    // equalities which can be shifted sideways to eliminate an equality.
    // e.g: A<ins>BA</ins>C -> <ins>AB</ins>AC
    var changes = false;
    pointer = 1;
    // Intentionally ignore the first and last element (don't need checking).
    while (pointer < (diffs.items.len - 1)) {
        if (diffs[pointer - 1].operation == .equal and
            diffs[pointer + 1].operation == .equal)
        {
            // This is a single edit surrounded by equalities.
            if (mem.endsWith(u8, diffs[pointer].text.items, diffs[pointer - 1].text.items)) {
                // Shift the edit over the previous equality.
                diffs[pointer].text = diffs[pointer - 1].text +
                    diffs[pointer].text.Substring(0, diffs[pointer].text.Length -
                    diffs[pointer - 1].text.Length);
                diffs[pointer + 1].text = diffs[pointer - 1].text + diffs[pointer + 1].text;
                try diffs.replaceRange(allocator, pointer - 1, 1, &.{});
                changes = true;
            } else if (mem.startsWith(u8, diffs[pointer].text.items, diffs[pointer + 1].text.items)) {
                // Shift the edit over the next equality.
                diffs[pointer - 1].text += diffs[pointer + 1].text;
                diffs[pointer].text =
                    diffs[pointer].text.Substring(diffs[pointer + 1].text.Length) + diffs[pointer + 1].text;
                try diffs.replaceRange(allocator, pointer + 1, 1, &.{});
                changes = true;
            }
        }
        pointer += 1;
    }
    // If shifts were made, the diff needs reordering and another shift sweep.
    if (changes) {
        this.diff_cleanupMerge(diffs);
    }
}

fn diffCleanupSemantic(allocator: std.mem.Allocator, diffs: ArrayListUnmanaged(Diff)) error{OutOfMemory}!void {
    var changes = false;
    // Stack of indices where equalities are found.
    var equalities = ArrayListUnmanaged(isize){};
    // Always equal to equalities[equalitiesLength-1][1]
    var lastEquality: ?[]const u8 = null;
    var pointer: isize = 0; // Index of current position.
    // Number of characters that changed prior to the equality.
    var length_insertions1: usize = 0;
    var length_deletions1: usize = 0;
    // Number of characters that changed after the equality.
    var length_insertions2: usize = 0;
    var length_deletions2: usize = 0;
    while (pointer < diffs.items.len) {
        if (diffs.items[pointer].operation == .equal) { // Equality found.
            equalities.Push(pointer);
            length_insertions1 = length_insertions2;
            length_deletions1 = length_deletions2;
            length_insertions2 = 0;
            length_deletions2 = 0;
            lastEquality = diffs.items[pointer].text;
        } else { // an insertion or deletion
            if (diffs.items[pointer].operation == .equal) {
                length_insertions2 += diffs.items[pointer].text.Length;
            } else {
                length_deletions2 += diffs.items[pointer].text.Length;
            }
            // Eliminate an equality that is smaller or equal to the edits on both
            // sides of it.
            if (lastEquality != null and (lastEquality.Length <= std.math.max(length_insertions1, length_deletions1)) and (lastEquality.length <= std.math.max(length_insertions2, length_deletions2))) {
                // Duplicate record.
                diffs.Insert(equalities.Peek(), Diff{ .operation = .delete, .text = lastEquality });
                // Change second copy to insert.
                diffs.items[equalities.Peek() + 1].operation = .insert;
                // Throw away the equality we just deleted.
                equalities.Pop();
                if (equalities.Count > 0) {
                    equalities.Pop();
                }
                pointer = if (equalities.items.len > 0) equalities.Peek() else -1;
                length_insertions1 = 0; // Reset the counters.
                length_deletions1 = 0;
                length_insertions2 = 0;
                length_deletions2 = 0;
                lastEquality = null;
                changes = true;
            }
        }
        pointer += 1;
    }

    // Normalize the diff.
    if (changes) {
        diffCleanupMerge(diffs);
    }
    diffCleanupSemanticLossless(diffs);

    // Find any overlaps between deletions and insertions.
    // e.g: <del>abcxxx</del><ins>xxxdef</ins>
    //   -> <del>abc</del>xxx<ins>def</ins>
    // e.g: <del>xxxabc</del><ins>defxxx</ins>
    //   -> <ins>def</ins>xxx<del>abc</del>
    // Only extract an overlap if it is as big as the edit ahead or behind it.
    pointer = 1;
    while (pointer < diffs.Count) {
        if (diffs.items[pointer - 1].operation == Operation.DELETE and
            diffs.items[pointer].operation == Operation.INSERT)
        {
            const deletion = diffs.items[pointer - 1].text.items;
            const insertion = diffs.items[pointer].text.items;
            var overlap_length1: isize = diff_commonOverlap(deletion, insertion);
            var overlap_length2: isize = diff_commonOverlap(insertion, deletion);
            if (overlap_length1 >= overlap_length2) {
                if (overlap_length1 >= deletion.Length / 2.0 or
                    overlap_length1 >= insertion.Length / 2.0)
                {
                    // Overlap found.
                    // Insert an equality and trim the surrounding edits.
                    diffs.Insert(pointer, Diff{ .operation = .equal, .text = insertion.Substring(0, overlap_length1) });
                    diffs.items[pointer - 1].text =
                        deletion.Substring(0, deletion.Length - overlap_length1);
                    diffs.items[pointer + 1].text = insertion.Substring(overlap_length1);
                    pointer += 1;
                }
            } else {
                if (overlap_length2 >= deletion.Length / 2.0 or
                    overlap_length2 >= insertion.Length / 2.0)
                {
                    // Reverse overlap found.
                    // Insert an equality and swap and trim the surrounding edits.
                    diffs.Insert(pointer, Diff{ .operation = .equal, .text = deletion.Substring(0, overlap_length2) });
                    diffs.items[pointer - 1].operation = Operation.INSERT;
                    diffs.items[pointer - 1].text =
                        insertion.Substring(0, insertion.Length - overlap_length2);
                    diffs.items[pointer + 1].operation = Operation.DELETE;
                    diffs.items[pointer + 1].text = deletion.Substring(overlap_length2);
                    pointer += 1;
                }
            }
            pointer += 1;
        }
        pointer += 1;
    }
}
