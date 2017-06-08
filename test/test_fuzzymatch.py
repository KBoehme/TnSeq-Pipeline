minbaseoffset = 4
maxbaseoffset = 6


def fuzzy_match_beginning(pattern, genome, mismatches):
    for i in range(minbaseoffset, maxbaseoffset+1):
        chunk = genome[i: i + len(pattern)]
        # now compare chunk with pattern to see if they match with at least mismatches.
        if (chunk == pattern):
            return i + len(pattern)
    return -1


def test_fuzzy_match_beginning():
    with open('../example_data/base_offset/extended_offset.fastq') as f:
        while True:
            name = f.readline().strip()[1:]  # Strip the @ sign infront of fastq names.
            seq = f.readline().strip().upper()
            plus = f.readline().strip()
            score = f.readline().strip()
            if not name or not seq or not plus or not score:
                break  # We are done, lets break out of the loop.
            match = fuzzy_match_beginning('TCGAGATGTGTATAAGAGACAG', seq, 0)
            assert match != -1, "Didn't find a match!"
