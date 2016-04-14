#!/usr/bin/env python

from utility import FormatError
import fileinput, sys, string, sequence, math, random

def read_prob_matrix(filehandle, alphabet):
    """Read a matrix from the specified input source, and return its
    transpose, as a list of lists. Each row of the input matrix must contain
    an entry for each letter in the alphabet, and the entries in a row must
    sum to 1, to within two decimal places (ie > 0.99 and < 1.01). In this
    manner, the matrix must adhere to the format of a PWM.
    Raises a SyntaxWarning if the file does not adhere to that format."""
    
    # alphabet must be protein, rna, or dna alphabet:
    assert ((alphabet == sequence.dnachar) or \
            (alphabet == sequence.rnachar) or (alphabet == sequence.protchar))

    matrix = []
    # For each line (ie row) of the file, convert it to a column in the final
    # list of lists:
    for line in filehandle:
        curr_col = []
        curr_col_sum = 0
        letter_probs = line.split()

        # Line must contain exactly one entry for each letter in alphabet:
        if (not (len(letter_probs) == len(alphabet))):
            raise SyntaxWarning(\
            "Length of line was not length of alphabet: " + line)
        for prob_token in letter_probs:
            curr_prob = None
            # Each line entry must be a probability:
            try:
                curr_prob = float(prob_token)
                # Must be a probability => >= 0 and <= 1
                if ((curr_prob < 0) or (curr_prob > 1.01)):
                    raise SyntaxWarning(\
                    "Probability was invalid: " + str(curr_prob))
            except ValueError:
                raise SyntaxWarning(\
                "Invalid probability token: " + prob_token)

            # Add the probability to the column:
            curr_col.append(curr_prob)
            curr_col_sum = curr_col_sum + curr_prob

        # Line entries must sum to 1, to within two decimal places:
        if (not ((0.99 < curr_col_sum) and (curr_col_sum < 1.01))):
            raise SyntaxWarning("Row did not sum to 1: " + line)

        # Re-normalise the column to sum to 1, incase the data in the
        # input file was not precise enough to sum to 1 exactly:
        curr_col_sum = sum(curr_col)
        curr_col = map(lambda x:x/curr_col_sum, curr_col)
        matrix.append(curr_col)

    return matrix
            
    
def read_prob_matrix_unrotated(filehandle, col_size, entry_type):
    """Read a matrix from the specified input source, and return it
    as a list of lists. Each line of the input file specifies a row of the
    matrix. Each column of the input file must be exactly <col_size> in length.
    Raises a SyntaxWarning if the file does not adhere to this format."""
    
    assert (col_size >= 0)

    n_cols = None # Number of columns in matrix will be determined when
                  # first line is parsed.
    matrix = None

    n_lines = 0
    
    # For each line (ie row) of the file, convert it to a row in the final
    # list of lists:
    for line in filehandle:
        if filehandle.isfirstline():
            # Determine the number of columns in the matrix by looking at the
            # first line:
            n_cols = len(line.split())
            
            # Initiate the matrix given that information:
            matrix = []
            for col_idx in range(n_cols):
                matrix.append([])

        n_lines += 1 # Count the number of lines in the file
        curr_row = []
        row_entries = line.split()

        # Line must contain exactly n_cols columns (ie all rows in the
        # file must have the same number of elements):
        if (not (len(row_entries) == n_cols)):
            raise SyntaxWarning(\
            "Rows in file do not have consistent number of columns.")

        curr_col_idx = 0
        for token in row_entries:
            try:
                mat_entry = entry_type(token)
            except ValueError:
                raise SyntaxWarning(\
                "Invalid matrix entry: " + token)

            # Add the matrix entry to the end of the appropriate column:
            curr_col = matrix[curr_col_idx]
            curr_col.append(mat_entry)
            curr_col_idx += 1

    # Total number of columns must equal the number of entries there are
    # supposed to be in each column:
    if (not (n_lines == col_size)):
        raise SyntaxWarning("Incorrect number of rows in file.")

    return matrix


def insert_site(site, seq, offset=None):
    """Inserts a sequence (represeting a site) into a larger sequence (which
    is a sequence object rather than a series of letters."""
    # inputs:
    #  site               The site to be inserted
    #  offsets            the offset where the site is to be inserted
    #  seq                The sequence into which the specified site is to
    #                     be implanted

    # get sequence info
    name = seq.getName()
    seq_data = seq.getSeq()

    assert ((offset == None) or ((offset >= 0) and \
                                 (offset <= (len(seq_data) - len(site)))))
    
    # select a random offset if none given:
    if (offset == None):
        # insert signal in a random position, from 0 up to m (= l - w)
        offset = random.randint(0,(len(seq_data) - len(site)))
        
    # insert the signal
    signal_seq = seq_data[:offset]+str(site)+seq_data[(offset + len(site)):]
    
    # create a modified sequence object to return
    new_seq = sequence.Seq(name, signal_seq)
    
    return new_seq


def get_alph_idx(letter, alph):
    """Return the index of the specified letter in the specified alphabet.
    Presently, does this by iterating over the letters in alph. Alteratively, I
    could restrict alph to certain types, and could then retrieve the letter
    through some kind of O(1) mapping."""
    assert (letter in alph)
    for lett_idx in range(len(alph)):
        if (letter == alph[lett_idx]):
            return lett_idx
    # As a precondition to this function, letter must be an element of alph =>
    # Not valid to reach the end of alph without finding letter:
    print >> sys.stderr, "INVALID LETTER FOR GET_ALPH_IDX: ", letter
    print "alph is: ", alph
    assert(False)


def get_alignment_counts(seq_align, alph):
    """Convert a list of aligned sequence (of identical length) into a matrix
    of letter counts (letters taken from the specified alphabet)."""

    seq_len = len(seq_align[0])
    count_matrix = []
    # Generate an empty letter count column for each column of the aligned
    # sequences:
    for nuc_idx in range(seq_len):
        col = []
        for lett in range(len(alph)):
            col.append(0)
        count_matrix.append(col)

    for seq in seq_align:
        assert (len(seq) == seq_len)

        # For each letter in the current sequence, increment the count of the
        # current letter at the current column of the count matrix:
        for nuc_idx in range(len(seq)):
            # Get the index of the current letter in the alphabet:
            lett_idx = get_alph_idx(seq[nuc_idx], alph)

            # Increment the count of that letter at the current index in the
            # count matrix:
            count_matrix[nuc_idx][lett_idx] += 1

    return count_matrix
            

def convert_to_freq(matrix):
    """Converts a matrix from a count matrix to a frequency matrix. Each
    entry in a given column will be divided by the sum of that column."""

    freq_matrix = []
    for col in matrix:
        freq_col = []

        # Convert the column from counts to frequencies...

        try:
            # Find the sum of the column entries:
            col_sum = reduce(lambda x, y: x+y, col)

            # For each entry in the count column, divide it by the sum
            # of the column and then append it to the current frequency
            # column:
            for count in col:
                freq = float(count)/col_sum
                freq_col.append(freq)
        except TypeError, e:
            raise SyntaxWarning("Invalid count matrix column")
        except ValueError, e:
            raise SyntaxWarning("Invalid count matrix column")

        # Add the current column of frequencies to the end of the accumulating
        # frequency matrix:
        freq_matrix.append(freq_col)
    return freq_matrix


def get_lett_freqs(seqs, alph):
    """Estimate the frequencies of the letters in alph, by scanning the specified
    list of sequences"""
    lett_counts = [0]*len(alph)
    total_count = 0
    for seq in seqs:
        for lett in seq:
            lett_idx = get_alph_idx(lett, alph)
            lett_counts[lett_idx] += 1
            total_count += 1
    lett_freqs = map(lambda x: x/float(total_count), lett_counts)
    return lett_freqs


def calc_IC(column, back_freqs=None):
    """Calculates the information content of a column of frequencies (which
    must sum to 1), using the specified background frequencies. Default
    background frequencies are uniform across all in the alphabet (ie over all
    column indeces)."""

    if (back_freqs == None):
        back_freqs = [1.0/len(column)]*len(column)

    alen = len(column)

    assert (alen >= 1)

    # Calculate IC using formula:
    # IC =
    #      Sum(lett = 0 to alen):
    #         frequency(lett)*log_2(frequency(lett)/back_freq(lett))

    IC = 0

    for lett_idx in range(len(column)):
        lett_freq = column[lett_idx]
        lett_IC = None
        if (lett_freq == 0):
            lett_IC = 0
        else:
            lett_IC = lett_freq*math.log((lett_freq/back_freqs[lett_idx]),2)
        IC = IC + lett_IC
    return IC


def parse_alignment(filehandle, alphabet):
    """Read an alignment of sequences in fasta format from the specified
    file. Return an array containing each of the sequences."""

    # FIXME: Currently, all sequences must be of the same length, and
    # provision for wildcard characters has not been made. ALSO: Currently,
    # I do not check to ensure that all letters in the sequences come from
    # the specified alphabet.
    # PERHAPS I SHOULD USE THE EXISTING FASTA-PARSING FUNCTION AVAILABLE IN
    # sequence.py?

    # alphabet must be protein, rna, or dna alphabet:
    assert ((alphabet == sequence.dnachar) or \
            (alphabet == sequence.rnachar) or (alphabet == sequence.protchar))

    # Parse the file:
    aligned_sites = []
    parsed_first_seq = False
    seq_length = None # Will record length of first sequence...
    for line in filehandle:
        if (line[0] != ">"):
            if (not parsed_first_seq):
                parsed_first_seq = True
                curr_seq = line.split()[0]
                seq_length = len(curr_seq) # Record length of this sequence
                aligned_sites.append(curr_seq)
            else:
                curr_seq = line.split()[0]
                if (not (len(curr_seq) == seq_length)):
                    raise FormatError(
                        "Sequences in alignment must all be same length")
                aligned_sites.append(curr_seq)
    return aligned_sites


def invert_seq(sequence):
    """Returns the inverse complement of the specified dna sequence."""

    ic = ""
    for lett in sequence:
        ic = invert_char(lett) + ic
    return ic


def invert_char(dnalett):
    """Returns the dna character complementary to that specified."""

    assert(dnalett in sequence.dnachar)
    if (dnalett == "A"):
        return "T"
    elif (dnalett == "C"):
        return "G"
    elif (dnalett == "G"):
        return "C"
    else:
        return "A"
