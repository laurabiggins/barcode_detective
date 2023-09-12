#!/usr/bin/env python
import re

# reverse complement a sequence
def reverseComplement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '_': '_'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

# check if it's a dual barcode
def isDual(seq):
    if "_" in seq:
        return True

# get a list of the 2 barcode lengths
def getDualLengths(seq):
    if not isDual(seq):
        print(f"This isn't a dual barcode - why are you running this function?")
    else:
        seq_lengths = list(map(len, seq.split("_")))
        return seq_lengths

def checkPolyG(observed_seq):
    if isDual(observed_seq):
        seq_length = min(getDualLengths(observed_seq)) # use minimum length of the dual barcodes
    else:
        seq_length = len(observed_seq)
    pattern = 'G{' + str(seq_length-1) + ',}' # searching for Gs the full length of barcode-1 
    #print(f"pattern: {pattern}") 
    p = re.compile(pattern)
    m = p.search(observed_seq)

    if m is not None:
        print (f"found polyG: {observed_seq}")


def checkPhiX(observed_seq):
    pattern = 'GTATGCCG'
    p = re.compile(pattern)
    m = p.search(observed_seq)

    if m is not None:
        print (f"found phiX: {observed_seq}")

    # Also check rev comp of sequence
    pattern = 'CGGCATAC'
    p = re.compile(pattern)
    m = p.search(observed_seq)

    if m is not None:
        print (f"found phiX reverse complemented: {observed_seq}")


def reverseSeq(seq):
    return seq[::-1]

def complementSeq(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '_': '_'}
    complemented_seq = "".join(complement.get(base, base) for base in seq)
    return complemented_seq

def isReverse(observed_seq, expected_seq):
    reversed_seq = reverseSeq(observed_seq)
    if expected_seq == reversed_seq:
        return (f"seq is reversed")

def isComplemented(observed_seq, expected_seq):
    complemented_seq = complementSeq(observed_seq)
    if expected_seq == complemented_seq:
        return (f"seq is complemented")

def isReverseComplemented(observed_seq, expected_seq):
    rev_comp = reverseComplement(observed_seq)
    if rev_comp == expected_seq:
        return (f"seq is reverse complemented") 

# Check if one sequence is found within another one
def isSubString(seq1, seq2):

    seq1in2 = re.search(seq1, seq2)
    seq2in1 = re.search(seq2, seq1)

    # if first sequence has been found within the second sequence
    if seq1in2 is not None:
        print (f"found seq1 {seq1} in seq2 {seq2}")
        m = seq1in2

        if m.span()[0] > 0:
            first_matched_pos = m.span()[0]
            extra_bases_at_start = seq2[:first_matched_pos]
            print(f"Extra bases at start of sequence, these are {extra_bases_at_start}")

        if m.span()[1] < len(seq2):
            last_matched_pos = m.span()[1]
            extra_bases_at_end = seq2[last_matched_pos:]
            print(f"Extra bases at end of observed sequence, these are {extra_bases_at_end}")

    # if second sequence has been found within the first sequence
    elif seq2in1 is not None:
        print (f"found seq2 {seq2} in seq1 {seq1}")
        m = seq2in1

        if m.span()[0] > 0:
            first_matched_pos = m.span()[0]
            extra_bases_at_start = seq1[:first_matched_pos]
            print(f"Extra bases at start, these are {extra_bases_at_start}")

        if m.span()[1] < len(seq1):
            last_matched_pos = m.span()[1]
            extra_bases_at_end = seq1[last_matched_pos:]
            print(f"Extra bases at end, these are {extra_bases_at_end}")


# Check if one sequence is found within another one
# Extended version to check revComp etc
def isSubStringExtended(shorter_seq, longer_seq):

    if (len(shorter_seq) > len(longer_seq)):
        print(f"The first sequence needs to be longer than the second")

    else:
        print("================")
        info_msg = ""
        m = re.search(shorter_seq, longer_seq) # this only works by finding the smaller seq within the longer seq  

        if m is not None:
            info_msg = f"found shorter_seq {shorter_seq} in longer_seq {longer_seq}"
            print (info_msg)
            # we didn't find a straightforward match, try some other options
        else: # try the reverse complement of the sequence
            m = re.search(reverseComplement(shorter_seq), longer_seq)  

            if m is not None:
                info_msg = f"found reverse comp of shorter_seq {shorter_seq} - {reverseComplement(shorter_seq)} in longer_seq {longer_seq}"
                print (info_msg)

            else: # try reversing the sequence
                m = re.search(reverseSeq(shorter_seq), longer_seq)  

                if m is not None:
                    info_msg = f"found reverse of shorter_seq {shorter_seq} - {reverseSeq(shorter_seq)} in longer_seq {longer_seq}"
                    print (info_msg)  

                else: # try reversing the sequence
                    m = re.search(complementSeq(shorter_seq), longer_seq)  

                if m is not None:
                    info_msg = f"found complement of shorter_seq {shorter_seq} - {complementSeq(shorter_seq)} in longer_seq {longer_seq}"
                    print (info_msg)  


        if m is not None:
            if m.span()[0] > 0:
                first_matched_pos = m.span()[0]
                extra_bases_at_start = longer_seq[:first_matched_pos]
                print(f"Extra bases at start of sequence, these are {extra_bases_at_start}")

            if m.span()[1] < len(longer_seq):
                last_matched_pos = m.span()[1]
                extra_bases_at_end = longer_seq[last_matched_pos:]
                print(f"Extra bases at end of sequence, these are {extra_bases_at_end}")
        else:
            print(f"Didn't find a match - this is still unexplained")



# checking for a single typo
def isOneMismatch(observed_seq, expected_seq):
    obs = list(observed_seq)
    exp = list(expected_seq)
    if len(obs) == len(exp):
        diff = False
        for i in range(0, len(obs), 1):
            if (obs[i] != exp[i]):

                # If first mismatch
                if (diff == False):
                    diff = True

                # Second mismatch
                else:
                    diff = False
                    break
        if diff == True:
            return(f"One mismatch between the 2 sequences - potential typo")




# checkPolyG(observed_seq1)
# checkPolyG(observed_seq2)

# checkPhiX(observed_seq1)
# checkPhiX(observed_seq2)
# checkPhiX(observed_seq3)
# checkPhiX(reverseComplement(observed_seq3))

print(isReverse("AGCTCGAC", "AGCGGTAC"))
print(isReverse("AGCTCGAC", "CAGCTCGA")) # this is reversed

# isReverseComplemented("AGCTCGAC", "GTCGAGCT") # this is reverse complemented
# isComplemented("AGCTCGAC", "GTCGAGCT")
# isComplemented("AGCTCGAC", "TCGAGCTG") #  this is complemented

#isOneMismatch("AGCTCGAC", "ATCGGTAC") # this is one mismatch

isSubStringExtended("TAACCAAG", "ACCAAG")
#isSubString("TAACCAAG", "AACCA")
isSubStringExtended("AACCA", "TAACCAAG")

isSubStringExtended("TGACGT", "TAACGTCAG")
isSubStringExtended("ATTGCAGT", "TAACGTCAG")
isSubStringExtended("ACTGCA", "TAACGTCAG")