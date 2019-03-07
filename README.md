# SequenceNormalizer
SequenceNormalizer.py script takes as an argument string containing proteins sequence and returns length of the minimal sub-sequence that needs to be substituted to have all amino acids equally represented.
SequenceNormalizer can be imported as a module, in this case rearrange_sequence() function should be used. 
rearrange_sequence() returns a tuple (start, end, length) of the shortest sub-sequence.

Important! start and end coordinates are given in 0-based form.
SequenceNormalizer.py normalizes counts of amino acids initially found in the input sequence without trying to fit all 20 canonical amino acids. 
