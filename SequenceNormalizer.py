__author__ = 'odergai'
"""SequenceNormalizer takes as an argument string containing proteins sequence
    and return length of minimal sub-sequence that needs to be substituted
    to have all amino acids equally represented.
    Args:
        :param inseq : input protein sequence
        :param verbose : controls output format
            If parameter verbose is True or T prints a table to STDOUT:
                start   end    length
                <int>   <int>   <int>
                otherwise prints just a length of substring to STDOUT.
    Prints:
        tuple (start, end, length) of the sun-sequence to be substituted.
    """

import re
import argparse
from operator import itemgetter

def get_shortest_subset(
        aminoacid='',
        aas_index_dict={},
        target_length=0):
    """Find the shortest sub-sequence which accommodates all excessive residues of a given type.
        Args:
            :param aminoacid : one letter code identifier (default '').
            :param aas_index_dict : dictionary: keys are amino acids and values are list of their indexes (default {}).
            :param target_length : indicates number or amino acids to be kept in protein (default 0).
        Returns:
                a tuple (start, end, length)
        """
    index_list = aas_index_dict[aminoacid]
    subseq_length = len(index_list) - target_length
    iterate_lim = len(index_list) - subseq_length
    distance_list = []
    for i in range(0, iterate_lim):
        subseq_distance = index_list[i + subseq_length - 1] - index_list[i] + 1
        distance_list.append((index_list[i], index_list[i + subseq_length - 1], subseq_distance))
    out = sorted(distance_list, key=itemgetter(2))
    return out[0]


def rearrange_sequence(inseq=''):
    """Find  the shortest sub-sequence which accommodates all excessive residues
    of all types
    Args:
        :param inseq : input protein sequence (default '').
    Returns:
        tuple (start, end, length) of the shortest sub-sequence
    """
    standard_aminoacids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    
    if len(inseq) == 0:
        print('Length of input sequence is 0, check input sequence')
        exit()
        
    # checks if length of input sequence can accommodate equal number of amino acids
    aas_list = list(set(inseq.upper()))
    for residue in aas_list:
        if residue not in standard_aminoacids:
            print('Check sequence: {} is non a standard amino acid notation'.format(residue))
            exit()
    if len(inseq) % len(aas_list) > 0:
        print('no possible rearrangement for protein of length {0}'.format(len(inseq)))
        exit()
    else:
        target_length = int(len(inseq) / len(aas_list))

    aas_index_dict = {}  # {'aa':[ind1 ... indN]}
    target_dict = {}  # {'aa':number of residues to be substituted,}
    subseq_dict = {}  # {'aa':(start, end, distance) }
    over_rep_aas_indexes = []  # list of indexes of overrepresented amino acids
    for aa in aas_list:
        if inseq.count(aa) > target_length:
            aas_index_dict[aa] = [ind.start() for ind in re.finditer(aa, inseq)]
            target_dict[aa] = len(aas_index_dict[aa]) - target_length
            over_rep_aas_indexes = over_rep_aas_indexes + aas_index_dict[aa]
            subseq_dict[aa] = get_shortest_subset(aminoacid=aa,
                                                  aas_index_dict=aas_index_dict,
                                                  target_length=target_length)

    min_dist_list = [subseq_dict[aa] for aa in sorted(subseq_dict,
                                                      key=lambda aa: subseq_dict[aa][2],
                                                      reverse=True)]
    min_dist = min_dist_list[0][2]
    threshold = len(aas_index_dict)  # number of overrepresented amino acids
    over_rep_subset = [i for i in over_rep_aas_indexes if i <= (len(inseq) - min_dist)]
    over_rep_subset = sorted(over_rep_subset)

    aa_to_search = [aa for aa in sorted(aas_index_dict,
                                        key=lambda aa: len(aas_index_dict[aa]),
                                        reverse=True)]


    # This block of code compares indexes of all overrepresented amino acids
    # with position of the shortest sub-sequence for the amino acid to balance which
    # one needs to substitute the longest sub-sequence

    keep_min_dist = True
    for k in aa_to_search:
        hit = len([i for i in aas_index_dict[k] if min_dist_list[0][0] <= i <= min_dist_list[0][1]])
        if hit < target_dict[k]:
            keep_min_dist = False
            break
    if keep_min_dist:               # if all other excessive amino acids will be substituted in this range
        return min_dist_list[0]     # return the shortest sub-sequence

    # Scans sequence with sliding window of length min_dist
    extra = 0
    min_seq_found = False
    while min_seq_found is False:
        for ind in over_rep_subset:
            fraglen = (ind, ind + min_dist - 1 + extra)
            hits_list = []
            for k in aa_to_search:
                hit = len([i for i in aas_index_dict[k] if fraglen[0] <= i <= fraglen[1]])
                if hit < target_dict[k]:
                    break
                elif hit >= target_dict[k]:
                    hits_list.append(k)
            if len(hits_list) >= threshold:
                min_seq_found = True
                out = (fraglen[0], fraglen[1], fraglen[1] - fraglen[0] + 1)
                break
        extra += 1
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("inseq", type=str, help="Enter input  protein sequence,"
                                                " Please enter protein sequence in single letter code")
    parser.add_argument('-v', '--verbose', choices=['True', 'T', 'False', 'F'],
                        default='True',
                        help='choice to output only length of sub-sequence or start, end and length')
    args = parser.parse_args()
    output = rearrange_sequence(args.inseq)
    if args.verbose in ['True', 'T']:

        print('start{0}end{0}length'.format('\t'))

        print('\t'.join(str(x) for x in output))
    else:
        print(output[2])

