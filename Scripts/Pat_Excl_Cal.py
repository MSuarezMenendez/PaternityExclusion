# Python script


import sys
sys.path.insert(0, './')


def read_tsv_genotype_file(data_file: str) -> tuple:
    """
    :param data_file: name and complete path for file with data in tsv format. First line with locus names \
    starting at index 1 and two labels per locus, i.e., locus_1(1)\tlocus_1(2)
    :returns: tuple with three lists; sample labels (in all upper case), marker names (also upper case), and
     diploid genotypes (list of two-index lists)
     """
    sample_ids, genotypes = list(), list()

    fp = open(data_file, 'r')
    loci = fp.readline().strip().split()[1::2]
    for index in range(len(loci)):
        loci[index] = loci[index].strip('(1)').upper()
    index = 0
    for line in fp:
        info = line.strip().split()
        sample_ids.append(info[0].upper())
        genotypes.append(list())
        for locus_index in range(1, 2 * len(loci), 2):
            genotypes[index].append([int(info[locus_index]), int(info[locus_index + 1])])
        index += 1
    fp.close()

    return sample_ids, loci, genotypes

def allele_counts(genotypes: list,
                  loci: list) -> dict:

    """Tallies the counts of each allele and its frequency
    :param genotypes: list of genotypes. Each index represents one locus containing a list
    of diploid genotypes (i.e., a two-index list, one index per allele). Order of loci the
    same as in the argument loci
    :param loci: list of marker names. Order of makers the same as in genotypes argument
    :return: dictionary indexed by marker names each comprising a dict for each allele
    with a two-index list; first index count of the allele, the second the frequency of
    that allele
    """
    allele_counts = dict()
    for locus in loci:
        allele_counts[locus] = dict()
    for locus in range(len(loci)):
        for sample in range(len(genotypes)):
            for allele in genotypes[sample][locus]:
                if allele not in allele_counts[loci[locus]] and allele not in [0, -1]:
                    allele_counts[loci[locus]][allele] = [0, 0.]
                allele_counts[loci[locus]][allele][0] += 1

    for locus in allele_counts:
        genecopy_count = 0
        for allele in allele_counts[locus]:
            genecopy_count += allele_counts[locus][allele][0]
        for allele in allele_counts[locus]:
            allele_counts[locus][allele][1] = allele_counts[locus][allele][0] / genecopy_count
        freq_sum = 0.
        for allele in allele_counts[locus]:
            freq_sum += allele_counts[locus][allele][1]

    return allele_counts

def calc_single_parent_pr_excl(allele_counts: dict) -> tuple:
    """Calculates the probability of excluding a non-parent for a sample of two individuals,
    i.e., one putative offspring and one of its putative parent

    :param allele_counts: dictionary indexed by marker names each comprising a dict for each allele
    with a two-index list; first index count of the allele, the second the frequency of
    that allele
    :return: two objects, the overall exclusion probability, and a dictionary indexed by
    marker names each with the exclusion probability of the locus
    """

    locus_excl_prob, all_loci_excl_prob = dict(), 1.
    for locus in allele_counts:
        jamieson_first, jamieson_third, jamieson_fourth = 0, 0, 0
        for allele in allele_counts[locus]:
            jamieson_first += allele_counts[locus][allele][1] ** 2
            jamieson_third += allele_counts[locus][allele][1] ** 3
            jamieson_fourth += allele_counts[locus][allele][1] ** 4

        locus_excl_prob[locus] = 1 - (4 * jamieson_first) + \
                                    (2 * jamieson_first ** 2) + \
                                    (4 * jamieson_third) - \
                                    (3 * jamieson_fourth)

        all_loci_excl_prob *= (1 - locus_excl_prob[locus])
    return ( all_loci_excl_prob), locus_excl_prob

def single_parent_excl_prob(genotypes: list=None,
                            loci: list=None,
                            data_file: str=None) -> tuple:
    """
    Calculates the probability of excluding a parent for putative offspring and single parent pairs. The input
     can be EITHER a list with genotype data and loci OR the path and name to a data file. In the first case BOTH
     lists are required. If a data file name is provided, them the genotypes and loci lists will be generated from
     the data file (and thus any genotype and/or loci list provided as an argument ignored)

    :param genotypes:  list of genotypes. Each index represents one locus containing a list
    of diploid genotypes (i.e., a two-index list, one index per allele). Order of loci the
    same as in the argument loci
    :param loci: list of marker names. Order of makers the same as in genotypes argument
    :param data_file: name and complete path for file with data in tsv format. First line with locus names \
    starting at index 1 and two labels per locus, i.e., locus_1(1)\tlocus_1(2)
    :return: two objects, the overall exclusion probability, and a dictionary indexed by
    marker names each with the exclusion probability of the locus
    """

    if data_file:
        sample_ids, loci, genotypes = read_tsv_genotype_file(data_file=data_file)

    frequencies = allele_counts(genotypes=genotypes, loci=loci)

    return calc_single_parent_pr_excl(allele_counts=frequencies)
print(single_parent_excl_prob(None,None,sys.argv[1]))
