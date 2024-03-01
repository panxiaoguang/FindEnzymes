import argparse
from pyfaidx import Fasta
from Bio.Data import IUPACData

def reverse_complement(dna_seq):
  complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
  return "".join([complement_dict[base] for base in dna_seq[::-1]])

def parse_coordinate(coord):
    chromosome, pos = coord.split(":")
    pos = int(pos)
    return chromosome,pos

def get_enzyme_seq(enzyme_file="REs.tsv"):
    enzyme_dct = {}
    with open(enzyme_file) as f:
        for line in f:
            if line.startswith("Name"):
                continue
            name, seq = line.strip().split("\t")
            enzyme_dct[name] = seq
    return enzyme_dct

def parse_meta(extrainfo):
    if extrainfo.startswith("p"):
        pos,seq,orient = extrainfo.split("|")
        pos = pos.replace("p","")
        pos = int(pos)
        return pos,"+"
    else:
        pos,seq,orient = extrainfo.split("|")
        pos = pos.replace("m","")
        pos = int(pos)
        return pos,"-"


def get_sequence_around_position(fasta_file, chromosome, position, upstream=120, downstream=120):
    genome = Fasta(fasta_file)
    start = max(0, position - upstream - 1)  # 修正为确保上游长度为120bp
    end = position + downstream  # 修正为确保下游长度为120bp
    sequence = str(genome[chromosome][start:end])
    return sequence

def get_aligned_sequence(center,pos,strand,enzymeSeq):
    seq = ["-" for _ in range(242)]
    sitelength = len(enzymeSeq)
    if strand == "+":
        cord = 120-(center-pos-1)
        seq[(cord-1):(cord+sitelength)] = enzymeSeq
    else:
        cord = 120-(center-pos-1)
        seq[(cord-sitelength):cord] = reverse_complement(enzymeSeq)
    return "".join(seq)

def are_bases_compatible(base1, base2):
    base1_set = set(IUPACData.ambiguous_dna_values[base1.upper()])
    base2_set = set(IUPACData.ambiguous_dna_values[base2.upper()])
    return not base1_set.isdisjoint(base2_set)

def annotate_diff(s1,s2):
    result = []
    for i in range(len(s1)):
        if s2[i] == "-":
            result.append(" ")
        elif are_bases_compatible(s1[i],s2[i]):
            result.append("|")
        else:
            result.append("*")
    return "".join(result)

def main():
    parser = argparse.ArgumentParser(description="Fetch a sequence around a given position with the position highlighted.")
    parser.add_argument("coordinate", help="Chromosome coordinate in a format like: chr2:123456")
    parser.add_argument("enzyme", help="The file of the enzyme")
    parser.add_argument("name", help="The name of the enzyme")
    parser.add_argument("meta", help="meta info")
    parser.add_argument("matched", help="The name of the enzyme")
    args = parser.parse_args()

    chrom, pos = parse_coordinate(args.coordinate)
    enzymeSeq = get_enzyme_seq(args.enzyme)
    pos2,orient = parse_meta(args.meta)
    sequence = get_sequence_around_position("../cast-seq-pipeline/annotations/human/bowtie2Index/genome.fa", chrom, pos)
    aligned_seq = get_aligned_sequence(pos,pos2,orient,enzymeSeq[args.name])
    diff = annotate_diff(sequence,aligned_seq)
    print(sequence)
    print(diff)
    print(aligned_seq)
    print("second enzyme seq:",enzymeSeq[args.matched])



if __name__ == "__main__":
    main()


