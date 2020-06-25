from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MuscleCommandline
import sys, subprocess, re, json, gzip
from multiprocessing.pool import ThreadPool

# shift SNP positions when read has an insertion, finds insertion position and count
def shift_positions(seq, positions):
    indel_pos = seq.find('-')
    indel_count = re.search('[ATCG]', seq[indel_pos:]).start()
    lower_pos = list(filter(lambda x: x < indel_pos, positions))
    higher_pos = list(filter(lambda x: x > indel_pos, positions))
    higher_pos_shifted = [i+indel_count for i in higher_pos]
    shifted_positions = lower_pos + higher_pos_shifted
    return shifted_positions

# getting the best match for the read, by aligning all ref-strains to the read
def get_best_target_match(query, targets, positions, abundances):
    # building MUSCLE commandline and issuing command
    muscle_cline = MuscleCommandline(clwstrict=True)
    child = subprocess.Popen(str(muscle_cline), stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
    stderr=subprocess.PIPE, universal_newlines=True, shell=(sys.platform!="win32"))
    input = targets + [query]
    SeqIO.write(input, child.stdin, "fasta")
    child.stdin.close()

    # reading the alignment and determine position of the read in the alignment
    align = AlignIO.read(child.stdout, "clustal")
    seq_names = [i.id for i in align]
    query_idx = seq_names.index(query.id[:32])
    idx_offset_det = 0
    if query_idx == 0:
        idx_offset_det = 1

    # getting baspair offset of the ref-strains as they align somewhere in the middle of the read
    offset = re.search('[ACTG]',  str(align[idx_offset_det].seq)).start()
    query_snp = [align[query_idx].seq[(offset+r)-1] for r in positions]
    for i in align:
        # test if there is an insertion in the read, which shifts SNP positions in the ref-strains
        if '-' in i[offset:offset+90] and (i.description in abundances):
            shifted_positions = shift_positions(str(i.seq)[offset:], positions)
            query_snp = [align[query_idx].seq[(offset+r)-1] for r in shifted_positions]
            target_snp = [i.seq[(offset+r)-1] for r in shifted_positions]
        else:
            target_snp = [i.seq[(offset+r)-1] for r in positions]
        if (query_snp == target_snp) and (i.description in abundances):
            return i.description
    editdist = len(query_snp) - len([i for i, j in zip(query_snp, target_snp) if i == j])
    if (editdist < 5) and (i.description in abundances):
        return 'similar'

# helper function for ThreadPool, getting best match for the read and increment abundance
def process(query):
    match = get_best_target_match(query, target_strains, target_positions, abundances)
    if match:
        abundances[match] += 1
    else:
        abundances['unknown'] += 1

if __name__ == '__main__':
    if not sys.argv[1:2]:
        print("USAGE: analyse_IGR.py merged_reads.fastq IGR_strains.fasta")
        exit()

    #load the data and assign SNP positions
    zipped_file = gzip.open(sys.argv[1], "rt")
    target_strains = list(SeqIO.parse(sys.argv[2], format="fasta"))
    query_reads = SeqIO.parse(zipped_file, "fastq")
    target_positions = [4,7,10,12,17,18,20,16,33,34,42,45,54,64,65,76,80,89]
    abundances = dict.fromkeys([i.name for i in target_strains] + ['unknown', 'similar'], 0)

    # initiate ThreadPool and process all reads
    pool = ThreadPool()
    pool.map(process, query_reads)

    # save abundances
    print(abundances)
    f = open(sys.argv[1] + "_abundance.json", "w")
    f.write(json.dumps(abundances))
    f.close()
