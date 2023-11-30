import sys
import collections
import pysam
from math import log2, log10
import argparse

# log likelihood of a given genotype
# g is number of reference alleles in the individual
# https://arxiv.org/abs/1203.6372
def pl_const_err(g, ref_count, alt_count, e=0.001):
    #        2^-k     lh of ref allele               lh of alt allele
    # return  -(rc+ac) + (rc*log2( (2-g)*e + g*(1-e) )) + (ac*log2( (2-g)*(1-e) + g*e ))
    log_sum = -(ref_count + alt_count)*log10(2)
    log_sum += ref_count * log10( (2-g)*e + g*(1-e) )
    log_sum += alt_count * log10( (2-g)*(1-e) + g*e )
    return round(-10*log_sum, 0)


def gl_to_pl(gl):
    return -10*log10(gl)

def read_fai(fai):
    id_to_contig = {}
    for i, line in  enumerate(open(fai, "r")):
        fields = line.strip().split()
        id_to_contig[i] = fields[0]
    return id_to_contig


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf")
    parser.add_argument("seeds")
    parser.add_argument("--fasta_fai")
    parser.add_argument("--sample_name", default="sample")
    args = parser.parse_args()
    vcf = pysam.VariantFile(args.vcf)
    vcf.subset_samples([])
    sample_name = args.sample_name
    m_counts = {}
    id_to_contig = None
    if args.fasta_fai:
        id_to_contig = read_fai(args.fasta_fai)
    for line in open(args.seeds):
        fields = line.strip().split()
        read = fields[0]
        # TODO: apply filters here?
        if fields[5] == ".":
            continue
        for marker in fields[5].split(","):
            s, p, ale = list(map(int, marker.split("/")))
            if args.fasta_fai:
                contig = id_to_contig[s]
                s = vcf.get_tid(contig)
            if (s, p) not in m_counts:
                m_counts[(s,p)] = {}
            if ale not in m_counts[(s,p)]:
                m_counts[(s,p)][ale] = 0
            m_counts[(s,p)][ale] += 1
    sorted_ms = sorted(m_counts.keys())
    # modify header to add/remove samples
    out_vcf = pysam.VariantFile("-", "w")
    out_vcf.header.add_sample(sample_name)
    for contig in vcf.header.contigs.keys():
        out_vcf.header.contigs.add(contig, length=vcf.header.contigs[contig].length)
    out_vcf.header.info.add("AC", "A", "Integer","Allele count in genotypes for each ALT allele, in the same order as listed")
    out_vcf.header.info.add("DP", 1, "Integer", "Total Depth")
    out_vcf.header.formats.add("GT", 1, "String", "Genotype")
    out_vcf.header.formats.add("GQ", "1", "Integer","Genotype quality (phred-scaled 2nd highest PL)")
    out_vcf.header.formats.add("PL", "G", "Integer","Phred-Scaled Genotype Likelihood")
    for rec in vcf:
        if (rec.rid, rec.start) not in m_counts:
            continue
        count_dict = m_counts[(rec.rid, rec.start)]
        alt_keys = [x for x in sorted(count_dict.keys()) if x != 0]
        alt_allele = alt_keys[0] if len(alt_keys) else 1
        if len(alt_keys) > 1: # more than one alt allele
            alt_allele = max([(k, count_dict[k]) for k in alt_keys], key=lambda x: x[1])[0]
        if 0 not in count_dict:
            count_dict[0] = 0
        if alt_allele not in count_dict:
            count_dict[alt_allele] = 0
        if alt_allele+1 > len(count_dict) or max(count_dict.keys()) > len(rec.alleles)-1:
            continue
        gs_pls = [(g, pl_const_err(g, count_dict[0], count_dict[alt_allele], e=0.001)) for g in range(3)]
        call = min(gs_pls, key=lambda x: x[1])[0]
        min_pl = min(gs_pls, key=lambda x: x[1])[1]
        scaled_pls = [pl-min_pl for g, pl in gs_pls]
        gq = min(p for p in scaled_pls if p != 0)
        new_rec = out_vcf.new_record(contig=rec.contig, id=rec.id, start=rec.start, stop=rec.stop, alleles=rec.alleles, filter=rec.filter)
        gt = None
        if call == 1:
            gt = (0, alt_allele)
        elif call == 0:
            gt = (alt_allele, alt_allele)
        elif call == 2:
            gt = (0,0)
        new_rec.samples[sample_name]["GT"] = gt
        new_rec.samples[sample_name]["GQ"] = gq
        pls = tuple([int(pl-min_pl) for g, pl in sorted(gs_pls, key=lambda x:x[0], reverse=True)])
        if len(new_rec.alleles) > 2:
            x = len(new_rec.alleles)-1
            x = int((x*(x+1)/2)+x)+1
            new_pls = [255]*x
            new_pls[0] = pls[0]
            new_pls[int((alt_allele*(alt_allele+1)/2)+0)] = pls[1]
            new_pls[int((alt_allele*(alt_allele+1)/2)+alt_allele)] = pls[2]
            pls = tuple(new_pls)
        new_rec.samples[sample_name]["PL"] = pls
        if len(rec.alleles) > len(count_dict):
            for i in range(len(rec.alleles)):
                if i not in count_dict:
                    count_dict[i] = 0
        ac = tuple(list(count_dict[a] for a in sorted(count_dict.keys()) if a != 0))
        new_rec.info["AC"] = ac
        new_rec.info["DP"] = sum(x[1] for x in count_dict.items())
        out_vcf.write(new_rec)
    vcf.close()
    out_vcf.close()
