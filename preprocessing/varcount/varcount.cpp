#include <cstdio>
#include <vector>
#include <string>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <getopt.h>
#include "hts_util.hpp"

struct VcntArgs {
    enum GT_TYPE {GTN, GTL, GTT, GTA};
    std::string vcf_fname = "";
    std::string sam_fname = "";
    std::string sample_name = "sample";
    int thres = 0;
    int gt = GT_TYPE::GTN;
    int keep = 0;
    int verbose = 0;
    int diploid = 0;
    int alt_default = 0;
    int min_ac = 0;
    int min_rc = 0;
    int min_c = 0;
    int gl = 0;
    double e = 0.004;
};

static inline std::tuple<std::string,std::string> truncate_str_pair(const std::string& s1, const std::string& s2) {
    return std::forward_as_tuple(s1.substr(s1.size()-1), s2.substr(s1.size()-1));
}


static inline bool var_match(const hts_util::Var& lv, const hts_util::Var& rv) {
    if (lv.type != rv.type) return false;
    switch (lv.type) {
        case hts_util::VTYPE::V_INS:
            return lv.pos+static_cast<int32_t>(lv.ref.size())-1 == rv.pos && truncate_str_pair(lv.ref, lv.alt) == truncate_str_pair(rv.ref, rv.alt) ;
        case hts_util::VTYPE::V_DEL:
            return lv.pos == rv.pos && truncate_str_pair(lv.alt, lv.ref) == truncate_str_pair(rv.alt, rv.ref);
        case hts_util::VTYPE::V_SNP:
            return lv.pos == rv.pos && lv.alt == rv.alt;
        default:
            fprintf(stderr, "no support for non SNPs & INDELs yet\n");
            return false;
    }
}


std::array<int32_t, 2> gt_by_threshold(const VcntArgs& args, const hts_util::Var& v) {
    std::array<int32_t, 2> gts;
    gts[0] = bcf_gt_missing;
    if (args.thres < 0 && v.ad[1]) {
        gts[0] = bcf_gt_phased(v.ad[1] != 0);
    } else if ((v.ad[1] || v.ad[0]) && std::abs(v.ad[0] - v.ad[1]) >= args.thres) {
        if (v.ad[0] < v.ad[1]) {
            gts[0] = bcf_gt_phased(1);
        } else if (v.ad[0] > v.ad[1])  {
            gts[0] = bcf_gt_phased(0);
        } else { // v.ad[0] == v.ad[1]
            gts[0] = bcf_gt_phased(args.alt_default);
        }
    }
    gts[1] = gts[0];
    return gts;
}

std::array<int32_t, 2> gt_by_alt_evidence(const VcntArgs& args, const hts_util::Var& v) {
    std::array<int32_t, 2> gts;
    if (v.ad[1]) {
        gts[0] = bcf_gt_phased(1);
    } else if (v.ad[0]) {
        gts[0] = bcf_gt_phased(0);
    } else {
        gts[0] = bcf_gt_missing;
    }
    gts[1] = gts[0]; // force homozygous
    return gts;
}

std::array<int32_t, 2> gt_by_likelihood(const VcntArgs& args, const std::array<double, 3>& gls) {
    std::array<int32_t, 2> gts;
    int max = gls[0], maxidx = 0;
    if (gls[1] > max) { max = gls[1]; maxidx = 1; }
    if (gls[2] > max) { max = gls[2]; maxidx = 2; }
    switch (maxidx) {
        case 0:
            gts[0] = bcf_gt_unphased(0); gts[1] = bcf_gt_unphased(0);
            break;
        case 1:
            gts[0] = bcf_gt_unphased(1); gts[1] = bcf_gt_unphased(0);
            break;
        case 2:
            gts[0] = bcf_gt_unphased(1); gts[1] = bcf_gt_unphased(1);
            break;
        default:
            fprintf(stderr, "something went wrong\n");
            exit(1);
            break;
    }
    return gts;
}

#include <htslib/sam.h>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include "mdparse.hpp"

void aln_varcount(const bam1_t* rec, hts_util::pos2var_map<hts_util::Var>* vmap) {
    uint32_t* cigar = bam_get_cigar(rec);
    int r = rec->core.pos;
    uint8_t* md_ptr = bam_aux_get(rec, "MD");
    char* md;
    if (md_ptr == NULL) {
        fprintf(stderr, "error: MD tag not found\n");
        exit(1);
    }
    else md = bam_aux2Z(md_ptr);
    auto mdv = md_parse(md);
    int mdi = 0;
    int nc = rec->core.n_cigar;
    // fprintf(stderr, "%s\n", bam_get_qname(rec));
    for (int i = 0, a=0; i < rec->core.n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int l = bam_cigar_oplen(cigar[i]);
        // load variants from p to p+l
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int nop = i+1 < rec->core.n_cigar ? bam_cigar_op(cigar[i+1]) : -1;
            for (int j = 0; j < l; ++j) {
                int apos = a + j; // apos is position within alignment
                int rpos = r + j; // rpos is position within reference
                bool mdm=0, vm=0;
                // skip deletions in md tag, we don't need them
                while (mdi<mdv.size() && (mdv[mdi].st == MD_DEL || mdv[mdi].p+rec->core.pos < rpos))
                    ++mdi;
                // check if SNV at this position
                char c = ' ';
                mdm = (mdv[mdi].st == MD_SNP && mdv[mdi].p+rec->core.pos == rpos);
                if (mdm) c = seq_nt16_str[bam_seqi(bam_get_seq(rec), apos)];
                // check var at this position
                auto found = vmap->find(rpos);
                if (found != vmap->end()) {
                    for (auto& vv: found->second) {
                        if (vv.type == hts_util::VTYPE::V_SNP) {
                            if (mdm) vv.ad[1] += (vv.alt[0] == c);
                            else ++vv.ad[0];
                        } else if (vv.type == hts_util::VTYPE::V_DEL) {
                            if (vv.alt.size() > 1)  {
                                fprintf(stderr, "%d warning: deletion alt allele has size > 1 (%d). Did you left-normalize?\n", vv.pos, vv.alt.size());
                            }
                            if (j != l-1 || (nop != BAM_CDEL && nop != BAM_CREF_SKIP )) {
                                ++vv.ad[0];
                            }
                        } else if (vv.type == hts_util::VTYPE::V_INS) {
                            if (vv.ref.size() > 1)  {
                                fprintf(stderr, "%d warning: insertion ref allele has size > 1 (%d). Did you left-normalize?\n", vv.pos, vv.ref.size());
                            }
                            if (j != l-1 || ( nop != BAM_CINS && nop != BAM_CSOFT_CLIP )) {
                                // fprintf(stderr, "ins ref allele at %d, %d, %d, %d\n", rpos, j, l, nop);
                                ++vv.ad[0];
                            }
                        } // let other op states handle j == l-1 case
                    }
                }
            }
            r += l;
            a += l;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            // fprintf(stderr, "size %d deletion found at %d\n", l, r);
            auto found = vmap->find(r-1);
            if (found != vmap->end()) {
                for (auto& vv: found->second) {
                    if (vv.type == hts_util::VTYPE::V_DEL) {
                        if (vv.alt.size() > 1)  {
                            fprintf(stderr, "%d warning: deletion alt allele has size > 1 (%d). Did you left-normalize?\n", vv.pos, vv.alt.size());
                        }
                        if (l == vv.ref.size() - vv.alt.size()) {
                            ++vv.ad[1];
                        }
                    }
                }
            }
            r += l;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            // fprintf(stderr, "size %d insertion found at %d\n", l, r);
            auto found = vmap->find(r-1);
            if (found != vmap->end()) {
                for (auto& vv: found->second) {
                    if (vv.type == hts_util::VTYPE::V_INS) {
                        if (vv.ref.size() > 1)  {
                            fprintf(stderr, "%d warning: insertion ref allele has size > 1 (%d). Did you left-normalize?\n", vv.pos, vv.ref.size());
                        }
                        if (l == vv.alt.size() - vv.ref.size()) {
                            ++vv.ad[1];
                        }
                    }
                }
            }
            a += l;
        }
    }
}

void varcount(const VcntArgs& args) {
    samFile* sam_fp = sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();

    vcfFile* vcf_fp = bcf_open(args.vcf_fname.data(), "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    if (bcf_hdr_set_samples(vcf_hdr, NULL, 0)) {
        fprintf(stderr, "error setting samples\n");
        exit(1);
    }

    // TODO: switch to an array implementation
    hts_util::contig2map_map<hts_util::Var> contig2vars(hts_util::bcf_to_map(vcf_fp, vcf_hdr));

    int32_t pid = -1;
    bam1_core_t* c = nullptr;
    hts_util::pos2var_map<hts_util::Var>* vmap = nullptr;
    while (sam_read1(sam_fp, sam_hdr, aln) >= 0) {
        c = &aln->core;
        // we only care about uniquely mapped reads
        if((c->flag & BAM_FUNMAP) == 0 && (c->flag & BAM_FSECONDARY) == 0) {
            if (c->tid != pid) {
                // load new hash table
                vmap = &(contig2vars[sam_hdr->target_name[c->tid]]);
            }
            pid = c->tid;
            aln_varcount(aln, vmap);
        }
    }

    bam_destroy1(aln);
    bam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);


    // output in VCF format!!
    vcfFile* out_vcf_fp = bcf_open("-", "w");
    bcf_hdr_t* out_vcf_hdr = bcf_hdr_dup(vcf_hdr);

    bcf_hdr_add_sample(out_vcf_hdr, args.sample_name.data());
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_INFO, NULL);
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_FMT, NULL);
    // cleanup from original VCF
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, "fileDate");
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, "bcftools_viewVersion");
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_GEN, "bcftools_viewCommand");
    bcf_hdr_remove(out_vcf_hdr, BCF_HL_STR, NULL);
    bcf_hdr_set_version(out_vcf_hdr, "VCFv4.3");
    // bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=ALTCNT,Number=1,Type=Integer,Description=\"Count of reads covering alt allele\">");
    // bcf_hdr_append(out_vcf_hdr, "##INFO=<ID=REFCNT,Number=1,Type=Integer,Description=\"Count of reads covering ref allele\">");
    bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic Depth\">");
    if (args.gt)
        bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if (args.gl)
        bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if (bcf_hdr_write(out_vcf_fp, out_vcf_hdr)) {
        fprintf(stderr, "bcf_hdr_write error\n");
        exit(1);
    }

    bcf1_t* out_vcf_rec = bcf_init();
    bcf_clear(out_vcf_rec);

    // TODO: sort the output?
    for (const auto& vs: contig2vars) { // vs: {seq, map}
        for (const auto& v: vs.second)  { // v: {pos, Varlist}
            auto it = v.second.begin();
            for (auto it = v.second.begin(); it != v.second.end(); ++it) {
                std::string allele_str(it->ref + "," + it->alt);
                out_vcf_rec->rid = bcf_hdr_name2id(out_vcf_hdr, vs.first.data());
                out_vcf_rec->pos = it->pos;
                out_vcf_rec->rlen = it->ref.size();
                bcf_update_alleles_str(out_vcf_hdr, out_vcf_rec, allele_str.data());
                bcf_update_id(out_vcf_hdr, out_vcf_rec, it->id.data());
                bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "ALTCNT", &(it->ad[1]), 1);
                bcf_update_info_int32(out_vcf_hdr, out_vcf_rec, "REFCNT", &(it->ad[0]), 1);
                // std::array<int32_t, 3> pls;
                std::array<double, 3> gls;
                std::array<int32_t, 2> gts;
                if (args.gl || (args.gt == VcntArgs::GT_TYPE::GTL)) {
                    /* TODO: handle indels smarter. Use individual base qualities */
                    // pls = hts_util::get_pls_naive_normalized(it->ad[0], it->ad[1], args.e);
                    gls = hts_util::get_gls_naive(it->ad[0], it->ad[1], args.e);
                }
                int gt_pass = 1;
                int count_pass = (it->ad[0] + it->ad[1] >= args.min_c &&
                                  it->ad[1] >= args.min_ac &&
                                  it->ad[0] >= args.min_rc);
                if (count_pass && args.gt) { // naive genotyping
                    if (args.gt == VcntArgs::GT_TYPE::GTA) {
                        gts = gt_by_alt_evidence(args, *it);
                    } else if (args.gt == VcntArgs::GT_TYPE::GTT) {
                        gts = gt_by_threshold(args, *it);
                    } else { // default: genotype by likelihood
                        // gts = gt_by_likelihood(args, pls);
                        gts = gt_by_likelihood(args, gls);
                    }
                    bcf_update_genotypes(out_vcf_hdr, out_vcf_rec, gts.data(), 2);
                    gt_pass = !bcf_gt_is_missing(gts[0]);
                }
                if (args.keep || ( count_pass  && gt_pass )) {
                    auto pls = hts_util::gl_to_pl(gls);
                    if (args.gl) bcf_update_format_int32(out_vcf_hdr, out_vcf_rec, "PL", pls.data(), 3);
                    bcf_update_format_int32(out_vcf_hdr, out_vcf_rec, "AD", it->ad.data(), 2);
                    if (bcf_write(out_vcf_fp, out_vcf_hdr, out_vcf_rec)) {
                        fprintf(stderr, "bcf_write error\n");
                        exit(1);
                    }
                }
                bcf_clear(out_vcf_rec);
            }
        }
    }
    bcf_hdr_destroy(out_vcf_hdr);
    bcf_close(out_vcf_fp);

    bcf_hdr_destroy(vcf_hdr);
    bcf_close(vcf_fp);
}

void print_help() {
fprintf(stderr,
"Description: \n\
\n\
Given a VCF and SAM file, calculate the alignment coverage over each ALT and\n\
REF allele in the VCF. Outputs in VCF format to stdout.\n\n\
Usage:\n\
\n\
./varcount [options] <vcf> <sam>\n\
\n\
<vcf>=STR               [bv]cf file name (required)\n\
<sam>=STR               [bs]sam file name (required)\n\
-s/--sample-name=STR    sample name in VCF output\n\
                        (default: sample)\n\
-g/--genotype           [likelihood, alt_sensitive, threshold]. 'predict' a genotype in the GT field\n\
                        likelihood: use crude gt likelihood from counts\n\
                        alt_sensitive: automatically call alleles with any alt evidence as alt/alt\n\
                        threshold: use a manual threshold of the difference between alt and ref count to determine gt. Use -c parameter to specifiy threshoold\n\
-c                      int>=0 (default: 0). for use with '-g threshold'. estabilish threshold of ref-alt for determing genotype.  \n\
-a/--min-alt-count      int>=0 (default: 0). filter loci by minimum depth of reads covering alt allele.\n\
-r/--min-ref-count      int>=0 (default: 0). filter loci by minimum depth of reads covering ref allele.\n\
-m/--min-total-count    int>=0 (default: 0). filter loci by minimum depth of total reads.\n\
-k/--keep               ignore filters, print all records regardless of coverage/genotype\n\
-v/--verbose            prints detailed logging information to stderr\n\
-h/--help               print this help message\n\
\n");
}

int main(int argc, char** argv) {
    VcntArgs args;
    static struct option long_options[] {
        {"sample-name", required_argument, 0, 's'},
        {"threshold", required_argument, 0, 'c'},
        {"genotype", required_argument, 0, 'g'},
        {"keep", no_argument, &args.keep, 1},
        {"alt-default", no_argument, &args.alt_default, 1},
        {"min-alt-count", required_argument, 0, 'a'},
        {"min-ref-count", required_argument, 0, 'r'},
        {"min-total-count", required_argument, 0, 'm'},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, &args.verbose, 1},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "-:s:m:a:r:g:c:vkhl", long_options, NULL)) != -1 ) {
        switch(ch) {
            case 0:
                break;
            case 1:
                if (argpos == 0) args.vcf_fname = std::string(optarg);
                else if (argpos == 1) args.sam_fname = std::string(optarg);
                else fprintf(stderr, "ignoring argument %s\n", optarg);
                ++argpos;
                break;
            case 2:
                break;
            case 's':
                args.sample_name = std::string(optarg);
                break;
            case 'c':
                args.thres = std::atoi(optarg);
                break;
            case 'a':
                args.min_ac = std::atoi(optarg);
                break;
            case 'r':
                args.min_rc = std::atoi(optarg);
                break;
            case 'm':
                args.min_c = std::atoi(optarg);
                break;
            case 'g':
                if (!strcmp(optarg, "threshold")) {
                    args.gt = VcntArgs::GT_TYPE::GTT;
                } else if (!strcmp(optarg, "alt_sensitive")) {
                    args.gt = VcntArgs::GT_TYPE::GTA;
                } else if (!strcmp(optarg, "") | !strcmp(optarg, "likelihood")) { // "likelihood"
                    args.gt = VcntArgs::GT_TYPE::GTL;
                } else {
                    args.gt = VcntArgs::GT_TYPE::GTL;
                }
                break;
            case 'k':
                args.keep = 1;
                break;
            case 'v':
                args.verbose = 1;
                break;
            case 'l':
                args.gl = 1;
                break;
            case 'e':
                args.e = std::atof(optarg);
                break;
            case 'h':
                print_help();
                exit(0);
                break;
            case '?':
                print_help();
                fprintf(stderr, "error: unknown option -%c\n", optopt);
                exit(1);
                break;
            default:
                print_help();
                exit(1);
        }
    }

    if (args.vcf_fname == "" || args.sam_fname == "") {
        print_help();
        fprintf(stderr, "error: vcf and sam are mandatory arguments\n");
        exit(1);
    }

    varcount(args);
}
