#include <cstdio>
#include <vector>
#include <string>
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <getopt.h>
#include "hts_util.hpp"

struct VHamArgs {
    std::string fn1 = "";
    std::string fn2 = "";
    std::string s1 = "";
    std::string s2 = "";
    int unphased = 0;
    int flip_gt = 0;
    int score_only = 0;
};

static inline int gta(int i) {
    return bcf_gt_allele(i);
}

struct VarGT {
    VarGT(hts_util::Var v, int g1, int g2) : var(v), gt1(g1), gt2(g2) {}
    int gt1 = 0;
    int gt2 = 0;
    hts_util::Var var;
    
    static std::vector<VarGT> from_bcf(bcf_hdr_t* hdr, bcf1_t* b) {
        std::vector<VarGT> rets;
        int ngt;
        int* gts = hts_util::get_genotype(hdr, b, &ngt);
        if (ngt < 2) { fprintf(stderr, "not diploid!\n"); exit(1); }
        auto vs = hts_util::Var::from_bcf(hdr, b);
        for (int i = 0; i < vs.size(); ++i) {
            /* this hack allows compatibility with vcfs that don't allow multi-allelic lines */
            for (int j = 0; j < 2; ++j) {
                if (gta(gts[j])) { // don't touch if ref allele
                    if (gta(gts[j]) == i+1) { // 1 for the alt allele the indv actually has
                        gts[j] = bcf_gt_is_phased(gts[j]) ? bcf_gt_phased(1) : bcf_gt_unphased(1);
                    } else { // 2 for the alt allele that they don't have
                        gts[j] = bcf_gt_is_phased(gts[j]) ? bcf_gt_phased(2) : bcf_gt_unphased(2);
                    }
                }
            }
            rets.push_back(VarGT(vs[i], gts[0], gts[1]));
        }
        free(gts);
        return rets;
    }
};

static inline bool var_match(const hts_util::Var lv, const hts_util::Var rv) {
    return lv.pos == rv.pos && lv.ref == rv.ref && lv.alt == rv.alt;
}

static inline bool varGT_match(const VarGT lv, const VarGT rv, int gt1, int gt2, int unphased) {
    if (!unphased)  {
        return var_match(lv.var, rv.var) && lv.gt1 == rv.gt1 && lv.gt2 == rv.gt2;
    }
    else {
        int la1 = bcf_gt_allele(lv.gt1);
        int la2 = bcf_gt_allele(lv.gt2);
        int ra1 = bcf_gt_allele(rv.gt1);
        int ra2 = bcf_gt_allele(rv.gt2);
        return var_match(lv.var, rv.var) && ( (la1 == ra1 && la2 == ra2) || (la2 == ra1 && la1 == ra2));
    }
}

void vcf_hamming(const VHamArgs& args) {
    vcfFile* fp1 = bcf_open(args.fn1.data(), "r");
    bcf_hdr_t* hdr1 = bcf_hdr_read(fp1);
    vcfFile* fp2 = bcf_open(args.fn2.data(), "r");
    bcf_hdr_t* hdr2 = bcf_hdr_read(fp2);
    if (args.s1 != "" && bcf_hdr_set_samples(hdr1, args.s1.data(), 0)) {
            fprintf(stderr, "error setting sample for ref vcf file\n");
            exit(1);
    }
    if (args.s2 != "" && bcf_hdr_set_samples(hdr2, args.s2.data(), 0)) {
            fprintf(stderr, "error setting sample for pred vcf file\n");
            exit(1);
    }

    hts_util::contig2map_map<VarGT> contig2vars(hts_util::bcf_to_map<VarGT>(fp1, hdr1));
    bcf1_t* rec = bcf_init();
    hts_util::pos2var_map<VarGT>* vmap = nullptr;
    int pid = -1;
    int hamm = 0;

    vcfFile* out_vcf_fp;
    bcf_hdr_t* out_vcf_hdr;
    if (!args.score_only) {
        out_vcf_fp = bcf_open("-", "w");
        out_vcf_hdr = bcf_hdr_dup(hdr2);
        // bcf_hdr_add_sample(out_vcf_hdr, args.s2 == "" ? "sample" : args.s2.data());
        bcf_hdr_set_version(out_vcf_hdr, "VCFv4.3");
        bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=TG,Number=1,Type=String,Description=\"True genotype\">");
        bcf_hdr_append(out_vcf_hdr, "##FORMAT=<ID=SC,Number=1,Type=Integer,Description=\"1 if 'correct', 0 otherwise\">");
        if (bcf_hdr_write(out_vcf_fp, out_vcf_hdr)) {
            fprintf(stderr, "bcf_hdr_write error\n");
            exit(1);
        }
    }
    while (!bcf_read(fp2, hdr2, rec)) {
        bcf_unpack(rec, BCF_UN_ALL);
        if (rec->rid != pid) {
            vmap = &(contig2vars[bcf_hdr_id2name(hdr2, rec->rid)]);
        }
        pid = rec->rid;
        int match = 0;
        char tgt[3]; 
        const auto found = vmap->find(rec->pos);
        if (found != vmap->end()) {
            const auto& v1s = found->second;
            auto v2s = VarGT::from_bcf(hdr2, rec);
            for (const auto& v2: v2s) {
                for (const auto& v1: v1s) {
                    tgt[0] = '0' + gta(v1.gt1);
                    tgt[1] = '|';
                    tgt[2] = '0' + gta(v1.gt2);
                    if (var_match(v1.var, v2.var)) {
                        int gt_match = 0;
                        if (!args.unphased) {
                            gt_match = args.flip_gt ? (gta(v1.gt1) == gta(v2.gt2) && gta(v1.gt2) == gta(v2.gt1)) : (gta(v1.gt1) == gta(v2.gt1) && gta(v1.gt2) == gta(v2.gt2));
                        } else {
                            gt_match =  gta(v1.gt1) == gta(v2.gt2) && gta(v1.gt2) == gta(v2.gt1);
                            gt_match |= gta(v1.gt1) == gta(v2.gt1) && gta(v1.gt2) == gta(v2.gt2);
                        }
                        if (gt_match) {
                            match = 1;
                            break;
                        }
                    }
                }
                if (!match) { hamm += 1; break;}
            }
        }
        if (match > 1) exit(1);
        if (!args.score_only) {
            bcf_update_format_char(out_vcf_hdr, rec, "TG", &tgt, 3);
            bcf_update_format_int32(out_vcf_hdr, rec, "SC", &match, 1);
            if (bcf_write(out_vcf_fp, out_vcf_hdr, rec)) {
                fprintf(stderr, "error writing record\n");
                exit(1);
            }
        }
    }
    if (args.score_only) {
        fprintf(stdout, "%d\n", hamm);
    } else {
        bcf_hdr_destroy(out_vcf_hdr);
        bcf_close(out_vcf_fp);
        fprintf(stderr, "%d\n", hamm);
    }

    bcf_hdr_destroy(hdr1);
    bcf_hdr_destroy(hdr2);
    bcf_close(fp1);
    bcf_close(fp2);
}

void print_help() {
fprintf(stderr,
        "help message in progress \n");
}

int main(int argc, char** argv) {
    VHamArgs args;
    static struct option long_options[] {
        {"ref_sample", required_argument, 0, 'r'},
        {"pred_sample", required_argument, 0, 'p'},
        {"unphased", no_argument, &args.unphased, 1},
        {"flip-gt", no_argument, &args.flip_gt, 1},
        {"score-only", no_argument, &args.score_only, 1},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "-:r:p:h", long_options, NULL)) != -1 ) {
        switch(ch) {
            case 0:
                break;
            case 1:
                if (argpos == 0) args.fn1 = std::string(optarg);
                else if (argpos == 1) args.fn2 = std::string(optarg);
                else fprintf(stderr, "ignoring argument %s\n", optarg);
                ++argpos;
                break;
            case 'r':
                args.s1 = optarg;
                break;
            case 'p':
                args.s2 = optarg;
                break;
            case 'f':
                args.flip_gt = 1;
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

    if (args.fn1 == "" || args.fn2 == "") {
        print_help();
        fprintf(stderr, "error: vcf and sam are mandatory arguments\n");
        exit(1);
    }

    vcf_hamming(args);
}

