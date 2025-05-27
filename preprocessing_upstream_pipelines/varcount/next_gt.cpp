#include <cstdio>
#include <getopt.h>
#include "hts_util.hpp"

struct CmpVcfArgs {
    std::string fname1;
    std::string fname2;
    std::string sample1;
    std::string sample2;
};

struct VarGT {
    /* Assuming one allele per line /shrug/ */
    VarGT(bcf_hdr_t* hdr, bcf1_t* r) {
        int ngt;
        rec = bcf_dup(r);
        // extract genotype
        int32_t* gts = hts_util::get_genotype(hdr, rec, &ngt);
        if (ngt < 2) { fprintf(stderr, "error getting genotype!\n"); exit(1);}
        gt1 = bcf_gt_allele(gts[0]);
        gt2 = bcf_gt_allele(gts[1]);
        free(gts);
    }
    ~VarGT() {
        bcf_destroy(rec);
    }
    static std::vector<VarGT> from_bcf(bcf_hdr_t* hdr, bcf1_t* rec) {
        std::vector<VarGT> v;
        v.push_back(VarGT(hdr, rec));
        return v;
    }
    int32_t gt1;
    int32_t gt2;
    bcf1_t* rec;
    // the following WON'T be used in comparators!
    bool rec_start = 0;
};

static inline bool var_match(const VarGT& lv, const VarGT& rv) {
    bcf1_t* l = lv.rec;
    bcf1_t* r = rv.rec;
    return l->pos == r->pos \
                   && !strcmp(l->d.allele[1] , r->d.allele[1]) \
                   && !strcmp(l->d.allele[0] , r->d.allele[0]) \
                   && lv.gt1 == rv.gt1 \
                   && lv.gt2 == rv.gt2;
}

void compare_vcfs(CmpVcfArgs args) {
    vcfFile* fp1 = bcf_open(args.fname1.data(), "r");
    bcf_hdr_t* hdr1 = bcf_hdr_read(fp1);
    if (args.sample1 != "" && bcf_hdr_set_samples(hdr1, args.sample1.data(), 0)) return;

    hts_util::contig2map_map<VarGT> contig2vars(hts_util::bcf_to_map<VarGT>(fp1, hdr1));

    vcfFile* fp2 = bcf_open(args.fname2.data(), "r");
    bcf_hdr_t* hdr2 = bcf_hdr_read(fp2);
    if (args.sample2 != "" && bcf_hdr_set_samples(hdr2, args.sample2.data(), 0)) return;
    bcf1_t* rec = bcf_init();

    vcfFile* out_fp = bcf_open("-", "w");
    bcf_hdr_t* out_hdr = bcf_hdr_dup(hdr2);
    // print the header of fp2
    if (!bcf_hdr_write(out_fp, out_hdr)) {
        fprintf(stderr, "error writing header\n"); exit(1);
    }

    int32_t pid = -1;
    hts_util::pos2var_map<VarGT>* vmap = nullptr;
    // We're going to assume one variant per line for both VCF files for now.
    while (!bcf_read(fp2, hdr2, rec)) {
        bcf_unpack(rec, BCF_UN_STR || BCF_UN_FMT);
        auto v2 = VarGT(hdr2, rec);
        if (rec->rid != pid) {
            vmap = &(contig2vars[bcf_hdr_id2name(hdr2, rec->rid)]);
        }
        std::string id(rec->d.id);
        auto vars1_iter = vmap->find(rec->pos);
        bool match = 0;
        if (vars1_iter != vmap->end()) {
            for (const auto& v1: vars1_iter->second) {
                if (var_match(v1, v2)) {
                    match = 1;
                    break;
                }
            }
            if (!match) { // write record here
                if (!bcf_write(out_fp, out_hdr, v2.rec)) {
                    fprintf(stderr, "error writing\n"); exit(1);
                }
            }
        }  // not found
        pid = rec->rid;
    }
   
    bcf_destroy(rec);
    bcf_hdr_destroy(out_hdr);
    bcf_hdr_destroy(hdr1);
    bcf_hdr_destroy(hdr2);
    bcf_close(out_fp);
    bcf_close(fp1);
    bcf_close(fp2);

}

void print_help() {
    fprintf(stderr, "usage: ./compare_vcfs -x sample1 -y sample2 <1.vcf> <2.vcf>...\n\
for all genotyped variants for sample2 in 2.vcf, checks if genotype matches\n\
that for sample 1 in 1.vcf. Ungenotyped variants in both files are ignored.\n\
\n\
output:\n\
matches(int) diffs(int) misses(int)\n");
}

int main(int argc, char** argv) {
    CmpVcfArgs args;
    static struct option long_options[] {
        {"s1", required_argument, NULL, 2},
        {"s2", required_argument, NULL, 3},
        {"help", no_argument, NULL, 'h'},
        {0,0,0,0}
    };

    int ch;
    int argpos = 0;
    while ( (ch = getopt_long(argc, argv, "-:h", long_options, NULL)) != -1 ) {
        switch(ch) {
            case 0:
                break;
            case 1:
                if (argpos == 0) args.fname1 = std::string(optarg);
                else if (argpos == 1) args.fname2 = std::string(optarg);
                else fprintf(stderr, "ignoring argument %s\n", optarg);
                ++argpos;
                break;
            case 2:
                args.sample1 = std::string(optarg);
                break;
            case 3:
                args.sample2 = std::string(optarg);
                break;
            case 'h':
                print_help();
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

    if (args.fname1 == "" || args.fname2 == "") {
        print_help();
        fprintf(stderr, "error: vcf files are mandatory arguments\n");
        exit(1);
    }

    compare_vcfs(args);
}

