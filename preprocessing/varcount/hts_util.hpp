#ifndef VARCOUNT_HPP
#define VARCOUNT_HPP

#include <cstdio>
#include <vector>
#include <array>
#include <list>
#include <string>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <utility>
#include "mdparse.hpp"
#include "flat_hash_map.hpp"
#include "util.hpp"

namespace hts_util {
    enum class VTYPE {V_UNK, V_SNP, V_INS, V_DEL};

    struct Var {
        template<typename Ts>
        Var(int32_t p, VTYPE t, std::string a, std::string r, std::string i, Ts&&...) : pos(p), type(t), alt(a), ref(r), id(i) {}
        Var(int32_t p, VTYPE t, std::string a, std::string r, std::string i) : pos(p), type(t), alt(a), ref(r), id(i) {}
        Var(int32_t p, VTYPE t, std::string a, std::string r) : Var(p,t,a,r,"") {}
        int32_t pos = 0;
        VTYPE type = VTYPE::V_UNK;
        std::string alt = "";
        std::string ref = "";
        // the following WON'T be used in comparators!
        std::string id = "";
        std::array<int32_t, 2> ad = {0, 0};
        // int32_t rc = 0;
        // int32_t ac = 0;
        bool rec_start = 0;

        static std::vector<Var> from_bcf(bcf_hdr_t* hdr, bcf1_t* b) {
            (void) hdr;
            std::vector<Var> vs;
            char* ref = b->d.allele[0];
            // to deal with multi-allele ids:
            std::vector<std::string> ids = parse_ids(b->d.id);
            for (uint32_t i = 1; i < b->n_allele; ++i) {
                char* alt = b->d.allele[i];
                std::string id = ids.size() == (b->n_allele - 1) ? ids[i-1] : std::string(b->d.id);
                if (alt[0] == '.') continue;
                if (strlen(alt)  < strlen(ref)) { // DEL
                    vs.push_back(Var(b->pos, VTYPE::V_DEL, alt, ref, id)); // don't need alt here
                } else if (strlen(alt) > strlen(ref)) { // INS
                    vs.push_back(Var(b->pos, VTYPE::V_INS, alt, ref, id)); // don't need ref here
                } else { // SNP
                    vs.push_back(Var(b->pos, VTYPE::V_SNP, alt, ref, id)); // don't need ref here
                }
            }
            vs[0].rec_start = 1;
            return vs;
        }

        static std::vector<Var> from_bam(bam1_t* aln) {
            std::vector<Var> vs;
            // look at md string and gather dels, snps
            char* md;
            std::list<std::pair<std::string, int32_t>> snps;
            std::list<std::pair<std::string, int32_t>> dels;
            if ((md = bam_aux2Z(bam_aux_get(aln, "MD")))) {
                std::vector<MDPos> mds = md_parse(md);
                for (auto m: mds) {
                    if (m.st == MD_DEL) {
                        std::string ref = "";
                        ref += m.str;
                        dels.push_back(std::pair<std::string, int32_t>(ref, aln->core.pos + m.p - 1));
                    }
                    else if (m.st == MD_SNP) {
                        snps.push_back(std::pair<std::string, int32_t>(m.str, aln->core.pos + m.p));
                    }
                }
            }

            // walk through the CIGAR string
            uint32_t* cs = bam_get_cigar(aln);
            int32_t qpos = 0;
            int32_t rpos = aln->core.pos;
            int32_t rlen, qlen;
            for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
                rlen = (bam_cigar_type(cs[i]) & 2) ? bam_cigar_oplen(cs[i]) : 0;
                qlen = (bam_cigar_type(cs[i]) & 1) ? bam_cigar_oplen(cs[i]) : 0;
                if (bam_cigar_op(cs[i]) == BAM_CINS) {
                    std::string ins_seq = "";
                    ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos-1)]; // adding the previous character here for compatibility with VCF format
                    for (uint32_t j = 0; j < bam_cigar_oplen(cs[i]); ++j) {
                        ins_seq += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos+j)];
                    }
                    vs.push_back(Var(rpos-1, VTYPE::V_INS, ins_seq, ins_seq.substr(0,1)));
                } else if (bam_cigar_op(cs[i]) == BAM_CMATCH) {
                    // check potential snps here
                    for (auto it = snps.begin(); it != snps.end(); ) {
                        int32_t s = it->second;
                        if (s >= rpos && s < rpos + rlen) {
                            std::string snp = "";
                            snp += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos + (s - rpos))];
                            vs.push_back(Var(s, VTYPE::V_SNP, snp, it->first));
                            it = snps.erase(it); // we do this so we don't recheck snps that have already been determined, but... maybe deleting it would actually be more expensive?
                        } else ++it;
                    }
                } else if (bam_cigar_op(cs[i]) == BAM_CDEL) {
                    for (auto it = dels.begin(); it != dels.end(); ) {
                        int32_t d = it->second;
                        if (d == rpos - 1) {
                            std::string alt = "";
                            alt += seq_nt16_str[bam_seqi(bam_get_seq(aln), qpos - 1)];
                            vs.push_back(Var(d, VTYPE::V_DEL, alt, alt + it->first));
                            it = dels.erase(it);
                        } else ++it;
                    }
                }
                rpos += rlen;
                qpos += qlen;
            }
            return vs;
        }

    };





    template <typename T=Var>
    using pos2var_map = ska::flat_hash_map< int32_t, std::vector<T> >;

    template <typename T=Var>
    using contig2map_map = ska::flat_hash_map<std::string, pos2var_map<T>>;


    template <typename T=Var>
    contig2map_map<T> bcf_to_map(vcfFile* vcf_fp, bcf_hdr_t* vcf_hdr) {
        bcf1_t* vcf_rec = bcf_init();
        contig2map_map<T> contig2vars;
        for (int32_t i = 0; i < vcf_hdr->n[BCF_DT_CTG]; ++i) {
            const char* seqk = vcf_hdr->id[BCF_DT_CTG][i].key;
            contig2vars.insert_or_assign(seqk, pos2var_map<T>());
        }

        int32_t pid = -1;
        int32_t ppos = -1;
        pos2var_map<T>* vmap = nullptr;
        std::vector<T> vs;
        while (!bcf_read(vcf_fp, vcf_hdr, vcf_rec)) {
            bcf_unpack(vcf_rec, BCF_UN_STR);
            if (vcf_rec->pos != ppos) {
                if (vmap && vs.size()) {
                    vmap->insert_or_assign(ppos, std::move(vs)); // vs is undefined after this
                }
                if (vcf_rec->rid != pid) { // change hash tables here
                    vmap = &(contig2vars[bcf_hdr_id2name(vcf_hdr, vcf_rec->rid)]);
                }
                vs.clear();
            }
            std::vector<T> more_vs = T::from_bcf(vcf_hdr, vcf_rec);
            vs.insert(vs.end(), more_vs.begin(), more_vs.end());
            ppos = vcf_rec->pos;
            pid = vcf_rec->rid;
        }

        if (vcf_rec->pos == ppos && vmap && vs.size()) {
            vmap->insert_or_assign(ppos, std::move(vs));
        }
        vs.clear();
        bcf_destroy(vcf_rec);
        return contig2vars;
    }


    void print_varlist(std::vector<Var> vs, FILE* out) {
        for (const auto& v: vs) {
            fprintf(out, "(%d %d %s %s", v.pos, static_cast<int>(v.type), v.alt.data(), v.ref.data());
            if (v.id.size()) fprintf(out, " %s", v.id.data());
            if (v.ad[0] + v.ad[1]) fprintf(out, " %d %d", v.ad[0], v.ad[1]);
            fprintf(out, ") ");
        } fprintf(out, "\n");
    }

    int32_t* get_genotype(bcf_hdr_t* hdr, bcf1_t* rec, int* ngt) {
        int32_t *gt_arr = NULL, ngt_arr = 0;
        *ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if ( ngt<=0 ) return NULL; // GT not present
        else return gt_arr;
    }

    /* input: ref allele read count, alt allele read count, uniform error rate
     * output: genotype likelihood for ref/ref, ref/alt, alt/alt genotypes, resp.
     * assumes uniform base qualities. Makes no distinction between SNPS and indels
     */
    inline std::array<int, 3> get_pls_naive(int rc, int ac, double e) {
        std::array<int, 3> pls;
        pls[0] = static_cast<int>( -10 * ((rc*std::log10(1-e)) + (ac*std::log10(e))) );
        pls[1] = static_cast<int>( 10 * (rc + ac) * std::log10(2) );
        pls[2] = static_cast<int>( -10 * ((ac*std::log10(1-e)) + (rc*std::log10(e))) );
        return pls;
    }

    inline std::array<int32_t, 3> get_pls_naive_normalized(int32_t rc, int32_t ac, double e) {
        std::array<int32_t, 3> pls;
        pls[0] = static_cast<int32_t>( -10 * ((rc*std::log10(1-e)) + (ac*std::log10(e))) );
        pls[1] = static_cast<int32_t>( 10 * (rc + ac) * std::log10(2) );
        pls[2] = static_cast<int32_t>( -10 * ((ac*std::log10(1-e)) + (rc*std::log10(e))) );
        int32_t min = pls[0];
        for (int i = 1; i < 3; ++i) {
            if (pls[i] < min) min = pls[i];
        }
        for (int i = 0; i < 3; ++i) {
            pls[i] -= min;
        }
        return pls;
    }

    inline std::array<double, 3> get_gls_naive(int rc, int ac, double e) {
        std::array<double, 3> gls;
        for (int i = 0; i < 3; ++i) {
            gls[2-i] = -(rc+ac) + (rc*std::log2( (2-i)*e + i*(1-e) )) + (ac*std::log2( (2-i)*(1-e) + i*e ));
        }
        return gls;
    }

    inline std::array<int32_t, 3> gl_to_pl(std::array<double, 3> gls) {
        std::array<int32_t, 3> pls;
        for (int i = 0; i < 3; ++i) {
            pls[i] = static_cast<int32_t>(-10 * gls[i]);
        }
        return pls;
    }
};

#endif // VARCOUNT_HPP
