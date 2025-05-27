#include <cstdio>
#include <deque>
#include <map>
#include <cmath>
#include <string>
extern "C" {
#include <htslib/vcf.h>
}

int fix_allele(int x) {
    if (bcf_gt_is_missing(x)) {
        return 0;
    } else {
        return bcf_gt_allele(x);
    }
}

struct Var {
        Var() {}
        Var(int32_t r, hts_pos_t p, int t1, int t2, int q1, int q2)
            : rid(r)
            , pos(p)
        {
            th1 = fix_allele(t1);
            th2 = fix_allele(t2);
            qh1 = fix_allele(q1);
            qh2 = fix_allele(q2);
        }
        int32_t rid;
        hts_pos_t pos;
        int th1;
        int th2;
        int qh1;
        int qh2;
};

int window_score(const std::deque<Var>& window) {
    int score = 0;
    int s11 = 1, s12 = 1;
    int s21 = 1, s22 = 1;
    for (auto it = window.begin(); it != window.end(); ++it) {
        if (it->th1 != it->qh1) {
            s11 = 0;
        }
        if (it->th1 != it->qh2) {
            s12 = 0;
        }
        if (it->th2 != it->qh1) {
            s22 = 0;
        }
        if (it->th2 != it->qh2) {
            s22 = 0;
        }
    }
    if (s11 || s12) {
        score += 1;
    }
    if (s21 || s22) {
        score += 1;
    }
    return score;
}

struct Score {
    int tps = 0;
    int fps = 0;
    int fns = 0;
    int tns = 0;
    int fhoms = 0;
    int fhets = 0;
    int thoms = 0;
    int thets = 0;
    Score& operator+=(const Score rhs) {
        tps += rhs.tps;
        tns += rhs.tns;
        fps += rhs.fps;
        fns += rhs.fns;
        fhoms += rhs.fhoms;
        fhets += rhs.fhets;
        thoms += rhs.thoms;
        thets += rhs.thets;
        return *this;
    }
};




Score variant_score(const Var v) {
    Score s;
    // allele identity
    if (v.qh1 == v.th1 || v.qh1 == v.th2) {
        if (v.qh1 == 0) s.tns += 1;
        else s.tps += 1;
    } else {
        if (v.qh1 == 0) s.fns += 1;
        else s.fps += 1;
    }
    if (v.qh2 == v.th1 || v.qh2 == v.th2) {
        if (v.qh2 == 0) s.tns += 1;
        else s.tps += 1;
    } else {
        if (v.qh2 == 0) s.fns += 1;
        else s.fps += 1;
    }
    // hom/het
    if (v.qh1 == v.qh2) {
        if (v.th1 == v.th2) s.thoms = 1;
        else s.fhoms  = 1;
    }
    else if (v.qh1 != v.qh2) {
        if (v.th1 != v.th2) s.thets = 1;
        else s.fhets = 1;
    }
    return s;
}

struct WindowScores {
    void update(int i, size_t nvars) {
        if (nvars < 6) {
            win_1_5 += i;
            total_win_1_5 += 2;
        } else if (nvars < 11) {
            win_6_10 += i;
            total_win_6_10 += 2;
        } else {
            win_11_n += i;
            total_win_11_n += 2;
        }
    }
    int win_1_5 = 0;
    int win_6_10 = 0;
    int win_11_n = 0;
    int total_win_1_5 = 0;
    int total_win_6_10 = 0;
    int total_win_11_n = 0;
};


enum VARTYPE {ANY, SNP, INDEL, SV};
VARTYPE get_variant_type(bcf1_t* rec) {
    int len = std::abs(int(strlen(rec->d.allele[0])) - int(strlen(rec->d.allele[1])));
    if (len == 0) {
        return SNP;
    } else if (len <  50) {
        return INDEL;
    } else return SV;
}

std::string vartype_to_string(VARTYPE vtype) {
    switch(vtype) {
        case ANY:
            return "any";
        case  SNP:
            return "snp";
        case INDEL:
            return "indel";
        case SV:
            return "sv";
        default:
            return "any";
    }
}


void print_score(const Score& score, VARTYPE vtype, size_t total, int wsize, const WindowScores& win_score) {
    fprintf(stdout, "{");
    fprintf(stdout, " \"var_type\": \"%s\",", vartype_to_string(vtype).data());
    fprintf(stdout, " \"total_vars\": %lu,", total);
    fprintf(stdout, " \"fps\": %d,", score.fps);
    fprintf(stdout, " \"fns\": %d,", score.fns);
    fprintf(stdout, " \"tps\": %d,", score.tps);
    fprintf(stdout, " \"tns\": %d,", score.tns);
    fprintf(stdout, " \"fhoms\": %d,", score.fhoms);
    fprintf(stdout, " \"fhets\": %d,", score.fhets);
    fprintf(stdout, " \"thoms\": %d,", score.thoms);
    fprintf(stdout, " \"thets\": %d,", score.thets);
    if (score.tps + score.fps)
        fprintf(stdout, " \"alt_precision\": %.5f,", float(score.tps) / float(score.tps + score.fps));
    else
        fprintf(stdout, " \"alt_precision\": null,");
    if (score.tps + score.fns)
        fprintf(stdout, " \"alt_recall\": %.5f,", float(score.tps) / float(score.tps + score.fns));
    else
        fprintf(stdout, " \"alt_recall\": null,");
    if (score.thets + score.fhets)
        fprintf(stdout, " \"het_precision\": %.5f,", float(score.thets) / float(score.thets + score.fhets));
    else
        fprintf(stdout, " \"het_precision\": null,");
    if  (score.thets + score.fhoms)
        fprintf(stdout, " \"het_recall\": %.5f,", float(score.thets) / float(score.thets + score.fhoms));
    else
        fprintf(stdout, " \"het_recall\": null,");
    if (vtype == ANY) {
        fprintf(stdout, " \"matching_%dbp_windows_1_5\": %d,", wsize, win_score.win_1_5);
        fprintf(stdout, " \"matching_%dbp_windows_6_10\": %d,", wsize, win_score.win_6_10);
        fprintf(stdout, " \"matching_%dbp_windows_11_n\": %d,", wsize, win_score.win_11_n);
        fprintf(stdout, " \"total_%dbp_windows_1_5\": %d,", wsize, win_score.total_win_1_5);
        fprintf(stdout, " \"total_%dbp_windows_6_10\": %d,", wsize, win_score.total_win_6_10);
        fprintf(stdout, " \"total_%dbp_windows_11_n\": %d,", wsize, win_score.total_win_11_n);
        fprintf(stdout, " \"perc_matching_%dbp_windows_1_5\": %.2f,", wsize, static_cast<float>(win_score.win_1_5)/win_score.total_win_1_5);
        fprintf(stdout, " \"perc_matching_%dbp_windows_6_10\": %.2f,", wsize, static_cast<float>(win_score.win_6_10)/win_score.total_win_6_10);
        fprintf(stdout, " \"perc_matching_%dbp_windows_11_n\": %.2f", wsize, static_cast<float>(win_score.win_11_n)/win_score.total_win_11_n);
    } else {
        fprintf(stdout, " \"matching_%dbp_windows_1_5\": null,", wsize);
        fprintf(stdout, " \"matching_%dbp_windows_6_10\": null,", wsize);
        fprintf(stdout, " \"matching_%dbp_windows_11_n\": null,", wsize);
        fprintf(stdout, " \"total_%dbp_windows_1_5\": null,", wsize);
        fprintf(stdout, " \"total_%dbp_windows_6_10\": null,", wsize);
        fprintf(stdout, " \"total_%dbp_windows_11_n\": null,", wsize);
        fprintf(stdout, " \"perc_matching_%dbp_windows_1_5\": null,", wsize);
        fprintf(stdout, " \"perc_matching_%dbp_windows_6_10\": null,", wsize);
        fprintf(stdout, " \"perc_matching_%dbp_windows_11_n\": null", wsize);
    }
    fprintf(stdout, " }\n");
}


int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "please specify a VCF file\n");
        exit(1);
    }
    int wsize = 200;
    if (argc >= 3) {
        wsize = std::atoi(argv[2]);
    }
    FILE* ofp = NULL;
    if (argc >= 4) {
        ofp = fopen(argv[3], "w");
    }
    vcfFile *fp = bcf_open(argv[1], "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    bcf1_t* rec = bcf_init();
    int err = 0;
    int ngt = bcf_hdr_nsamples(hdr);
    int32_t *gt_arr = NULL, ngt_arr = 0;
    std::deque<Var> window;
    WindowScores win_score;
    std::map<VARTYPE, Score> scores_by_variant;
    scores_by_variant[ANY] = Score();
    scores_by_variant[SNP] = Score();
    scores_by_variant[INDEL] = Score();
    scores_by_variant[SV] = Score();
    std::map<VARTYPE, size_t> totals;
    totals[ANY] = 0;
    totals[SNP] = 0;
    totals[INDEL] = 0;
    totals[SV] =  0;
    while (!(err = bcf_read(fp, hdr, rec))) {
        bcf_unpack(rec, BCF_UN_FMT | BCF_UN_STR);
        // record window when we change contig or reach a certain size
        if (window.size() && (rec->rid != window.front().rid || rec->pos - window.front().pos > wsize)) {
            win_score.update(window_score(window), window.size());
            if (rec->rid != window.front().rid) {
                window.clear();
            }
            // TODO: instead of using the difference, count up base tota
            while (window.size() && rec->pos - window.front().pos > wsize) {
                window.pop_front();
            }
        }
        ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
        if ( ngt<=0 ) continue; // GT not present
        if (ngt != 4) {
            fprintf(stderr, "bad line!\n");
            exit(1);
        }
        if (ofp) {
            //fprintf(ofp, "%s %d %ld %d/%d %d/%d\n", bcf_seqname(hdr, rec), rec->rid, rec->pos, fix_allele(gt_arr[0]), fix_allele(gt_arr[1]), fix_allele(gt_arr[2]), fix_allele(gt_arr[3]));
            fprintf(ofp, "%s %d %lld %d/%d %d/%d\n", bcf_seqname(hdr, rec), rec->rid, rec->pos, fix_allele(gt_arr[0]), fix_allele(gt_arr[1]), fix_allele(gt_arr[2]), fix_allele(gt_arr[3]));
        }
        auto var = Var(rec->rid, rec->pos, gt_arr[0], gt_arr[1], gt_arr[2], gt_arr[3]);
        window.push_back(var);
        Score s = variant_score(var);
        scores_by_variant[ANY] += s;
        totals[ANY] += 1;
        VARTYPE vtype = get_variant_type(rec);
        scores_by_variant[vtype] += s;
        totals[vtype] += 2;
    }
    if (window.size()) {
        win_score.update(window_score(window), window.size());
    }
    fprintf(stdout, "[ ");
    print_score(scores_by_variant[ANY], ANY, totals[ANY], wsize, win_score);
    fprintf(stdout, ", ");
    print_score(scores_by_variant[SNP], SNP, totals[SNP], wsize, win_score);
    fprintf(stdout, ", ");
    print_score(scores_by_variant[INDEL], INDEL, totals[INDEL], wsize, win_score);
    fprintf(stdout, ", ");
    print_score(scores_by_variant[SV], SV, totals[SV], wsize, win_score);
    fprintf(stdout, " ]\n");
    if (ofp != NULL) {
        fclose(ofp);
    }
    return 0;
}
