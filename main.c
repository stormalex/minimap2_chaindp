#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"

#include "fpga_chaindp.h"

#ifdef HAVE_GETOPT
#include <getopt.h>
#else
#include "getopt.h"
#endif

#define MM_VERSION "2.10-r761"

#ifdef __linux__
#include <sys/resource.h>
#include <sys/time.h>
void liftrlimit()
{
	struct rlimit r;
	getrlimit(RLIMIT_AS, &r);
	r.rlim_cur = r.rlim_max;
	setrlimit(RLIMIT_AS, &r);
}
#else
void liftrlimit() {}
#endif

#include <sys/time.h>
#include<time.h>
static double realtime_msec(void)
{
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return tp.tv_sec*1000 + tp.tv_nsec*1e-6;
}

static struct option long_options[] = {
	{ "bucket-bits",    required_argument, 0, 0 },
	{ "mb-size",        required_argument, 0, 'K' },
	{ "seed",           required_argument, 0, 0 },
	{ "no-kalloc",      no_argument,       0, 0 },
	{ "print-qname",    no_argument,       0, 0 },
	{ "no-self",        no_argument,       0, 'D' },
	{ "print-seeds",    no_argument,       0, 0 },
	{ "max-chain-skip", required_argument, 0, 0 },
	{ "min-dp-len",     required_argument, 0, 0 },
	{ "print-aln-seq",  no_argument,       0, 0 },
	{ "splice",         no_argument,       0, 0 },
	{ "cost-non-gt-ag", required_argument, 0, 'C' },
	{ "no-long-join",   no_argument,       0, 0 },
	{ "sr",             no_argument,       0, 0 },
	{ "frag",           required_argument, 0, 0 },
	{ "secondary",      required_argument, 0, 0 },
	{ "cs",             optional_argument, 0, 0 },
	{ "end-bonus",      required_argument, 0, 0 },
	{ "no-pairing",     no_argument,       0, 0 },
	{ "splice-flank",   required_argument, 0, 0 },
	{ "idx-no-seq",     no_argument,       0, 0 },
	{ "end-seed-pen",   required_argument, 0, 0 },   // 21
	{ "for-only",       no_argument,       0, 0 },   // 22
	{ "rev-only",       no_argument,       0, 0 },   // 23
	{ "heap-sort",      required_argument, 0, 0 },   // 24
	{ "all-chain",      no_argument,       0, 'P' },
	{ "dual",           required_argument, 0, 0 },   // 26
	{ "max-clip-ratio", required_argument, 0, 0 },   // 27
	{ "min-occ-floor",  required_argument, 0, 0 },   // 28
	{ "MD",             no_argument,       0, 0 },   // 29
	{ "help",           no_argument,       0, 'h' },
	{ "max-intron-len", required_argument, 0, 'G' },
	{ "version",        no_argument,       0, 'V' },
	{ "min-count",      required_argument, 0, 'n' },
	{ "min-chain-score",required_argument, 0, 'm' },
	{ "mask-level",     required_argument, 0, 'M' },
	{ "min-dp-score",   required_argument, 0, 's' },
	{ "sam",            no_argument,       0, 'a' },
	{ 0, 0, 0, 0}
};

static inline int64_t mm_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(optarg, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

static inline void yes_or_no(mm_mapopt_t *opt, int flag, int long_idx, const char *arg, int yes_to_set)
{
	if (yes_to_set) {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag |= flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag &= ~flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	} else {
		if (strcmp(arg, "yes") == 0 || strcmp(arg, "y") == 0) opt->flag &= ~flag;
		else if (strcmp(arg, "no") == 0 || strcmp(arg, "n") == 0) opt->flag |= flag;
		else fprintf(stderr, "[WARNING]\033[1;31m option '--%s' only accepts 'yes' or 'no'.\033[0m\n", long_options[long_idx].name);
	}
}

int max_task = 0;

double result_time[100];
double send_task1[100];
double send_task2[100];
double process_result[100];

double soft_chaindp_time[100];
int soft_chaindp_num = 0;

typedef struct _params{
    mm_idx_reader_t* idx_rdr;
    mm_idx_t* mi[2];
    mm_mapopt_t* opt;
    int n_threads;
    pthread_mutex_t mutex;
    pthread_cond_t cond_r;
    pthread_cond_t cond_w;
    int counter;
    int exit;
    
    char *rg;
    int argc;
    char **argv;
}params_t;

void *read_task_thread(void* args)
{
    int mi_index = 0;
    params_t* params = (params_t*)args;
    mm_idx_reader_t *idx_rdr = params->idx_rdr;
    int n_threads = params->n_threads;
    
    while(params->exit) {
        pthread_mutex_lock(&params->mutex);
        mm_idx_t* mi = params->mi[mi_index];
        if(mi == NULL) {
            mi = mm_idx_reader_read(idx_rdr, n_threads);
            params->mi[mi_index] = mi;
            pthread_cond_signal(&params->cond_r);
            pthread_mutex_unlock(&params->mutex);
            mi_index++;
            if(mi_index == 2) {
                mi_index = 0;
            }
        }
        else {
            fprintf(stderr, "wait idle mi\n");
            pthread_cond_wait(&params->cond_w, &params->mutex);
            pthread_mutex_unlock(&params->mutex);
            continue;
        }
    }
    
    
    
    return NULL;
}
void load_index_buf(idx_buf_t* idx_buf, int type);
void idx_buf_destroy(idx_buf_t* idx_buf);
void *map_task_thread(void* args)
{
    int i = 0;
    int mi_index = 0;
    params_t* params = (params_t*)args;
    mm_idx_reader_t *idx_rdr = params->idx_rdr;
    mm_mapopt_t* opt = params->opt;
    char *rg = params->rg;
    char **argv = params->argv;
    int argc = params->argc;
    int n_threads = params->n_threads;
    
    while(params->exit) {
        pthread_mutex_lock(&params->mutex);
        mm_idx_t* mi = params->mi[mi_index];
        pthread_mutex_unlock(&params->mutex);
        if(mi != NULL) {
            
            //load index
            load_index_buf(mi->b_idx, TYPE_INDEX_B);
            load_index_buf(mi->h_idx, TYPE_INDEX_H);
            load_index_buf(mi->v_idx, TYPE_INDEX_V);
            load_index_buf(mi->p_idx, TYPE_INDEX_P);
            idx_buf_destroy(mi->b_idx);
            idx_buf_destroy(mi->h_idx);
            idx_buf_destroy(mi->v_idx);
            idx_buf_destroy(mi->p_idx);
            mi->b_idx = NULL;
            mi->h_idx = NULL;
            mi->v_idx = NULL;
            mi->p_idx = NULL;
            
            if ((opt->flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
                fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
                mm_idx_destroy(mi);
                mm_idx_reader_close(idx_rdr);
                return NULL;
            }
            
            if ((opt->flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
                if (mm_idx_reader_eof(idx_rdr)) {
                    mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
                } else {
                    mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
                    if (mm_verbose >= 2)
                        fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted.\033[0m\n");
                }
            }
            
            if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
            if (argc != optind + 1) mm_mapopt_update(opt, mi);

            if (mm_verbose >= 3) mm_idx_stat(mi);

            //bw,  is_cdna,  max_skip,  min_sc,  flag,  max_occ
            fprintf(stderr, "fpga params:bw:0x%x, is_cdna:0x%x, max_skip:0x%x, min_sc:0x%x, flag:0x%x, max_occ:0x%x\n", opt->bw, !!(opt->flag & MM_F_SPLICE), opt->max_chain_skip, opt->min_chain_score, opt->flag, opt->mid_occ);
            fpga_set_params(opt->bw, !!(opt->flag & MM_F_SPLICE), opt->max_chain_skip, opt->min_chain_score, opt->flag, opt->mid_occ);  //data1

            if (!(opt->flag & MM_F_FRAG_MODE)) {
                for (i = optind + 1; i < argc; ++i)
                    mm_map_file(mi, argv[i], opt, n_threads);
            } else {
                mm_map_file_frag(mi, argc - (optind + 1), (const char**)&argv[optind + 1], opt, n_threads);
            }

            mm_idx_destroy(mi);
            pthread_mutex_lock(&params->mutex);
            params->mi[mi_index] = 0;
            pthread_mutex_unlock(&params->mutex);
            pthread_cond_signal(&params->cond_w);       //唤醒写线程
            mi_index++;
            if(mi_index == 2) {
                mi_index = 0;
            }
        }
        else {
            fprintf(stderr, "wait usable mi\n");
            pthread_cond_wait(&params->cond_r, &params->mutex);
            pthread_mutex_unlock(&params->mutex);
        }
    }
    return NULL;
}

int main(int argc, char *argv[])
{
	const char *opt_str = "2aSDw:k:K:t:r:f:Vv:g:G:I:d:XT:s:x:Hcp:M:n:z:A:B:O:E:m:N:Qu:R:hF:LC:y";
	mm_mapopt_t opt;
	mm_idxopt_t ipt;
	int i, c, n_threads = 3, long_idx;
	char *fnw = 0, *rg = 0, *s;
	FILE *fp_help = stderr;
	mm_idx_reader_t *idx_rdr;
	mm_idx_t *mi;
    double t1, t2;
    t1 = realtime_msec();
    
    memset(result_time, 0, sizeof(result_time));
    memset(send_task1, 0, sizeof(send_task1));
    memset(send_task2, 0, sizeof(send_task2));
    memset(process_result, 0, sizeof(process_result));
    memset(soft_chaindp_time, 0, sizeof(soft_chaindp_time));
    
    fprintf(stderr, "sizeof(struct new_seed)=%ld\n", sizeof(struct new_seed));
    fprintf(stderr, "sizeof(collect_task_t)=%ld\n", sizeof(collect_task_t));
    assert(sizeof(collect_task_t) == 64);
    fprintf(stderr, "sizeof(chaindp_sndhdr_t)=%ld\n", sizeof(chaindp_sndhdr_t));
    assert(sizeof(chaindp_sndhdr_t) == 64);
    fprintf(stderr, "sizeof(collect_result_t)=%ld\n", sizeof(collect_result_t));
    assert(sizeof(collect_result_t) == 64);
    
	mm_verbose = 3;
	liftrlimit();
	mm_realtime0 = realtime();
	mm_set_opt(0, &ipt, &opt);

	while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) >= 0) // apply option -x/preset first
		if (c == 'x') {
			if (mm_set_opt(optarg, &ipt, &opt) < 0) {
				fprintf(stderr, "[ERROR] unknown preset '%s'\n", optarg);
				return 1;
			}
			break;
		}
	optind = 0; // for musl getopt, optind=0 has the same effect as optreset=1; older libc doesn't have optreset

	while ((c = getopt_long(argc, argv, opt_str, long_options, &long_idx)) >= 0) {
		if (c == 'w') ipt.w = atoi(optarg);
		else if (c == 'k') ipt.k = atoi(optarg);
		else if (c == 'H') ipt.flag |= MM_I_HPC;
		else if (c == 'd') fnw = optarg; // the above are indexing related options, except -I
		else if (c == 'r') opt.bw = (int)mm_parse_num(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') mm_verbose = atoi(optarg);
		else if (c == 'g') opt.max_gap = (int)mm_parse_num(optarg);
		else if (c == 'G') mm_mapopt_max_intron_len(&opt, (int)mm_parse_num(optarg));
		else if (c == 'F') opt.max_frag_len = (int)mm_parse_num(optarg);
		else if (c == 'N') opt.best_n = atoi(optarg);
		else if (c == 'p') opt.pri_ratio = atof(optarg);
		else if (c == 'M') opt.mask_level = atof(optarg);
		else if (c == 'c') opt.flag |= MM_F_OUT_CG | MM_F_CIGAR;
		else if (c == 'D') opt.flag |= MM_F_NO_DIAG;
		else if (c == 'P') opt.flag |= MM_F_ALL_CHAINS;
		else if (c == 'X') opt.flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN; // -D -P --no-long-join --dual=no
		else if (c == 'a') opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
		else if (c == 'Q') opt.flag |= MM_F_NO_QUAL;
		else if (c == 'Y') opt.flag |= MM_F_SOFTCLIP;
		else if (c == 'L') opt.flag |= MM_F_LONG_CIGAR;
		else if (c == 'y') opt.flag |= MM_F_COPY_COMMENT;
		else if (c == 'T') opt.sdust_thres = atoi(optarg);
		else if (c == 'n') opt.min_cnt = atoi(optarg);
		else if (c == 'm') opt.min_chain_score = atoi(optarg);
		else if (c == 'A') opt.a = atoi(optarg);
		else if (c == 'B') opt.b = atoi(optarg);
		else if (c == 's') opt.min_dp_max = atoi(optarg);
		else if (c == 'C') opt.noncan = atoi(optarg);
		else if (c == 'I') ipt.batch_size = mm_parse_num(optarg);
		else if (c == 'K') opt.mini_batch_size = (int)mm_parse_num(optarg);
		else if (c == 'R') rg = optarg;
		else if (c == 'h') fp_help = stdout;
		else if (c == '2') opt.flag |= MM_F_2_IO_THREADS;
		else if (c == 0 && long_idx == 0) ipt.bucket_bits = atoi(optarg); // --bucket-bits
		else if (c == 0 && long_idx == 2) opt.seed = atoi(optarg); // --seed
		else if (c == 0 && long_idx == 3) mm_dbg_flag |= MM_DBG_NO_KALLOC; // --no-kalloc
		else if (c == 0 && long_idx == 4) mm_dbg_flag |= MM_DBG_PRINT_QNAME; // --print-qname
		else if (c == 0 && long_idx == 6) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_SEED, n_threads = 1; // --print-seed
		else if (c == 0 && long_idx == 7) opt.max_chain_skip = atoi(optarg); // --max-chain-skip
		else if (c == 0 && long_idx == 8) opt.min_ksw_len = atoi(optarg); // --min-dp-len
		else if (c == 0 && long_idx == 9) mm_dbg_flag |= MM_DBG_PRINT_QNAME | MM_DBG_PRINT_ALN_SEQ, n_threads = 1; // --print-aln-seq
		else if (c == 0 && long_idx ==10) opt.flag |= MM_F_SPLICE; // --splice
		else if (c == 0 && long_idx ==12) opt.flag |= MM_F_NO_LJOIN; // --no-long-join
		else if (c == 0 && long_idx ==13) opt.flag |= MM_F_SR; // --sr
		else if (c == 0 && long_idx ==17) opt.end_bonus = atoi(optarg); // --end-bonus
		else if (c == 0 && long_idx ==18) opt.flag |= MM_F_INDEPEND_SEG; // --no-pairing
		else if (c == 0 && long_idx ==20) ipt.flag |= MM_I_NO_SEQ; // --idx-no-seq
		else if (c == 0 && long_idx ==21) opt.anchor_ext_shift = atoi(optarg); // --end-seed-pen
		else if (c == 0 && long_idx ==22) opt.flag |= MM_F_FOR_ONLY; // --for-only
		else if (c == 0 && long_idx ==23) opt.flag |= MM_F_REV_ONLY; // --rev-only
		else if (c == 0 && long_idx ==27) opt.max_clip_ratio = atof(optarg); // --max-clip-ratio
		else if (c == 0 && long_idx ==28) opt.min_mid_occ = atoi(optarg); // --min-occ-floor
		else if (c == 0 && long_idx ==29) opt.flag |= MM_F_OUT_MD; // --MD
		else if (c == 0 && long_idx == 14) { // --frag
			yes_or_no(&opt, MM_F_FRAG_MODE, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 15) { // --secondary
			yes_or_no(&opt, MM_F_NO_PRINT_2ND, long_idx, optarg, 0);
		} else if (c == 0 && long_idx == 16) { // --cs
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR;
			if (optarg == 0 || strcmp(optarg, "short") == 0) {
				opt.flag &= ~MM_F_OUT_CS_LONG;
			} else if (strcmp(optarg, "long") == 0) {
				opt.flag |= MM_F_OUT_CS_LONG;
			} else if (strcmp(optarg, "none") == 0) {
				opt.flag &= ~MM_F_OUT_CS;
			} else if (mm_verbose >= 2) {
				fprintf(stderr, "[WARNING]\033[1;31m --cs only takes 'short' or 'long'. Invalid values are assumed to be 'short'.\033[0m\n");
			}
		} else if (c == 0 && long_idx == 19) { // --splice-flank
			yes_or_no(&opt, MM_F_SPLICE_FLANK, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 24) { // --heap-sort
			yes_or_no(&opt, MM_F_HEAP_SORT, long_idx, optarg, 1);
		} else if (c == 0 && long_idx == 26) { // --dual
			yes_or_no(&opt, MM_F_NO_DUAL, long_idx, optarg, 0);
		} else if (c == 'S') {
			opt.flag |= MM_F_OUT_CS | MM_F_CIGAR | MM_F_OUT_CS_LONG;
			if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING]\033[1;31m option -S is deprecated and may be removed in future. Please use --cs=long instead.\033[0m\n");
		} else if (c == 'V') {
			puts(MM_VERSION);
			return 0;
		} else if (c == 'f') {
			double x;
			char *p;
			x = strtod(optarg, &p);
			if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
			else opt.mid_occ = (int)(x + .499);
			if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
		} else if (c == 'u') {
			if (*optarg == 'b') opt.flag |= MM_F_SPLICE_FOR|MM_F_SPLICE_REV; // both strands
			else if (*optarg == 'f') opt.flag |= MM_F_SPLICE_FOR, opt.flag &= ~MM_F_SPLICE_REV; // match GT-AG
			else if (*optarg == 'r') opt.flag |= MM_F_SPLICE_REV, opt.flag &= ~MM_F_SPLICE_FOR; // match CT-AC (reverse complement of GT-AG)
			else if (*optarg == 'n') opt.flag &= ~(MM_F_SPLICE_FOR|MM_F_SPLICE_REV); // don't try to match the GT-AG signal
			else {
				fprintf(stderr, "[ERROR]\033[1;31m unrecognized cDNA direction\033[0m\n");
				return 1;
			}
		} else if (c == 'z') {
			opt.zdrop = opt.zdrop_inv = strtol(optarg, &s, 10);
			if (*s == ',') opt.zdrop_inv = strtol(s + 1, &s, 10);
		} else if (c == 'O') {
			opt.q = opt.q2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.q2 = strtol(s + 1, &s, 10);
		} else if (c == 'E') {
			opt.e = opt.e2 = strtol(optarg, &s, 10);
			if (*s == ',') opt.e2 = strtol(s + 1, &s, 10);
		}
	}
	if ((opt.flag & MM_F_SPLICE) && (opt.flag & MM_F_FRAG_MODE)) {
		fprintf(stderr, "[ERROR]\033[1;31m --splice and --frag should not be specified at the same time.\033[0m\n");
		return 1;
	}
	if (!fnw && !(opt.flag&MM_F_CIGAR))
		ipt.flag |= MM_I_NO_SEQ;
	if (mm_check_opt(&ipt, &opt) < 0)
		return 1;

	if (argc == optind || fp_help == stdout) {
		fprintf(fp_help, "Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]\n");
		fprintf(fp_help, "Options:\n");
		fprintf(fp_help, "  Indexing:\n");
		fprintf(fp_help, "    -H           use homopolymer-compressed k-mer\n");
		fprintf(fp_help, "    -k INT       k-mer size (no larger than 28) [%d]\n", ipt.k);
		fprintf(fp_help, "    -w INT       minizer window size [%d]\n", ipt.w);
		fprintf(fp_help, "    -I NUM       split index for every ~NUM input bases [4G]\n");
		fprintf(fp_help, "    -d FILE      dump index to FILE []\n");
		fprintf(fp_help, "  Mapping:\n");
		fprintf(fp_help, "    -f FLOAT     filter out top FLOAT fraction of repetitive minimizers [%g]\n", opt.mid_occ_frac);
		fprintf(fp_help, "    -g NUM       stop chain enlongation if there are no minimizers in INT-bp [%d]\n", opt.max_gap);
		fprintf(fp_help, "    -G NUM       max intron length (effective with -xsplice; changing -r) [200k]\n");
		fprintf(fp_help, "    -F NUM       max fragment length (effective with -xsr or in the fragment mode) [800]\n");
		fprintf(fp_help, "    -r NUM       bandwidth used in chaining and DP-based alignment [%d]\n", opt.bw);
		fprintf(fp_help, "    -n INT       minimal number of minimizers on a chain [%d]\n", opt.min_cnt);
		fprintf(fp_help, "    -m INT       minimal chaining score (matching bases minus log gap penalty) [%d]\n", opt.min_chain_score);
//		fprintf(fp_help, "    -T INT       SDUST threshold; 0 to disable SDUST [%d]\n", opt.sdust_thres); // TODO: this option is never used; might be buggy
		fprintf(fp_help, "    -X           skip self and dual mappings (for the all-vs-all mode)\n");
		fprintf(fp_help, "    -p FLOAT     min secondary-to-primary score ratio [%g]\n", opt.pri_ratio);
		fprintf(fp_help, "    -N INT       retain at most INT secondary alignments [%d]\n", opt.best_n);
		fprintf(fp_help, "  Alignment:\n");
		fprintf(fp_help, "    -A INT       matching score [%d]\n", opt.a);
		fprintf(fp_help, "    -B INT       mismatch penalty [%d]\n", opt.b);
		fprintf(fp_help, "    -O INT[,INT] gap open penalty [%d,%d]\n", opt.q, opt.q2);
		fprintf(fp_help, "    -E INT[,INT] gap extension penalty; a k-long gap costs min{O1+k*E1,O2+k*E2} [%d,%d]\n", opt.e, opt.e2);
		fprintf(fp_help, "    -z INT[,INT] Z-drop score and inversion Z-drop score [%d,%d]\n", opt.zdrop, opt.zdrop_inv);
		fprintf(fp_help, "    -s INT       minimal peak DP alignment score [%d]\n", opt.min_dp_max);
		fprintf(fp_help, "    -u CHAR      how to find GT-AG. f:transcript strand, b:both strands, n:don't match GT-AG [n]\n");
		fprintf(fp_help, "  Input/Output:\n");
		fprintf(fp_help, "    -a           output in the SAM format (PAF by default)\n");
		fprintf(fp_help, "    -Q           don't output base quality in SAM\n");
		fprintf(fp_help, "    -L           write CIGAR with >65535 ops at the CG tag\n");
		fprintf(fp_help, "    -R STR       SAM read group line in a format like '@RG\\tID:foo\\tSM:bar' []\n");
		fprintf(fp_help, "    -c           output CIGAR in PAF\n");
		fprintf(fp_help, "    --cs[=STR]   output the cs tag; STR is 'short' (if absent) or 'long' [none]\n");
		fprintf(fp_help, "    --MD         output the MD tag\n");
		fprintf(fp_help, "    -Y           use soft clipping for supplementary alignments\n");
		fprintf(fp_help, "    -t INT       number of threads [%d]\n", n_threads);
		fprintf(fp_help, "    -K NUM       minibatch size for mapping [500M]\n");
//		fprintf(fp_help, "    -v INT       verbose level [%d]\n", mm_verbose);
		fprintf(fp_help, "    --version    show version number\n");
		fprintf(fp_help, "  Preset:\n");
		fprintf(fp_help, "    -x STR       preset (always applied before other options) []\n");
		fprintf(fp_help, "                 map-pb: -Hk19 (PacBio vs reference mapping)\n");
		fprintf(fp_help, "                 map-ont: -k15 (Oxford Nanopore vs reference mapping)\n");
		fprintf(fp_help, "                 asm5: -k19 -w19 -A1 -B19 -O39,81 -E3,1 -s200 -z200 (asm to ref mapping; break at 5%% div.)\n");
		fprintf(fp_help, "                 asm10: -k19 -w19 -A1 -B9 -O16,41 -E2,1 -s200 -z200 (asm to ref mapping; break at 10%% div.)\n");
		fprintf(fp_help, "                 ava-pb: -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 (PacBio read overlap)\n");
		fprintf(fp_help, "                 ava-ont: -k15 -Xw5 -m100 -g10000 -r2000 --max-chain-skip 25 (ONT read overlap)\n");
		fprintf(fp_help, "                 splice: long-read spliced alignment (see minimap2.1 for details)\n");
		fprintf(fp_help, "                 sr: short single-end reads without splicing (see minimap2.1 for details)\n");
		fprintf(fp_help, "\nSee `man ./minimap2.1' for detailed description of command-line options.\n");
		return fp_help == stdout? 0 : 1;
	}

	if ((opt.flag & MM_F_SR) && argc - optind > 3) {
		fprintf(stderr, "[ERROR] incorrect input: in the sr mode, please specify no more than two query files.\n");
		return 1;
	}
	idx_rdr = mm_idx_reader_open(argv[optind], &ipt, fnw);
	if (idx_rdr == 0) {
		fprintf(stderr, "[ERROR] failed to open file '%s'\n", argv[optind]);
		return 1;
	}
	if (!idx_rdr->is_idx && fnw == 0 && argc - optind < 2) {
		fprintf(stderr, "[ERROR] missing input: please specify a query file to map or option -d to keep the index\n");
		mm_idx_reader_close(idx_rdr);
		return 1;
	}
	if (opt.best_n == 0 && (opt.flag&MM_F_CIGAR) && mm_verbose >= 2)
		fprintf(stderr, "[WARNING]\033[1;31m `-N 0' reduces alignment accuracy. Please use --secondary=no to suppress secondary alignments.\033[0m\n");
	
#if FPGA_ON
    fpga_init(BLOCK);
#endif
    
    pthread_t send_tid, recv_tid;
    init_fpga_task_array();
    init_fpga_result_array();
    pthread_create(&send_tid, NULL, send_task_thread, NULL);
    pthread_create(&recv_tid, NULL, recv_task_thread, NULL);
    
    t2 = realtime_msec();
    fprintf(stderr, "init time=%.3f\n", t2 - t1);
    
    pthread_t read_tid, map_tid;
    //pthread_create(&read_tid, NULL, read_task_thread, );
    //pthread_create(&map_tid, NULL, map_task_thread, );
    
    while ((mi = mm_idx_reader_read(idx_rdr, n_threads)) != 0) {
        double t1, t2;
        t1 = realtime_msec();
		if ((opt.flag & MM_F_CIGAR) && (mi->flag & MM_I_NO_SEQ)) {
			fprintf(stderr, "[ERROR] the prebuilt index doesn't contain sequences.\n");
			mm_idx_destroy(mi);
			mm_idx_reader_close(idx_rdr);
			return 1;
		}
        t2 = realtime_msec();
        fprintf(stderr, "handle 1 time:%.3f\n", t2 - t1);
        
        t1 = realtime_msec();
		if ((opt.flag & MM_F_OUT_SAM) && idx_rdr->n_parts == 1) {
			if (mm_idx_reader_eof(idx_rdr)) {
				mm_write_sam_hdr(mi, rg, MM_VERSION, argc, argv);
			} else {
				mm_write_sam_hdr(0, rg, MM_VERSION, argc, argv);
				if (mm_verbose >= 2)
					fprintf(stderr, "[WARNING]\033[1;31m For a multi-part index, no @SQ lines will be outputted.\033[0m\n");
			}
		}
        t2 = realtime_msec();
        fprintf(stderr, "handle 2 time:%.3f\n", t2 - t1);
		//fprintf(stderr,"finished work_post\n");
        t1 = realtime_msec();
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
					__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
		if (argc != optind + 1) mm_mapopt_update(&opt, mi);
        
        t2 = realtime_msec();
        fprintf(stderr, "handle 3 time:%.3f\n", t2 - t1);
        
        t1 = realtime_msec();

		if (mm_verbose >= 3) mm_idx_stat(mi);
#if FPGA_ON
        //bw,  is_cdna,  max_skip,  min_sc,  flag,  max_occ
        fprintf(stderr, "fpga params:bw:0x%x, is_cdna:0x%x, max_skip:0x%x, min_sc:0x%x, flag:0x%x, max_occ:0x%x\n", opt.bw, !!(opt.flag & MM_F_SPLICE), opt.max_chain_skip, opt.min_chain_score, opt.flag, opt.mid_occ);
        fpga_set_params(opt.bw, !!(opt.flag & MM_F_SPLICE), opt.max_chain_skip, opt.min_chain_score, opt.flag, opt.mid_occ);  //data1
#endif
        t2 = realtime_msec();
        fprintf(stderr, "handle 4 time:%.3f\n", t2 - t1);
        
        t1 = realtime_msec();
		if (!(opt.flag & MM_F_FRAG_MODE)) {
			for (i = optind + 1; i < argc; ++i)
				mm_map_file(mi, argv[i], &opt, n_threads);
		} else {
			mm_map_file_frag(mi, argc - (optind + 1), (const char**)&argv[optind + 1], &opt, n_threads);
		}

		mm_idx_destroy(mi);
        t2 = realtime_msec();
        fprintf(stderr, "handle 5 time:%.3f\n", t2 - t1);
	}
    t1 = realtime_msec();
	mm_idx_reader_close(idx_rdr);

    stop_fpga_send_thread();
    stop_fpga_recv_thread();
#if FPGA_ON
    fpga_exit_block();
#endif
    pthread_join(send_tid, NULL);
    pthread_join(recv_tid, NULL);
#if FPGA_ON
    fpga_set_block();
    fpga_finalize();
#endif
	if (fflush(stdout) == EOF) {
		fprintf(stderr, "[ERROR] failed to write the results\n");
		exit(EXIT_FAILURE);
	}

	if (mm_verbose >= 3) {
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, MM_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[M::%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - mm_realtime0, cputime());
	}
    
    double send_task1_total = 0;
    for(i = 0; i < sizeof(send_task1)/sizeof(send_task1[0]); i++) {
        if(send_task1[i] == 0)
            break;
        send_task1_total += send_task1[i];
    }
    fprintf(stderr, "send task 1 time:      %.3f msec i=%d\n", send_task1_total/i, i);
    
    double send_task2_total = 0;
    for(i = 0; i < sizeof(send_task2)/sizeof(send_task2[0]); i++) {
        if(send_task2[i] == 0)
            break;
        send_task2_total += send_task2[i];
    }
    fprintf(stderr, "send task 2 time:      %.3f msec i=%d\n", send_task2_total/i, i);
    
    double process_result_total = 0;
    for(i = 0; i < sizeof(process_result)/sizeof(process_result[0]); i++) {
        if(process_result[i] == 0)
            break;
        process_result_total += process_result[i];
    }
    fprintf(stderr, "process result time:      %.3f msec i=%d\n", process_result_total/i, i);
    
    double soft_chaindp_total = 0;
    int soft_chaindp_tid = 0;
    for(i = 0; i < sizeof(soft_chaindp_time)/sizeof(soft_chaindp_time[0]); i++) {
        if(soft_chaindp_time[i] != 0) {
            soft_chaindp_total += soft_chaindp_time[i];
            soft_chaindp_tid++;
        }
    }
    fprintf(stderr, "soft chaindp total time:      %.3f msec, counter=%d, avg time:%.3f\n", soft_chaindp_total, soft_chaindp_num, soft_chaindp_total/soft_chaindp_num);
    
    fprintf(stderr, "\nmax_task=%d\n", max_task);
    t2 = realtime_msec();
    fprintf(stderr, "finilize time=%.3f\n", t2 - t1);
	return 0;
}
