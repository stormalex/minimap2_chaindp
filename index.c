#include <stdlib.h>
#include <assert.h>
#if defined(WIN32) || defined(_WIN32)
#include <io.h> // for open(2)
#else
#include <unistd.h>
#endif
#include <fcntl.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#define __STDC_LIMIT_MACROS
#include "kthread.h"
#include "bseq.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kvec.h"
#include "khash.h"

#include "fpga_chaindp.h"

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
typedef khash_t(idx) idxhash_t;

KHASH_MAP_INIT_STR(str, uint32_t)

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

extern struct mm_idx_bucket_s *g_B;
extern int32_t g_b;

#if FPGA_ON
#include <sys/time.h>
#include<time.h>
static double realtime_msec(void)
{
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return tp.tv_sec*1000 + tp.tv_nsec*1e-6;
}
#endif

idx_buf_t* idx_buf_create()
{
    idx_buf_t* idx_buf = (idx_buf_t*)malloc(sizeof(idx_buf_t));
    if(idx_buf) {
        idx_buf->buf = (char*)malloc(INDEX_BUF_SIZE);
        idx_buf->size = INDEX_BUF_SIZE;
        idx_buf->pos = 0;
    }
    return idx_buf;
}

int idx_buf_write64(idx_buf_t* idx_buf, uint8_t* data, int n)
{
    int i = 0;
    if(idx_buf->pos + n >= idx_buf->size) {
        idx_buf->size += INDEX_BUF_SIZE;
        idx_buf->buf = (char*)realloc(idx_buf->buf, idx_buf->size);
        if(idx_buf->buf == NULL) {
            fprintf(stderr, "realloc failed\n");
            exit(1);
        }
    }
    for(i = 0; i < n; i++) {
        idx_buf->buf[idx_buf->pos++] = data[i];
    }
    return n;
}

void idx_buf_write_to_file(idx_buf_t* idx_buf, char* file_name)
{
    FILE* fp = fopen(file_name, "a+");
    if(fp) {
        int ret = fwrite(idx_buf->buf, idx_buf->pos, 1, fp);
        if(ret != 1) {
            fprintf(stderr, "write %s failed:%s\n", file_name, strerror(errno));
            exit(1);
        }
        fclose(fp);
        idx_buf->pos = 0;
    }
    return;
}

void idx_buf_destroy(idx_buf_t* idx_buf)
{
    free(idx_buf->buf);
    free(idx_buf);
}

void load_index_buf(idx_buf_t* idx_buf, int type)
{
    unsigned long long size = idx_buf->pos;
    unsigned long long i = 0;
    
    while(size > INDEX_LOAD_SIZE) {
        fpga_load_index(&idx_buf->buf[i], INDEX_LOAD_SIZE, type);
        size -= INDEX_LOAD_SIZE;
        i += INDEX_LOAD_SIZE;
    }
    if(size > 0) {
        unsigned long long tmp_size = ADDR_ALIGN(size, 64);
        fpga_load_index(&idx_buf->buf[i], tmp_size, type);
        i += size;
    }
    fprintf(stderr, "load size:%lld buf type:%d\n", i, type);
    return;
}

void load_index(char* file_name, int type)
{
    int ret = 0;
    struct stat stat;
    char* buf = NULL;
    long size = 0;
    
    fprintf(stderr, "load index:%s\n", file_name);
    FILE* fp = fopen(file_name, "r");
    if(fp == NULL) {
        perror("fopen failed\n");
        exit(1);
    }
    
    ret = lstat(file_name, &stat);
    if(ret) {
        perror("ERROR:fstat failed!");
        exit(1);
    }
    //size = ((stat.st_size + 64 - 1) & (~(64 - 1)));
    size = stat.st_size;
    while(size > INDEX_LOAD_SIZE) {
        buf = (char*)malloc(INDEX_LOAD_SIZE);
        memset(buf, 0, INDEX_LOAD_SIZE);
        ret = fread(buf, 1, INDEX_LOAD_SIZE, fp);
        if(ret != INDEX_LOAD_SIZE) {
            perror("1.fread failed");
            exit(1);
        }
        fpga_load_index(buf, INDEX_LOAD_SIZE, type);
        free(buf);
        size -= INDEX_LOAD_SIZE;
    }
    if(size > 0) {
        long tmp_size = ADDR_ALIGN(size, 64);
        buf = (char*)malloc(tmp_size);
        memset(buf, 0, tmp_size);
        ret = fread(buf, 1, size, fp);
        if(ret != size) {
            perror("2.fread failed");
            exit(1);
        }
        fpga_load_index(buf, tmp_size, type);
        free(buf);
    }
    
    fclose(fp);
    return;
}

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
    
    mi->b_idx = idx_buf_create();
    mi->h_idx = idx_buf_create();
    mi->v_idx = idx_buf_create();
    mi->p_idx = idx_buf_create();
    
	return mi;
}

void mm_idx_destroy(mm_idx_t *mi)
{
	int i;
	if (mi == 0) return;
    
    if(mi->b_idx)
        idx_buf_destroy(mi->b_idx);
    if(mi->h_idx)
        idx_buf_destroy(mi->h_idx);
    if(mi->v_idx)
        idx_buf_destroy(mi->v_idx);
    if(mi->p_idx)
        idx_buf_destroy(mi->p_idx);
    
    if(mi->rever_rid)
        free(mi->rever_rid);
    if(mi->rname_rid)
        free(mi->rname_rid);
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	for (i = 0; i < 1<<mi->b; ++i) {
		free(mi->B[i].p);
		free(mi->B[i].a.a);
		kh_destroy(idx, (idxhash_t*)mi->B[i].h);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(uint64_t minier, int *n)
{
	int mask = (1<<g_b) - 1;
	khint_t k;
	mm_idx_bucket_t *b = &g_B[minier&mask];
	idxhash_t *h = (idxhash_t*)b->h;
	*n = 0;
	if (h == 0) return 0;
	k = kh_get(idx, h, minier>>g_b<<1);
	if (k == kh_end(h)) return 0;
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int i, n = 0, n1 = 0;
	uint64_t sum = 0, len = 0;
    double t1, t2;
    t1 = realtime_msec();
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum);
    t2 = realtime_msec();
    fprintf(stderr, "mm_idx_stat time:%.3f\n", t2 - t1);
}

int mm_idx_index_name(mm_idx_t *mi)
{
	khash_t(str) *h;
	uint32_t i;
	int has_dup = 0, absent;
	if (mi->h) return 0;
	h = kh_init(str);
	for (i = 0; i < mi->n_seq; ++i) {
		khint_t k;
		k = kh_put(str, h, mi->seq[i].name, &absent);
		if (absent) kh_val(h, k) = i;
		else has_dup = 1;
	}
	mi->h = h;
	if (has_dup && mm_verbose >= 2)
		fprintf(stderr, "[WARNING] some database sequences have identical sequence names\n");
	return has_dup;
}

int mm_idx_name2id(const mm_idx_t *mi, const char *name)
{
	khash_t(str) *h = (khash_t(str)*)mi->h;
	khint_t k;
	if (h == 0) return -2;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

int mm_idx_getseq(const mm_idx_t *mi, uint32_t rid, uint32_t st, uint32_t en, uint8_t *seq)
{
	uint64_t i, st1, en1;
	if (rid >= mi->n_seq || st >= mi->seq[rid].len) return -1;
	if (en > mi->seq[rid].len) en = mi->seq[rid].len;
	st1 = mi->seq[rid].offset + st;
	en1 = mi->seq[rid].offset + en;
	for (i = st1; i < en1; ++i)
		seq[i - st1] = mm_seq4_get(mi->S, i);
	return en - st;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

int cmp(const void *r0, const void *r1)
{
    rname_rid_t *a = (rname_rid_t *)r0;
    rname_rid_t *b = (rname_rid_t *)r1;
    return strcmp(a->rname, b->rname);
}
/*********************************
 * Sort and generate hash tables *
 *********************************/

static void worker_post(void *g, long i, int tid)
{
	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];
	if (b->a.n == 0) return;

	// sort by minimizer
	radix_sort_128x(b->a.a, b->a.a + b->a.n);

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;
			if (n > 1) b->n += n;
			n = 1;
		} else ++n;
	}
	h = kh_init(idx);
	kh_resize(idx, h, n_keys);
	b->p = (uint64_t*)calloc(b->n, 8);

	//modified by LQX

	// create the hash table
	/*
	n=1 val ----63...43|42...22| 21   |20...0
				 refid |refpos |strand|rankid
	n>1 p[k]----63...43|42...22| 21   |20...0
				 refid |refpos |strand|rankid
	*/
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
			assert(absent && j - start_a == n);
			if (n == 1) {
				kh_key(h, itr) |= 1;
				//uint64_t val = p->y;
				uint64_t refid = (p->y)>>32;
				uint64_t refpos_strand = (p->y)&0xFFFFFFFF;
				uint32_t rankid = mi->rever_rid[refid];
				uint64_t val = ((refid&0x1FFFFF)<<43)|((refpos_strand&0x3FFFFF)<<21)|((rankid&0x1FFFFF));
				kh_val(h, itr) = val;			
				//fprintf(stderr,"finished work_post_n_1 \n");
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				for (k = 0; k < n; ++k){
					uint64_t refid = (b->p[start_p+k])>>32;
					uint64_t refpos_strand = (b->p[start_p+k])&0xFFFFFFFF;
					uint32_t rankid = mi->rever_rid[refid];
					uint64_t pk = ((refid&0x1FFFFF)<<43)|((refpos_strand&0x3FFFFF)<<21)|((rankid&0x1FFFFF));
					b->p[start_p+k] = pk;
				}
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;
				start_p += n;
				//fprintf(stderr,"finished work_post\n");
				//modified by LQX
				
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
	//fprintf(stderr,"finished work_post\n");
}
 
static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

/******************
 * Generate index *
 ******************/

#include <string.h>
#include <zlib.h>
#include "bseq.h"

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	int i, mask = (1<<mi->b) - 1;
	for (i = 0; i < n; ++i) {
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);
	}
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;
        s = (step_t*)calloc(1, sizeof(step_t));
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						mm_seq4_set(p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: compute sketch
        step_t *s = (step_t*)in;
		for (i = 0; i < s->n_seq; ++i) {
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a);
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // dispatch sketch to buckets
        step_t *s = (step_t*)in;
		mm_idx_add(p->mi, s->a.n, s->a.a);
		kfree(0, s->a.a); free(s);
	}
    return 0;
}


mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	pipeline_t pl;
    double t1, t2;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
    t1 = realtime_msec();
	memset(&pl, 0, sizeof(pipeline_t));
	pl.mini_batch_size = mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
    t2 = realtime_msec();
    fprintf(stderr, "mm_idx_gen 1:%.3f msec\n", t2 - t1);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));
	#if 1
    t1 = realtime_msec();
	{
		mm_idx_t *mi = pl.mi;
		mi->rname_rid = (rname_rid_t *)calloc(mi->n_seq, sizeof(rname_rid_t));
		mi->rever_rid = (int *)calloc(mi->n_seq, sizeof(int));

		//fprintf(stderr, "---------------------\n");
		for (int k=0; k<mi->n_seq; k++)
			{
				const mm_idx_seq_t *s = &mi->seq[k];
				strncpy(mi->rname_rid[k].rname, s->name, 255);
				mi->rname_rid[k].rname[strlen(s->name)] = 0;
				mi->rname_rid[k].rid = k;
			}
		
		qsort(mi->rname_rid, mi->n_seq, sizeof(rname_rid_t), cmp);
		/*
		fprintf(stderr, "---------------------\n");
		for (int k=0; k<mi->n_seq; k++)
			{
				fprintf(stderr, "%s\n", mi->rname_rid[k].rname);
			}
		*/
		for (int k=0; k<mi->n_seq; k++)
			{
				mi->rever_rid[mi->rname_rid[k].rid] = k;
			}
		/*
		for (int k=0; k<mi->n_seq; k++)
			{
				fprintf(stderr, "table,%d:%d\n", k, mi->rever_rid[k]);
			}
		*/
	}
    t2 = realtime_msec();
    fprintf(stderr, "mm_idx_gen 2:%.3f msec\n", t2 - t1);
	#endif
    
    t1 = realtime_msec();
    
	mm_idx_post(pl.mi, n_threads);
    
    t2 = realtime_msec();
    fprintf(stderr, "mm_idx_gen 3:%.3f msec\n", t2 - t1);
#if FPGA_ON

    //system("rm *.dat -f");

    double start, end;
    mm_idx_t *mi = pl.mi;
    
    start = realtime_msec();
    
    uint64_t allh = 0;
    uint64_t allp = 0;
    int i;
    for(i = 0;i < (1<<mi->b);++i){
        kh_idx_t* h = (idxhash_t*)mi->B[i].h;
        uint64_t B[2] = {0};
        if (mi->B[i].h) {
            uint64_t values;
            //Get B Array
            uint64_t n_buckets;
            n_buckets = kh_end((idxhash_t*)mi->B[i].h);
            B[1] = allh << 28 | allp >> 8;
            B[0] =(((allp&0xFF)<<56) | (n_buckets<<24));
            //if(idx_buf_write64(b_buf, (uint8_t *)B, 16) == 0) {
                //idx_buf_write_to_file(b_buf, "idxb.dat");
                idx_buf_write64(mi->b_idx, (uint8_t *)B, 16);
            //}

            allh += n_buckets;
            allp += mi->B[i].n;
            //将 n_buckets h p 合并到128bit中，实际是拆分两个64bit 顺序分别为h-36bit|p-28bit|p-8bit|n_buckets-32bit|padding
            for(int j = 0;j < n_buckets;++j){
                uint64_t flags;
                uint64_t keys;
                flags = h->flags[j>>4];
                keys = h->keys[j];
                uint64_t a[2]={0};

                //Get H Array
                a[1] |= flags;
                a[1] = (a[1] << 32) |((keys&0xffffffffffffUL)>>16);
                a[0] |= keys&0xffff;
                a[0] = a[0]<<48;
                //if(idx_buf_write64(h_buf, (uint8_t *)a, 16) == 0) {
                    //idx_buf_write_to_file(h_buf, "idxh.dat");
                    idx_buf_write64(mi->h_idx, (uint8_t *)a, 16);
                //}

                //Get Values Array
                values = h->vals[j];//value是后期可能要优化为48bit    21-bit|22-bit|padding
                //if(idx_buf_write64(v_buf, (uint8_t *)&values, 8) == 0) {
                    //idx_buf_write_to_file(v_buf, "idxv.dat");
                    idx_buf_write64(mi->v_idx, (uint8_t *)&values, 8);
                //}
            }

            //Get p Array
            int n = mi->B[i].n;
            for(int j = 0;j < n;++j){
                values = mi->B[i].p[j];
                //if(idx_buf_write64(p_buf, (uint8_t *)&values, 8) == 0) {
                    //idx_buf_write_to_file(p_buf, "idxp.dat");
                    idx_buf_write64(mi->p_idx, (uint8_t *)&values, 8);
                //}
            }
        } else {
            //if(idx_buf_write64(b_buf, (uint8_t *)B, 16) == 0) {
                //idx_buf_write_to_file(b_buf, "idxb.dat");
                idx_buf_write64(mi->b_idx, (uint8_t *)B, 16);
            //}
        }
    }
    //idx_buf_write_to_file(b_buf, "idxb.dat");
    //idx_buf_write_to_file(h_buf, "idxh.dat");
    //idx_buf_write_to_file(p_buf, "idxp.dat");
    //idx_buf_write_to_file(v_buf, "idxv.dat");
    
    
    //load_index("idxb.dat", TYPE_INDEX_B);
    //load_index("idxh.dat", TYPE_INDEX_H);
    //load_index("idxv.dat", TYPE_INDEX_V);
    //load_index("idxp.dat", TYPE_INDEX_P);
    
    //load_index_buf(mi->b_idx, TYPE_INDEX_B);
    //load_index_buf(mi->h_idx, TYPE_INDEX_H);
    //load_index_buf(mi->v_idx, TYPE_INDEX_V);
    //load_index_buf(mi->p_idx, TYPE_INDEX_P);
    
    end = realtime_msec();
    fprintf(stderr, "create idx time:%.3f\n", end - start);
#endif
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_build(const char *fn, int w, int k, int flag, int n_threads) // a simpler interface; deprecated
{
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
	fp = mm_bseq_open(fn);
	if (fp == 0) return 0;
	mi = mm_idx_gen(fp, w, k, 14, flag, 1<<18, n_threads, UINT64_MAX);
	mm_bseq_close(fp);
	return mi;
}

mm_idx_t *mm_idx_str(int w, int k, int is_hpc, int bucket_bits, int n, const char **seq, const char **name)
{
	uint64_t sum_len = 0;
	mm128_v a = {0,0,0};
	mm_idx_t *mi;
	int i, flag = 0;
	if (n <= 0) return 0;
	for (i = 0; i < n; ++i) // get the total length
		sum_len += strlen(seq[i]);
	if (is_hpc) flag |= MM_I_HPC;
	if (name == 0) flag |= MM_I_NO_NAME;
	if (bucket_bits < 0) bucket_bits = 14;
	mi = mm_idx_init(w, k, bucket_bits, flag);
	mi->n_seq = n;
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, n, sizeof(mm_idx_seq_t)); // ->seq is allocated from km
	mi->S = (uint32_t*)calloc((sum_len + 7) / 8, 4);
	for (i = 0, sum_len = 0; i < n; ++i) {
		const char *s = seq[i];
		mm_idx_seq_t *p = &mi->seq[i];
		uint32_t j;
		if (name && name[i]) {
			p->name = (char*)kmalloc(mi->km, strlen(name[i]) + 1);
			strcpy(p->name, name[i]);
		}
		p->offset = sum_len;
		p->len = strlen(s);
		for (j = 0; j < p->len; ++j) {
			int c = seq_nt4_table[(uint8_t)s[j]];
			uint64_t o = sum_len + j;
			mm_seq4_set(mi->S, o, c);
		}
		sum_len += p->len;
		if (p->len > 0) {
			a.n = 0;
			mm_sketch(0, s, p->len, w, k, i, is_hpc, &a);
			mm_idx_add(mi, a.n, a.a);
		}
	}
	free(a.a);
	mm_idx_post(mi, 1);
	return mi;
}

/*************
 * index I/O *
 *************/

void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5];
	int i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		l = strlen(mi->seq[i].name);
		fwrite(&l, 1, 1, fp);
		fwrite(mi->seq[i].name, 1, l, fp);
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	int i;
	char magic[4];
	uint32_t x[5];
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 5, fp) != 5) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		s->name = (char*)kmalloc(mi->km, l + 1);
		fread(s->name, 1, l, fp);
		s->name[l] = 0;
		fread(&s->len, 4, 1, fp);
		s->offset = sum_len;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
	return mi;
}

int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	off_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	is_idx = mm_idx_is_idx(fn);
	if (is_idx < 0) return 0; // failed to open the index
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
    double t1, t2;
	mm_idx_t *mi;
    t1 = realtime_msec();
	if (r->is_idx) {
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	if (mi) {
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		++r->n_parts;
	}
    t2 = realtime_msec();
    fprintf(stderr, "mm_idx_reader_read time:%.3f\n", t2 - t1);
	return mi;
}

int mm_idx_reader_eof(const mm_idx_reader_t *r) // TODO: in extremely rare cases, mm_bseq_eof() might not work
{
	return r->is_idx? (feof(r->fp.idx) || ftell(r->fp.idx) == r->idx_size) : mm_bseq_eof(r->fp.seq);
}
