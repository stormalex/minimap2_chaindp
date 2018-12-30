#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include "kthread.h"
#include "kvec.h"
#include "kalloc.h"
#include "sdust.h"
#include "mmpriv.h"
#include "bseq.h"
#include "khash.h"

#include "fpga_chaindp.h"

void* result_thread(void* args);


struct mm_idx_bucket_s;
struct mm_idx_bucket_s *g_B = NULL;
int32_t g_b = 0;

static int g_parm_flag = 0;
static int g_parm_midocc = 0;
static int g_parm_bw = 0;
static int g_parm_skip = 0;
static int g_parm_issplic = 0;
static int g_parm_score = 0;


struct mm_tbuf_s {
	void *km;
};

mm_tbuf_t *mm_tbuf_init(void)
{
	mm_tbuf_t *b;
	b = (mm_tbuf_t*)calloc(1, sizeof(mm_tbuf_t));
	if (!(mm_dbg_flag & 1)) b->km = km_init();
	return b;
}

void mm_tbuf_destroy(mm_tbuf_t *b)
{
	if (b == 0) return;
	km_destroy(b->km);
	free(b);
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (uint32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && dreg[u]>>32 < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && dreg[v]>>32 < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > dreg[v]>>32? s : dreg[v]>>32;
				int ee = e < (uint32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, j, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

#include "ksort.h"
#define heap_lt(a, b) ((a).x > (b).x)
KSORT_INIT(heap, mm128_t, heap_lt)

typedef struct {
	uint32_t n;
	uint32_t q_pos, q_span;
	uint32_t seg_id:31, is_tandem:1;
	const uint64_t *cr;
} mm_match_t;

static mm_match_t *collect_matches(int *_n_m, int max_occ, mm128_t *mv_a, size_t mv_n, int64_t *n_a, int *rep_len, int *n_mini_pos, uint64_t **mini_pos)
{
	int i, rep_st = 0, rep_en = 0, n_m;
	mm_match_t *m;
	*n_mini_pos = 0;
	*mini_pos = (uint64_t*)malloc(mv_n * sizeof(uint64_t));
	m = (mm_match_t*)malloc(mv_n * sizeof(mm_match_t));
	for (i = n_m = 0, *rep_len = 0, *n_a = 0; i < mv_n; ++i) {
		const uint64_t *cr;
		mm128_t *p = &mv_a[i];
		uint32_t q_pos = (uint32_t)p->y, q_span = p->x & 0xff;
		int t;
		cr = mm_idx_get(p->x>>8, &t);
		if (t >= max_occ) {
			int en = (q_pos >> 1) + 1, st = en - q_span;
			if (st > rep_en) {
				*rep_len += rep_en - rep_st;
				rep_st = st, rep_en = en;
			} else rep_en = en;
		} else {
			mm_match_t *q = &m[n_m++];
			q->q_pos = q_pos, q->q_span = q_span, q->cr = cr, q->n = t, q->seg_id = p->y >> 32;
			q->is_tandem = 0;
			if (i > 0 && p->x>>8 == mv_a[i - 1].x>>8) q->is_tandem = 1;
			if (i < mv_n - 1 && p->x>>8 == mv_a[i + 1].x>>8) q->is_tandem = 1;
			*n_a += q->n;
			(*mini_pos)[(*n_mini_pos)++] = (uint64_t)q_span<<32 | q_pos>>1;
		}
	}
	*rep_len += rep_en - rep_st;
	*_n_m = n_m;
	return m;
}

static inline int skip_seed(int flag, uint64_t r, const mm_match_t *q, unsigned int bid, int qlen, int *is_self)
{
	*is_self = 0;

	if (1 & flag & (MM_F_NO_DIAG|MM_F_NO_DUAL)) {
        int cmp = 0;
#if 1
		uint32_t rankID = (uint32_t)r&0x1FFFFF;
        int flg = (bid &0x80000000)>>31;
        int val = bid  &0x7fffffff;
        if (val > rankID) cmp = 1;
        else if (val < rankID) cmp = -1;
        else if (val == rankID) {
            if (flg) cmp = 0;
            else cmp = -1;
        }
#endif
		//int cmp2 = strcmp(qname, s->name);
        //if (cmp2 > 0) cmp2 = 1;
        //else if(cmp2 < 0) cmp2 = -1;
        //if (cmp != cmp2) fprintf(stderr, "xxxxxxxxxxxxxx,cmp %d, cmp2 %d, qname %s,rname %s\n", qname, s->name);

        //   cmp = cmp2;
		
		if ((flag&MM_F_NO_DIAG) && cmp == 0) {
			if ((((r>>22) & 0x1fffff)) == (q->q_pos>>1)) return 1; // avoid the diagnonal anchors
			if ((r&MM_P_STRAND)>>21 == (q->q_pos&1)) *is_self = 1; // this flag is used to avoid spurious extension on self chain
		}
		if ((flag&MM_F_NO_DUAL) && cmp > 0) // all-vs-all mode: map once
			return 1;
	}
	if (flag & (MM_F_FOR_ONLY|MM_F_REV_ONLY)) {
		if ((r&MM_P_STRAND)>>21 == (q->q_pos&1)) { // forward strand
			if (flag & MM_F_REV_ONLY) return 1;
		} else {
			if (flag & MM_F_FOR_ONLY) return 1;
		}
	}
	return 0;
}

static mm128_t *collect_seed_hits(int flag, int max_occ, mm128_t *mv_a, size_t mv_n, unsigned int bid, int qlen, int64_t *n_a, int *rep_len,
								  int *n_mini_pos, uint64_t **mini_pos)
{
	int i, k, n_m;
	mm_match_t *m;
	mm128_t *a;
	//fprintf(stderr,"collect_matches begin\n");
	m = collect_matches(&n_m, max_occ, mv_a, mv_n, n_a, rep_len, n_mini_pos, mini_pos);
	
	a = (mm128_t*)malloc(*n_a * sizeof(mm128_t));
	for (i = 0, *n_a = 0; i < n_m; ++i) {
		mm_match_t *q = &m[i];
		const uint64_t *r = q->cr;
		for (k = 0; k < q->n; ++k) {
			int32_t is_self;
			//modified by LQX
			int32_t rpos = (r[k]>>22) & 0x1fffff;
			//fprintf(stderr,"get rpos OK\n");
			mm128_t *p;
			if (skip_seed(flag, r[k], q, bid, qlen, &is_self)) continue;
			//fprintf(stderr,"skip_seed OK\n");
			p = &a[(*n_a)++];
			//modified by LQX
			/*
			n=1 val ----63...43|42...22| 21   |20...0
						 refid |refpos |strand|rankid
			n>1 p[k]----63...43|42...22| 21   |20...0
						 refid |refpos |strand|rankid
			*/
			//if ((r[k]&1) == (q->q_pos&1)) { // forward strand
			if((r[k]&MM_P_STRAND)>>21 == (q->q_pos&1)){// forward strand
				//p->x = (r[k]&0xffffffff00000000ULL) | rpos;
				p->x = ((r[k]&0xfffff80000000000ULL)>>11) | rpos;
				p->y = (uint64_t)q->q_span << 32 | q->q_pos >> 1;				
			} else { // reverse strand
				//p->x = 1ULL<<63 | (r[k]&0xffffffff00000000ULL) | rpos;
				p->x = 1ULL<<63 | ((r[k]&0xfffff80000000000ULL)>>11)  | rpos;
				p->y = (uint64_t)q->q_span << 32 | (qlen - ((q->q_pos>>1) + 1 - q->q_span) - 1);
			}
			//fprintf(stderr,"seed assemble OK\n");
			p->y |= (uint64_t)q->seg_id << MM_SEED_SEG_SHIFT;
			if (q->is_tandem) p->y |= MM_SEED_TANDEM;
			if (is_self) p->y |= MM_SEED_SELF;
		}
	}
	free(m);
	radix_sort_128x(a, a + (*n_a));
	//fprintf(stderr,"seed num %d\n",*n_a);
	return a;
}

static void chain_post(const mm_mapopt_t *opt, int max_chain_gap_ref, const mm_idx_t *mi, void *km, int qlen, int n_segs, const int *qlens, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		if (n_segs <= 1) mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		else mm_select_sub_multi(km, opt->pri_ratio, 0.2f, 0.7f, max_chain_gap_ref, mi->k*2, opt->best_n, n_segs, qlens, n_regs, regs);
		if (!(opt->flag & (MM_F_SPLICE|MM_F_SR|MM_F_NO_LJOIN))) // long join not working well without primary chains
			mm_join_long(km, opt, qlen, n_regs, regs, a);
	}
}

static mm_reg1_t *align_regs(const mm_mapopt_t *opt, const mm_idx_t *mi, void *km, int qlen, const char *seq, int *n_regs, mm_reg1_t *regs, mm128_t *a)
{
	if (!(opt->flag & MM_F_CIGAR)) return regs;
	regs = mm_align_skeleton(km, opt, mi, qlen, seq, n_regs, regs, a); // this calls mm_filter_regs()
	if (!(opt->flag & MM_F_ALL_CHAINS)) { // don't choose primary mapping(s)
		mm_set_parent(km, opt->mask_level, *n_regs, regs, opt->a * 2 + opt->b);
		mm_select_sub(km, opt->pri_ratio, mi->k*2, opt->best_n, n_regs, regs);
		mm_set_sam_pri(*n_regs, regs);
	}
	return regs;
}

#include "minimap.h"
unsigned int dichotomy_sort(const char *qname, rname_rid_t* ref_name, int ref_name_size)
{
    unsigned int start=0, end = ref_name_size-1, mid;
    int cmp;
    while(start < end)
        {
            mid = (start + end)>>1;
            cmp = strcmp(qname, ref_name[mid].rname);
            if(cmp == 0)
                return (mid | 1UL<<31);
            else if( cmp < 0 )
                end = mid;
            else
                start = mid+1;
        }
    if (start == end) {
        cmp = strcmp(qname, ref_name[start].rname);
        if(cmp == 0)
            return (start | 1UL<<31);
        else if (cmp > 0) return start+1;
    }
    return start;
}

static void* package_task(collect_task_t** tasks, int num, int size, int* buf_size)
{
    int data_size = size + sizeof(chaindp_sndhdr_t);
    chaindp_sndhdr_t* head = NULL;      //包头指针
    collect_task_t* sub_head = NULL;    //子头指针
    char* p;
    int i = 0;
    if(num == 0)
        return NULL;
    
    void * buf = malloc(data_size);
    p = (char*)buf;
    
    head = (chaindp_sndhdr_t*)buf;
    head->size = data_size;
    head->num = num;
    head->type = 3;
    
    p += sizeof(chaindp_sndhdr_t);      //跳过包头
    sub_head = (collect_task_t*)p;
    
    for(i = 0; i < num; i++) {
        *sub_head = *tasks[i];
        p += sizeof(collect_task_t);    //移动p跳过子包头
        memcpy(p, tasks[i]->mv_a, tasks[i]->seednum * sizeof(mm128_t));   
        
        int tmp_size = (tasks[i]->seednum * sizeof(mm128_t));
        tmp_size = ADDR_ALIGN(tmp_size, 64);
        
        p += tmp_size;                  //移动p跳过数据
        sub_head = (collect_task_t*)p;
    }
    
    assert(p == (buf + data_size));
    
    *buf_size = data_size;
    
    return buf;
}

void mm_map_frag(const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, int *n_regs, mm_reg1_t **regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname, long read_id, user_params_t* params, int tid)
{
	int i, qlen_sum;
	int max_chain_gap_qry, max_chain_gap_ref, is_splice = !!(opt->flag & MM_F_SPLICE), is_sr = !!(opt->flag & MM_F_SR);
	uint32_t hash;
	mm128_v mv = {0,0,0};
    km_stat_t kmst;

    if(read_id < 0 || tid < 0) {
        fprintf(stderr, "invalid read_id(%ld)\n", read_id);
        assert(read_id >= 0);
        assert(tid >= 0);
    }
    
	for (i = 0, qlen_sum = 0; i < n_segs; ++i)
		qlen_sum += qlens[i], n_regs[i] = 0, regs[i] = 0;

	if (qlen_sum == 0 || n_segs <= 0 || n_segs > MM_MAX_SEG) return;

	hash  = qname? __ac_X31_hash_string(qname) : 0;
	hash ^= __ac_Wang_hash(qlen_sum) + __ac_Wang_hash(opt->seed);
	hash  = __ac_Wang_hash(hash);

   
    unsigned int bid = dichotomy_sort(qname, mi->rname_rid, mi->n_seq);
    
	collect_minimizers(b->km, opt, mi, n_segs, qlens, seqs, &mv);
    if(mv.n == 0) {
        __sync_sub_and_fetch(&params->read_num, 1);
        return;
    }
    // set max chaining gap on the query and the reference sequence
	if (is_sr)
		max_chain_gap_qry = qlen_sum > opt->max_gap? qlen_sum : opt->max_gap;
	else max_chain_gap_qry = opt->max_gap;
	if (opt->max_gap_ref > 0) {
		max_chain_gap_ref = opt->max_gap_ref; // always honor mm_mapopt_t::max_gap_ref if set
	} else if (opt->max_frag_len > 0) {
		max_chain_gap_ref = opt->max_frag_len - qlen_sum;
		if (max_chain_gap_ref < opt->max_gap) max_chain_gap_ref = opt->max_gap;
	} else max_chain_gap_ref = opt->max_gap;
    
    //为每一条read创建context_t
    assert(params->read_contexts[read_id] == NULL);
    context_t* context = (context_t*)malloc(sizeof(context_t));
    memset(context, 0, sizeof(context_t));
    params->read_contexts[read_id] = context;
    context->b = b;
    context->opt = opt;
    context->mi = mi;
    context->bid = bid;
    context->qlen_sum = qlen_sum;
    context->max_chain_gap_ref = max_chain_gap_ref;
    context->max_chain_gap_qry = max_chain_gap_qry;
    context->is_splice = is_splice;
    context->n_segs = n_segs;
    context->hash = hash;
    context->qlens = (int*)malloc(MM_MAX_SEG * sizeof(int));
    memcpy(context->qlens, qlens, MM_MAX_SEG * sizeof(int));
    context->is_sr = is_sr;
    context->seqs = (char **)malloc(MM_MAX_SEG * sizeof(char *));
    memcpy(context->seqs, seqs, MM_MAX_SEG * sizeof(char *));
    context->n_regs = n_regs;
    context->regs = regs;
    context->qname = qname;
    
    collect_task_t *task = (collect_task_t*)malloc(sizeof(collect_task_t));
    memset(task, 0, sizeof(collect_task_t));
    task->gap_qry = max_chain_gap_qry;
    task->gap_ref = max_chain_gap_ref;
    task->seednum = mv.n;
    task->qlensum = qlen_sum;
    task->read_id = read_id;
    task->bid = bid;
    task->n_segs = n_segs;
    task->b = mi->b;
    task->mv_a = (mm128_t*)malloc(mv.n * sizeof(mm128_t));
    memcpy(task->mv_a, mv.a, mv.n * sizeof(mm128_t));
    
    if(mv.n == 0) {
        fprintf(stderr, "WARNING:seed num = 0, read_id=%ld\n", read_id);
    }
    
    kfree(b->km, mv.a);
    
    assert(params->tasks[read_id] == NULL);
    params->tasks[read_id] = task;
    
    params->send_task[tid].tasks[params->send_task[tid].num] = task;
    params->send_task[tid].num++;
    params->send_task[tid].data_size += sizeof(collect_task_t);
    int data_size = mv.n * sizeof(mm128_t);
    data_size = ADDR_ALIGN(data_size, 64);
    params->send_task[tid].data_size += data_size;
    
    //打包
    if(params->send_task[tid].num >= 8) {
        int size = 0;
        void* buf = package_task(params->send_task[tid].tasks, params->send_task[tid].num, params->send_task[tid].data_size, &size);
        if(buf) {
            params->send_task[tid].num = 0;
            params->send_task[tid].data_size = 0;
        }
        //发送到fpga发送线程
        buf_info_t buf_info;
        buf_info.buf = buf;
        buf_info.size = size;
        
        while(send_fpga_task(buf_info));
    }
    if (b->km) {
        km_stat(b->km, &kmst);
        if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
            fprintf(stderr, "QM\t%s\t%d\tcap=%ld,nCore=%ld,largest=%ld\n", qname, qlen_sum, kmst.capacity, kmst.n_cores, kmst.largest);
        assert(kmst.n_blocks == kmst.n_cores); // otherwise, there is a memory leak
        if (kmst.largest > 1U<<28) {
            km_destroy(b->km);
            b->km = km_init();
        }
    }
}
void last_send(void *data, int tid)
{
    user_params_t* params = (user_params_t*)data;
    if(params->send_task[tid].num > 0) {
        int size = 0;
        void* buf = package_task(params->send_task[tid].tasks, params->send_task[tid].num, params->send_task[tid].data_size, &size);
        if(buf) {
            params->send_task[tid].num = 0;
            params->send_task[tid].data_size = 0;
        }
        //发送到fpga发送线程
        buf_info_t buf_info;
        buf_info.buf = buf;
        buf_info.size = size;
        
        while(send_fpga_task(buf_info));
    }
}

void* fpga_work(void* buf, int size, int* result_size)
{
    int i = 0;
    
    chaindp_sndhdr_t* head = (chaindp_sndhdr_t*)buf;
    collect_task_t* sub_head = (collect_task_t*)(buf + sizeof(chaindp_sndhdr_t));
    char* p = (char*)sub_head;
    assert(head->type == 3);
    void* out_buf = malloc(16*1024*1024);
    chaindp_sndhdr_t* out_head = (chaindp_sndhdr_t*)out_buf;
    collect_result_t* out_sub_head = (collect_result_t*)(out_buf + sizeof(chaindp_sndhdr_t));
    char* p2 = (char*)out_sub_head;
    
    //输出包的头信息设置
    out_head->magic = head->magic;
    out_head->tid = head->tid;
    out_head->num = head->num;
    out_head->type = head->type;
    out_head->lat = head->lat;
    
    for(i = 0; i < head->num; i++) {
        int tmp_size;
        mm128_t *mv_a = (mm128_t*)(p + sizeof(collect_task_t));
        size_t mv_n = sub_head->seednum;
        unsigned int bid = sub_head->bid;
        int qlen_sum = sub_head->qlensum;
        int64_t n_a = 0;
        int rep_len = 0;
        int n_mini_pos = 0;
        uint64_t* mini_pos = NULL;
        
        int max_chain_gap_ref = sub_head->gap_ref;
        int max_chain_gap_qry = sub_head->gap_qry;
        int n_segs = sub_head->gap_qry;
        uint32_t new_i = 0;
        
        int out_size = 0;
        
        //on fpga
        mm128_t *a = collect_seed_hits(g_parm_flag, g_parm_midocc, mv_a, mv_n, bid, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
        //a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
        struct new_seed* fpga_a = mm_chain_dp_fpga(max_chain_gap_ref, max_chain_gap_qry, g_parm_bw, g_parm_skip, g_parm_score, g_parm_issplic, n_segs, n_a, a, &new_i);
        
        out_sub_head->err_flag = 0;
        out_sub_head->read_id = sub_head->read_id;
        out_sub_head->n_a = new_i;
        out_sub_head->n_minipos = n_mini_pos;
        out_sub_head->rep_len = rep_len;
        
        out_size += sizeof(collect_result_t);
        
        //复制结果的fpga_a数组
        struct new_seed* out_fpga_a = (struct new_seed*)(p2 + sizeof(collect_result_t));
        memcpy(out_fpga_a, fpga_a, new_i * sizeof(struct new_seed));
        free(fpga_a);
        tmp_size = (new_i * sizeof(struct new_seed));
        tmp_size = ADDR_ALIGN(tmp_size, 64);
        p2 += (tmp_size + sizeof(collect_result_t));
        
        out_size += tmp_size;
        
        //复制结果的minipos数组
        uint64_t* out_mini_pos = (uint64_t*)p2;
        memcpy(out_mini_pos, mini_pos, n_mini_pos * sizeof(uint64_t));
        free(mini_pos);
        tmp_size = (n_mini_pos * sizeof(uint64_t));
        tmp_size = ADDR_ALIGN(tmp_size, 64);
        p2 += tmp_size;
        
        out_size += tmp_size;
        out_sub_head->sub_size = out_size;
        
        //移动任务的指针
        tmp_size = (mv_n * sizeof(mm128_t));
        tmp_size = ADDR_ALIGN(tmp_size, 64);
        sub_head = (collect_task_t*)(p + sizeof(collect_task_t) + tmp_size);
        p = (char*)sub_head;
        
        //移动结果的指针
        out_sub_head = (collect_result_t*)p2;
    }
    
    *result_size = ((void*)p2 - out_buf);
    return out_buf;
}

mm_reg1_t *mm_map(const mm_idx_t *mi, int qlen, const char *seq, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt, const char *qname)
{
	mm_reg1_t *regs;
	mm_map_frag(mi, 1, &qlen, &seq, n_regs, &regs, b, opt, qname, -1, NULL, -1);
	return regs;
}

/**************************
 * Multi-threaded mapping *
 **************************/

typedef struct {
	int mini_batch_size, n_processed, n_threads, n_fp;
	const mm_mapopt_t *opt;
	mm_bseq_file_t **fp;
	const mm_idx_t *mi;
	kstring_t str;
} pipeline_t;

typedef struct {
	const pipeline_t *p;
    int n_seq, n_frag;
	mm_bseq1_t *seq;
	int *n_reg, *seg_off, *n_seg;
	mm_reg1_t **reg;
	mm_tbuf_t **buf;
} step_t;

static void worker_for(void *_data, long i, int tid, void* _params) // kt_for() callback
{
    user_params_t* params = (user_params_t*)_params;
    step_t *s = (step_t*)_data;
	int qlens[MM_MAX_SEG], j, off = s->seg_off[i], pe_ori = s->p->opt->pe_ori;
	const char *qseqs[MM_MAX_SEG];
	mm_tbuf_t *b = s->buf[tid];
	assert(s->n_seg[i] <= MM_MAX_SEG);
	if (mm_dbg_flag & MM_DBG_PRINT_QNAME)
		fprintf(stderr, "QR\t%s\t%d\t%d\n", s->seq[off].name, tid, s->seq[off].l_seq);
	for (j = 0; j < s->n_seg[i]; ++j) {
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
			mm_revcomp_bseq(&s->seq[off + j]);
		qlens[j] = s->seq[off + j].l_seq;
		qseqs[j] = s->seq[off + j].seq;
	}
	if (s->p->opt->flag & MM_F_INDEPEND_SEG) {
		for (j = 0; j < s->n_seg[i]; ++j)
			mm_map_frag(s->p->mi, 1, &qlens[j], &qseqs[j], &s->n_reg[off+j], &s->reg[off+j], b, s->p->opt, s->seq[off+j].name, i, params, tid);
	} else {
		mm_map_frag(s->p->mi, s->n_seg[i], qlens, qseqs, &s->n_reg[off], &s->reg[off], b, s->p->opt, s->seq[off].name, i, params, tid);
	}
	for (j = 0; j < s->n_seg[i]; ++j) // flip the query strand and coordinate to the original read strand
		if (s->n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1)))) {
			int k, t;
			mm_revcomp_bseq(&s->seq[off + j]);
			for (k = 0; k < s->n_reg[off + j]; ++k) {
				mm_reg1_t *r = &s->reg[off + j][k];
				t = r->qs;
				r->qs = qlens[j] - r->qe;
				r->qe = qlens[j] - t;
				r->rev = !r->rev;
			}
		}
}

void last_send(void *data, int tid);
void* result_thread(void* args);

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i, j, k;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
		int with_qual = (!!(p->opt->flag & MM_F_OUT_SAM) && !(p->opt->flag & MM_F_NO_QUAL));
		int with_comment = !!(p->opt->flag & MM_F_COPY_COMMENT);
		int frag_mode = (p->n_fp > 1 || !!(p->opt->flag & MM_F_FRAG_MODE));
        step_t *s;
        s = (step_t*)calloc(1, sizeof(step_t));
		if (p->n_fp > 1) s->seq = mm_bseq_read_frag2(p->n_fp, p->fp, p->mini_batch_size, with_qual, with_comment, &s->n_seq);
		else s->seq = mm_bseq_read3(p->fp[0], p->mini_batch_size, with_qual, with_comment, frag_mode, &s->n_seq);
		if (s->seq) {
			s->p = p;
			for (i = 0; i < s->n_seq; ++i)
				s->seq[i].rid = p->n_processed++;
			s->buf = (mm_tbuf_t**)calloc(p->n_threads, sizeof(mm_tbuf_t*));
			for (i = 0; i < p->n_threads; ++i)
				s->buf[i] = mm_tbuf_init();
			s->n_reg = (int*)calloc(3 * s->n_seq, sizeof(int));
			s->seg_off = s->n_reg + s->n_seq; // seg_off and n_seg are allocated together with n_reg
			s->n_seg = s->seg_off + s->n_seq;
			s->reg = (mm_reg1_t**)calloc(s->n_seq, sizeof(mm_reg1_t*));
			for (i = 1, j = 0; i <= s->n_seq; ++i)
				if (i == s->n_seq || !frag_mode || !mm_qname_same(s->seq[i-1].name, s->seq[i].name)) {
					s->n_seg[s->n_frag] = i - j;
					s->seg_off[s->n_frag++] = j;
					j = i;
				}
			return s;
		} else free(s);
    } else if (step == 1) { // step 1: map
        int i = 0;
        user_params_t params;
        memset(&params, 0, sizeof(params));
        params.read_num = ((step_t*)in)->n_frag;
        params.read_contexts = (context_t**)malloc(params.read_num * sizeof(context_t*));
        memset(params.read_contexts, 0, params.read_num * sizeof(context_t*));
        params.read_is_complete = (char*)malloc(params.read_num);
        memset(params.read_is_complete, 0, params.read_num);
        params.read_results = (read_result_t*)malloc(params.read_num * sizeof(read_result_t));
        memset(params.read_results, 0, params.read_num * sizeof(read_result_t));
        params.send_task = (send_task_t*)malloc(p->n_threads * sizeof(send_task_t));
        memset(params.send_task, 0, p->n_threads * sizeof(send_task_t));
        for(i = 0; i < p->n_threads; i++) {
            params.send_task[i].num = 0;
            params.send_task[i].size = SEND_ARRAY_MAX;
            params.send_task[i].data_size = 0;
            params.send_task[i].tasks = (collect_task_t**)malloc(SEND_ARRAY_MAX * sizeof(collect_task_t*));
        }
        params.tasks = (collect_task_t**)malloc(params.read_num * sizeof(collect_task_t*));
        memset(params.tasks, 0, params.read_num * sizeof(collect_task_t*));
        
        params.exit = 1;
        
        pthread_t result_tid[10];
        for(i = 0; i < 10; i++) {
            pthread_create(&result_tid[i], NULL, result_thread, &params);
        }

        g_B = ((step_t*)in)->p->mi->B;      //软件模拟fpga的时候，将B设置到全局变量上
        g_b = ((step_t*)in)->p->mi->b;

        g_parm_flag = p->opt->flag;
        g_parm_midocc = p->opt->mid_occ;
        g_parm_bw = p->opt->bw;
        g_parm_skip = p->opt->max_chain_skip;
        g_parm_issplic = !!(p->opt->flag & MM_F_SPLICE);
        g_parm_score = p->opt->min_chain_score;

		//kt_for_map(p->n_threads, worker_for, in, ((step_t*)in)->n_frag, (void *)&params, last_send, NULL);
        kt_for_map(p->n_threads - 10, worker_for, in, ((step_t*)in)->n_frag, (void *)&params, last_send, result_thread);
        
        for(i = 0; i < 10; i++)
            pthread_join(result_tid[i], NULL);
        
        for(i = 0; i < p->n_threads; i++) {
            free(params.send_task[i].tasks);
        }
        free(params.tasks);
        free(params.send_task);
        free(params.read_results);
        free(params.read_is_complete);
        free(params.read_contexts);
        
		return in;
    } else if (step == 2) { // step 2: output
		void *km = 0;
        step_t *s = (step_t*)in;
		const mm_idx_t *mi = p->mi;
		for (i = 0; i < p->n_threads; ++i) mm_tbuf_destroy(s->buf[i]);
		free(s->buf);
		if ((p->opt->flag & MM_F_OUT_CS) && !(mm_dbg_flag & MM_DBG_NO_KALLOC)) km = km_init();
		for (k = 0; k < s->n_frag; ++k) {
			int seg_st = s->seg_off[k], seg_en = s->seg_off[k] + s->n_seg[k];
			for (i = seg_st; i < seg_en; ++i) {
				mm_bseq1_t *t = &s->seq[i];
				for (j = 0; j < s->n_reg[i]; ++j) {
					mm_reg1_t *r = &s->reg[i][j];
					assert(!r->sam_pri || r->id == r->parent);
					if ((p->opt->flag & MM_F_NO_PRINT_2ND) && r->id != r->parent)
						continue;
					if (p->opt->flag & MM_F_OUT_SAM)
						mm_write_sam2(&p->str, mi, t, i - seg_st, j, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag);
					else
						mm_write_paf(&p->str, mi, t, r, km, p->opt->flag);
					mm_err_puts(p->str.s);
				}
				if (s->n_reg[i] == 0 && (p->opt->flag & MM_F_OUT_SAM)) { // write an unmapped record
					mm_write_sam2(&p->str, mi, t, i - seg_st, -1, s->n_seg[k], &s->n_reg[seg_st], (const mm_reg1_t*const*)&s->reg[seg_st], km, p->opt->flag);
					mm_err_puts(p->str.s);
				}
			}
			for (i = seg_st; i < seg_en; ++i) {
				for (j = 0; j < s->n_reg[i]; ++j) free(s->reg[i][j].p);
				free(s->reg[i]);
				free(s->seq[i].seq); free(s->seq[i].name);
				if (s->seq[i].qual) free(s->seq[i].qual);
			}
		}
		free(s->reg); free(s->n_reg); free(s->seq); // seg_off and n_seg were allocated with reg; no memory leak here
		km_destroy(km);
		if (mm_verbose >= 3)
			fprintf(stderr, "[M::%s::%.3f*%.2f] mapped %d sequences\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), s->n_seq);
		free(s);
	}
    return 0;
}

int mm_map_file_frag(const mm_idx_t *idx, int n_segs, const char **fn, const mm_mapopt_t *opt, int n_threads)
{
	int i, j, pl_threads;
	pipeline_t pl;
	if (n_segs < 1) return -1;
	memset(&pl, 0, sizeof(pipeline_t));
	pl.n_fp = n_segs;
	pl.fp = (mm_bseq_file_t**)calloc(n_segs, sizeof(mm_bseq_file_t*));
	for (i = 0; i < n_segs; ++i) {
		pl.fp[i] = mm_bseq_open(fn[i]);
		if (pl.fp[i] == 0) {
			if (mm_verbose >= 1)
				fprintf(stderr, "ERROR: failed to open file '%s'\n", fn[i]);
			for (j = 0; j < i; ++j)
				mm_bseq_close(pl.fp[j]);
			free(pl.fp);
			return -1;
		}
	}
	pl.opt = opt, pl.mi = idx;
	pl.n_threads = n_threads > 1? n_threads : 1;
	pl.mini_batch_size = opt->mini_batch_size;
	pl_threads = n_threads == 1? 1 : (opt->flag&MM_F_2_IO_THREADS)? 3 : 2;
	kt_pipeline(pl_threads, worker_pipeline, &pl, 3);
	free(pl.str.s);
	for (i = 0; i < n_segs; ++i)
		mm_bseq_close(pl.fp[i]);
	free(pl.fp);
	return 0;
}

int mm_map_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads)
{
	return mm_map_file_frag(idx, 1, &fn, opt, n_threads);
}

int read_result_handle(context_t* context, int n_minipos, uint64_t* mini_pos, struct new_seed* fpga_a, uint32_t new_i, int rep_len, unsigned int read_id)
{
    int i, j;
    mm128_t *a = NULL;
    int n_regs0;
    uint64_t *u;
    mm_reg1_t *regs0;
    const mm_mapopt_t *opt = context->opt;
    int n_segs = context->n_segs;
    //struct mm_tbuf_s *b = context->b;
    uint32_t hash = context->hash;
    int qlen_sum = context->qlen_sum;
    const mm_idx_t *mi = context->mi;
    int max_chain_gap_ref = context->max_chain_gap_ref;
    const int *qlens = context->qlens;
    int is_sr = context->is_sr;
    char **seqs = context->seqs;
    int *n_regs = context->n_regs;
    mm_reg1_t **regs = context->regs;
    //const char *qname = context->qname;
    
    a = mm_chain_dp_bottom(opt->min_cnt, opt->min_chain_score, n_segs, &n_regs0, &u, NULL, fpga_a, new_i);

    if (opt->max_occ > opt->mid_occ && rep_len > 0) {
        assert(0);
        /*int rechain = 0;
        if (n_regs0 > 0) { // test if the best chain has all the segments
            int n_chained_segs = 1, max = 0, max_i = -1, max_off = -1, off = 0;
            for (i = 0; i < n_regs0; ++i) { // find the best chain
                if (max < u[i]>>32) max = u[i]>>32, max_i = i, max_off = off;
                off += (uint32_t)u[i];
            }
            for (i = 1; i < (uint32_t)u[max_i]; ++i) // count the number of segments in the best chain
                if ((a[max_off+i].y&MM_SEED_SEG_MASK) != (a[max_off+i-1].y&MM_SEED_SEG_MASK))
                    ++n_chained_segs;
            if (n_chained_segs < n_segs)
                rechain = 1;
        } else rechain = 1;
        if (rechain) { // redo chaining with a higher max_occ threshold
            //if (opt->flag & MM_F_HEAP_SORT) a = collect_seed_hits_heap(b->km, opt, opt->max_occ, mi, qname, &mv, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
            //else 
            a = collect_seed_hits(opt->flag, opt->max_occ, &mv, bid, qlen_sum, &n_a, &rep_len, &n_mini_pos, &mini_pos);
            a = mm_chain_dp(max_chain_gap_ref, max_chain_gap_qry, opt->bw, opt->max_chain_skip, opt->min_cnt, opt->min_chain_score, is_splice, n_segs, n_a, a, &n_regs0, &u, b->km);
        }*/
    }

    regs0 = mm_gen_regs(NULL, hash, qlen_sum, n_regs0, u, a);
    //fprintf(stderr,"finish mm_gen_regs one read!\n");
    if (mm_dbg_flag & MM_DBG_PRINT_SEED)
        for (j = 0; j < n_regs0; ++j)
            for (i = regs0[j].as; i < regs0[j].as + regs0[j].cnt; ++i)
                fprintf(stderr, "CN\t%d\t%s\t%d\t%c\t%d\t%d\t%d\n", j, mi->seq[a[i].x<<1>>33].name, (int32_t)a[i].x, "+-"[a[i].x>>63], (int32_t)a[i].y, (int32_t)(a[i].y>>32&0xff),
                        i == regs0[j].as? 0 : ((int32_t)a[i].y - (int32_t)a[i-1].y) - ((int32_t)a[i].x - (int32_t)a[i-1].x));

    chain_post(opt, max_chain_gap_ref, mi, NULL, qlen_sum, n_segs, qlens, &n_regs0, regs0, a);
    
    if (!is_sr) mm_est_err(mi, qlen_sum, n_regs0, regs0, a, n_minipos, mini_pos);

    if (n_segs == 1) { // uni-segment
        regs0 = align_regs(opt, mi, NULL, qlens[0], seqs[0], &n_regs0, regs0, a);
        mm_set_mapq(NULL, n_regs0, regs0, opt->min_chain_score, opt->a, rep_len, is_sr);
        n_regs[0] = n_regs0, regs[0] = regs0;
    } else { // multi-segment
        mm_seg_t *seg;
        seg = mm_seg_gen(NULL, hash, n_segs, qlens, n_regs0, regs0, n_regs, regs, a); // split fragment chain to separate segment chains
        free(regs0);
        for (i = 0; i < n_segs; ++i) {
            mm_set_parent(NULL, opt->mask_level, n_regs[i], regs[i], opt->a * 2 + opt->b); // update mm_reg1_t::parent
            regs[i] = align_regs(opt, mi, NULL, qlens[i], seqs[i], &n_regs[i], regs[i], seg[i].a);
            mm_set_mapq(NULL, n_regs[i], regs[i], opt->min_chain_score, opt->a, rep_len, is_sr);
        }
        mm_seg_free(NULL, n_segs, seg);
        if (n_segs == 2 && opt->pe_ori >= 0 && (opt->flag&MM_F_CIGAR))
            mm_pair(NULL, max_chain_gap_ref, opt->pe_bonus, opt->a * 2 + opt->b, opt->a, qlens, n_regs, regs); // pairing
    }
    free(context->qlens);
    free(context->seqs);
    if(a != NULL)
        free(a);
    if(u != NULL)
        free(u);
    return 0;
}

void* result_thread(void* args)
{
    int ret;
    int i = 0;
    user_params_t *params = (user_params_t *)args;
    
    buf_info_t result;
    
    while(params->exit) {
        ret = get_fpga_result(&result);
        if(ret == 1) {
            continue;
        }
        
        chaindp_sndhdr_t* out_head = (chaindp_sndhdr_t*)result.buf;
        collect_result_t* out_sub_head = (collect_result_t*)(result.buf + sizeof(chaindp_sndhdr_t));
    
    
        for(i = 0; i < out_head->num; i++) {
            char* p = (char*)out_sub_head;
            int n_minipos = out_sub_head->n_minipos;
            uint32_t new_i = out_sub_head->n_a;
            struct new_seed* fpga_a = (struct new_seed*)(p + sizeof(collect_result_t));
            int tmp_size = new_i * sizeof(struct new_seed);
            tmp_size = ADDR_ALIGN(tmp_size, 64);
            p = p + sizeof(collect_result_t) + tmp_size;
            uint64_t* mini_pos = (uint64_t*)p;
            tmp_size = n_minipos * sizeof(uint64_t);
            tmp_size = ADDR_ALIGN(tmp_size, 64);
            p += tmp_size;
            unsigned int read_id = out_sub_head->read_id;
            int rep_len = out_sub_head->rep_len;
            
            if(out_sub_head->err_flag == 1) {   //硬件处理不了的数据，软件处理
                fprintf(stderr, "soft process, read_id=%u\n", read_id);
                collect_task_t* task = params->tasks[read_id];
                int64_t n_a = 0;
                //on fpga
                mm128_t *a = collect_seed_hits(g_parm_flag, g_parm_midocc, task->mv_a, task->seednum, task->bid, task->qlensum, &n_a, &rep_len, &n_minipos, &mini_pos);
                fpga_a = mm_chain_dp_fpga(task->gap_ref, task->gap_qry, g_parm_bw, g_parm_skip, g_parm_score, g_parm_issplic, task->n_segs, n_a, a, &new_i);
            }
            
            context_t* context = params->read_contexts[read_id];
            
            ret = read_result_handle(context, n_minipos, mini_pos, fpga_a, new_i, rep_len, read_id);
            if(ret == 0) {  //正确处理，判断所有read是否处理完成
                if(out_sub_head->err_flag == 1) {   //如果是软件处理，则将生成的结果释放
                    free(fpga_a);
                    free(mini_pos);
                }
                //销毁上下文
                free(context);
                params->read_contexts[read_id] = NULL;
                
                //销毁任务
                free(params->tasks[read_id]->mv_a);
                free(params->tasks[read_id]);
                params->tasks[read_id] = NULL;
                
                long read_num_tmp = __sync_sub_and_fetch(&params->read_num, 1);
                if(read_num_tmp == 0) { //所有read处理完成
                    free(result.buf);   //释放结果buf
                    params->exit = 0;
                    fprintf(stderr, "1.exit result_thread\n");
                    return NULL;
                }
            }
            if(out_sub_head->err_flag == 1)
                out_sub_head++;
            else
                out_sub_head = (collect_result_t*)p;
        }
        free(result.buf);   //释放结果buf
    }
    fprintf(stderr, "2.exit result_thread\n");
    return NULL;
}
