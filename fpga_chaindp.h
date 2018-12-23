#ifndef __FPGA_CHAINDP_H__
#define __FPGA_CHAINDP_H__

#define SEND_ARRAY_MAX  1024

typedef struct _context{
    void *km;
    const mm_mapopt_t *opt;
    const mm_idx_t *mi;
    unsigned int bid;
    int qlen;
    int max_chain_gap_qry;
    int max_chain_gap_ref;
    int is_splice;
    int n_segs;
}context_t;

typedef struct _read_result{
    void *km;
}read_result_t;

typedef struct _send_task{
    int num;
    int size;
}send_task_t;

#endif //__FPGA_CHAINDP_H__