#ifndef __FPGA_CHAINDP_H__
#define __FPGA_CHAINDP_H__

#define SEND_ARRAY_MAX  1024

#define ADDR_ALIGN(addr, align)   (((addr)+(align)-1)&(~((align)-1)))

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



typedef struct __attribute__((__packed__)) _collect_task{
    int gap_qry;
    int gap_ref;
    int seednum;
    int qlensum;
    int read_id;    //文档里的ctxpos，这里填read id
    unsigned int bid;
    int16_t n_segs;
    char b;      //mi->b
    char reserve1[1];
    mm128_t *mv_a;
    char reserve2[28];
}collect_task_t;

typedef struct _send_task{
    int num;
    int size;
    int data_size;
    collect_task_t** tasks;
}send_task_t;


//发给fpga的包格式
typedef struct __attribute__((__packed__)) _chaindp_sndhdr{
    unsigned int magic;
    unsigned int size;
    uint16_t tid;
    uint16_t num;       //包里read的数量
    uint8_t  type;      //3
    uint8_t lat;
    uint8_t reserve1[50];
}chaindp_sndhdr_t;

#endif //__FPGA_CHAINDP_H__