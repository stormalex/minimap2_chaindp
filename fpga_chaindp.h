#ifndef __FPGA_CHAINDP_H__
#define __FPGA_CHAINDP_H__

#include <stdint.h>
#include "minimap.h"

#define FPGA_ON 1
#define DUMP_FILE   0
#if FPGA_ON
#include "fpga.h"
#endif

#define SEND_ARRAY_MAX  1024

#define INDEX_BUF_SIZE  1024*1024*1024    //256M
#define INDEX_LOAD_SIZE  1024*1024*1024    //4M

#define ADDR_ALIGN(addr, align)   (((addr)+(align)-1)&(~((align)-1)))

struct mm_tbuf_s;
typedef struct _context{
    struct mm_tbuf_s *b;
    const mm_mapopt_t *opt;
    const mm_idx_t *mi;
    unsigned int bid;
    int qlen_sum;
    int max_chain_gap_qry;
    int max_chain_gap_ref;
    int is_splice;
    int n_segs;
    uint32_t hash;
    int *qlens;
    int is_sr;
    char **seqs;
    int *n_regs;
    mm_reg1_t **regs;
    const char *qname;
}context_t;

typedef struct _read_result{
    void *km;
}read_result_t;



typedef struct __attribute__((__packed__)) _collect_task{
    int gap_qry;
    int gap_ref;
    int seednum;
    int qlensum;
    unsigned int read_id;    //文档里的ctxpos，这里填read id
    unsigned int bid;
    int16_t n_segs;
    char b;      //mi->b
    char reserve1[1];
    mm128_t *mv_a;
    char reserve2[28];
}collect_task_t;

typedef struct __attribute__((__packed__)) _collect_result{
    unsigned int err_flag;
    unsigned int read_id;
    unsigned int sub_size;
    unsigned int n_a;
    unsigned int n_minipos;
    unsigned int rep_len;
    char reserve1[40];
}collect_result_t;

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

typedef struct _buf_info{
    void* buf;
    int size;
    int type;
}buf_info_t;

typedef struct _idx_buf idx_buf_t;

int send_fpga_task(buf_info_t task);
int get_fpga_task(buf_info_t* task);

void init_fpga_task_array();
void stop_fpga_send_thread();
void* send_task_thread(void* arg);

void init_fpga_result_array();
void stop_fpga_recv_thread();
void* recv_task_thread(void* arg);

int send_fpga_task(buf_info_t task);
int get_fpga_task(buf_info_t* task);
int send_fpga_result(buf_info_t result);
int get_fpga_result(buf_info_t* result);
#endif //__FPGA_CHAINDP_H__