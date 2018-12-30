#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#include "fpga_chaindp.h"

//保存发送任务的队列
#define QUEUE_MAX    4096
static buf_info_t tasks[QUEUE_MAX];
static pthread_mutex_t tasks_mutex = PTHREAD_MUTEX_INITIALIZER;
static int tasks_head = 0;
static int tasks_tail = 0;

static int fpga_send_task_stop = 1;

void* fpga_work(void* buf, int size, int* result_size);

void init_fpga_task_array()
{
    memset(tasks, 0, sizeof(tasks));
    tasks_head = 0;
    tasks_tail = 0;
    fpga_send_task_stop = 1;
    return;
}

void stop_fpga_send_thread()
{
    fpga_send_task_stop = 0;
}

static int fpga_task_is_full()
{
    return ((tasks_tail + 1) % QUEUE_MAX == tasks_head);
}
static int fpga_task_is_empty()
{
    return (tasks_head == tasks_tail);
}

int send_fpga_task(buf_info_t task)
{
    pthread_mutex_lock(&tasks_mutex);
    if (fpga_task_is_full()) {
        pthread_mutex_unlock(&tasks_mutex);
        return -1;
    }
    tasks[tasks_tail] = task;
    tasks_tail = (tasks_tail + 1) % QUEUE_MAX;
        
    pthread_mutex_unlock(&tasks_mutex);
    return 0;
}

int get_fpga_task(buf_info_t* task)
{
    pthread_mutex_lock(&tasks_mutex);
    if (fpga_task_is_empty()) {
        pthread_mutex_unlock(&tasks_mutex);
        return 1;
    }
    *task = tasks[tasks_head];
    tasks_head = (tasks_head + 1) % QUEUE_MAX;
    pthread_mutex_unlock(&tasks_mutex);
    return 0;
}

void* send_task_thread(void* arg)
{
    int magic_idx = 0;
    buf_info_t task;
    
    while(fpga_send_task_stop) {
        if(get_fpga_task(&task)) {
            continue;
        }
        chaindp_sndhdr_t* head = (chaindp_sndhdr_t*)task.buf;      //包头指针
        head->magic = magic_idx;
        magic_idx++;
#if FPGA_ON
        void* fpga_buf;
        fpga_buf = fpga_get_writebuf(task.size, BUF_TYPE_SW);
        if(fpga_buf == NULL) {
            fprintf(stderr, "fpga_get_writebuf return NULL\n");
            exit(1);
        }
        memcpy(fpga_buf, task.buf, task.size);
        fpga_writebuf_submit(fpga_buf, task.size, TYPE_CD);
#else
        
        int out_size;
        buf_info_t result;
#if DUMP_FILE
        collect_task_t* sub_head = (collect_task_t*)(task.buf + sizeof(chaindp_sndhdr_t));    //子头指针
        if(sub_head->read_id == 33) {
            FILE* fp = fopen("33_collect_in.dat", "w+");
            if(fp) {
                int ret = fwrite(task.buf, task.size, 1, fp);
                if(ret != 1) {
                    fprintf(stderr, "fwrite failed\n");
                    exit(1);
                }
                fclose(fp);
            }
        }
#endif
        void* out_buf = fpga_work(task.buf, task.size, &out_size);
#if DUMP_FILE
        sub_head = (collect_task_t*)(out_buf + sizeof(chaindp_sndhdr_t));    //子头指针
        if(sub_head->read_id == 33) {
            FILE* fp = fopen("33_collect_out.dat", "w+");
            if(fp) {
                int ret = fwrite(task.buf, task.size, 1, fp);
                if(ret != 1) {
                    fprintf(stderr, "fwrite failed\n");
                    exit(1);
                }
                fclose(fp);
            }
            exit(1);
        }
#endif
        result.buf = out_buf;
        result.size = out_size;
        while(send_fpga_result(result));
#endif
        free(task.buf);
    }
    return NULL;
}

static buf_info_t results[QUEUE_MAX];
static pthread_mutex_t result_mutex = PTHREAD_MUTEX_INITIALIZER;
static int result_head = 0;
static int result_tail = 0;

static int fpga_recv_result_stop = 1;

void init_fpga_result_array()
{
    memset(results, 0, sizeof(results));
    result_head = 0;
    result_tail = 0;
    fpga_recv_result_stop = 1;
    return;
}

void stop_fpga_recv_thread()
{
    fpga_recv_result_stop = 0;
}

static int fpga_result_is_full()
{
    return ((result_tail + 1) % QUEUE_MAX == result_head);
}
static int fpga_result_is_empty()
{
    return (result_head == result_tail);
}

int send_fpga_result(buf_info_t result)
{
    pthread_mutex_lock(&result_mutex);
    if (fpga_result_is_full()) {
        pthread_mutex_unlock(&result_mutex);
        return -1;
    }
    results[result_tail] = result;
    result_tail = (result_tail + 1) % QUEUE_MAX;
        
    pthread_mutex_unlock(&result_mutex);
    return 0;
}

int get_fpga_result(buf_info_t* result)
{
    pthread_mutex_lock(&result_mutex);
    if (fpga_result_is_empty()) {
        pthread_mutex_unlock(&result_mutex);
        return 1;
    }
    *result = results[result_head];
    result_head = (result_head + 1) % QUEUE_MAX;
    pthread_mutex_unlock(&result_mutex);
    return 0;
}
void* recv_task_thread(void* arg)
{
    while(fpga_recv_result_stop) {
#if FPGA_ON
        buf_info_t result;
        int fpga_len;
        void* fpga_buf;
        fpga_buf = fpga_get_retbuf(&fpga_len, RET_TYPE_CS);
        if(fpga_len == 0) {
            fprintf(stderr, "exit recv fpga thread\n");
            return NULL;
        }
        if(fpga_buf == NULL) {
            fprintf(stderr, "fpga_get_retbuf return NULL\n");
            exit(1);
        }
        result.buf = malloc(fpga_len);
        result.size = fpga_len;
        memcpy(result.buf, fpga_buf, result.size);
        fpga_release_retbuf(fpga_buf);
        while(send_fpga_result(result));
#else
        usleep(500000);
#endif
    }
    return NULL;
}