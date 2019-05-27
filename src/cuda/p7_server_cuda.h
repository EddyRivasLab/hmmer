#ifndef __P7SERVERCUDA_H
#define __P7SERVERCUDA_H
#ifdef __cplusplus // magic to make C++ compilers happy
extern "C" {
#endif



void *p7_server_cuda_worker_thread(void *worker_argument);

#ifdef __cplusplus // magic to make C++ compilers happy
}
#endif
#endif