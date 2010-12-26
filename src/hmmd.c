/* hmmd: hmmer deamon to run hmmer programs.
 * 
 * MSF, Fri Aug 06, 2010 [Janelia]
 * SVN $Id: hmmsearch.c 3324 2010-07-07 19:30:12Z wheelert $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include <ctype.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <pthread.h>


#define MAX_ARGS       512
#define MAX_BUFFER_LEN 2048
#define MAX_SEQ        (512*1024)
#define MAX_PENDING    5

#define SERVER_PORT      41139

pthread_mutex_t  SEARCH_MUTEX;
pthread_cond_t   SEARCH_COND;

int              SEARCH_BUSY;

typedef struct {
  int     pfd1;
  int     pfd2;

  int     sock;

} THREAD_ARGS;

static void
sig_pipe(int signo)
{
  fprintf(stderr, "Exiting due to SIGPIPE\n");
  exit(1);
}

void *
thread_main(void *args)
{
  int           i;
  int           eod;
  int           total;

  int32_t       n;

  int           clnt_sock;

  char         *seq;
  char          buffer[MAX_BUFFER_LEN];

  THREAD_ARGS  *info = (THREAD_ARGS *)args;

  /* Guarantees that thread resources are deallocated upon return */
  pthread_detach(pthread_self()); 

  /* Extract socket file descriptor from argument */
  clnt_sock = info->sock;

  seq = malloc(MAX_SEQ);
  if (seq == NULL) {
    fprintf(stderr, "%08X: malloc seq error\n", (unsigned int)pthread_self());
    exit(1);
  }

  /* Receive message from client */
  if ((n = recv(clnt_sock, seq, MAX_SEQ-1, 0)) < 0) {
    fprintf(stderr, "%08X: recv error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
    exit(1);
  }
  seq[n] = 0;

  while (n > 0) {
    char cache[32];

    printf("%08X: Received %d\n", (unsigned int)pthread_self(), n);
    //printf("%s", seq);

    /* wait for a phmmer to be free */
    if (pthread_mutex_lock (&SEARCH_MUTEX) != 0) {
      fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    while (SEARCH_BUSY) {
      if (pthread_cond_wait (&SEARCH_COND, &SEARCH_MUTEX) != 0) {
        fprintf(stderr, "%08X: cond wait error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
        exit(1);
      }
    }

    SEARCH_BUSY = 1;

    if (pthread_mutex_unlock (&SEARCH_MUTEX) != 0) {
      fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    printf("%08X: Starting search %d\n", (unsigned int)pthread_self(), n);

    n = strlen(seq);
    if (write(info->pfd1, seq, n) != n) {
      fprintf(stderr, "%08X: write (size %d) error %d - %s\n", (unsigned int)pthread_self(), n, errno, strerror(errno));
      exit(1);
    }

    eod = 0;
    total = 0;
    while (!eod) {
      if ((n = read(info->pfd2, buffer, MAX_BUFFER_LEN-1)) == -1) {
        fprintf(stderr, "%08X: read error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
        exit(1);
      }

      if (n == 0) {
        fprintf(stderr, "%08X: child pipe closed\n", (unsigned int)pthread_self());
        break;
      }

      buffer[n] = 0;
      total += n;

      //printf("%6d\n", n);
      //i = n - 30;
      //if (i < 0) i = 0;
      //fprintf(stdout, "%8d: %s\n", total, buffer + i);
      //fflush(stdout);

      /* Send the output to the client */
      if (send(clnt_sock, buffer, n, 0) == -1) {
        fprintf(stderr, "%08X: send (size %d) error %d - %s\n", (unsigned int)pthread_self(), n, errno, strerror(errno));
        exit(1);
      }

      if (n < sizeof(cache)) {
        int s = sizeof(cache)-n;
        memmove(cache, cache+n, s);
        memcpy(cache+s, buffer, n);
      } else {
        memcpy(cache, buffer+n-sizeof(cache), sizeof(cache));
      }

      //for (i = 0; i < sizeof(cache); ++i) fputc(cache[i], stdout);
      //fputc('\n', stdout);
      //fflush(stdout);

      /* scan backwards looking for the end of the line */
      i = sizeof(cache) - 1;
      while (i > 1 && isspace(cache[i])) i--;

      /* scan backwards looking for the start of the line */
      while (i > 1 && cache[i] != '\n' && cache[i] != '\r') i--;
      
      eod = (cache[i+1] == '/' && cache[i+2] == '/');
    }

    printf ("%08X: Report sent: %d bytes\n", (unsigned int)pthread_self(), total);

    /* signal that phmmer is available */
    if (pthread_mutex_lock (&SEARCH_MUTEX) != 0) {
      fprintf(stderr, "%08X: mutex lock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    SEARCH_BUSY = 0;

    if (pthread_mutex_unlock (&SEARCH_MUTEX) != 0) {
      fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    if (pthread_cond_broadcast (&SEARCH_COND) != 0) {
      fprintf(stderr, "%08X: mutex unlock error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }

    /* Receive message from client */
    if ((n = recv(clnt_sock, seq, MAX_SEQ-1, 0)) < 0) {
      fprintf(stderr, "%08X: recv error %d - %s\n", (unsigned int)pthread_self(), errno, strerror(errno));
      exit(1);
    }
    seq[n] = 0;
  }

  close(clnt_sock);

  free(args);
  free(seq);

  return (NULL);
}

int
main(int argc, char **argv)
{
  int           i;
  int           pfd1[2];
  int           pfd2[2];

  pthread_t     threadID;

  int32_t       n;

  pid_t         pid;

  int                serv_sock;
  int                clnt_sock;
  unsigned short     serv_port;
  struct sockaddr_in serv_addr;
  struct sockaddr_in clnt_addr;

  THREAD_ARGS *args;

  /* make sure the agguement list does not overflow */
  if (argc < 3 || (strcmp(argv[1], "-p") == 0 && argc < 5)) {
    fprintf(stderr, "Usage: %s [-p port] hmmpgm [options...] db\n", argv[0]);
    fprintf(stderr, "    port    - port for daemon to listen on\n");
    fprintf(stderr, "    hmmpgm  - hmmer program to run (i.e. phmmer)\n");
    fprintf(stderr, "    options - options to pass to the hmmer program\n");
    fprintf(stderr, "    db      - database to use\n");
    exit(1);
  }

  if (argc >= MAX_ARGS) {
    fprintf(stderr, "%s: too many command line arguements (max %d)\n", argv[0], MAX_ARGS);
    exit(1);
  }

  serv_port = (strcmp(argv[1], "-p") == 0) ? atoi(argv[2]) : SERVER_PORT;

  /* Create socket for incoming connections */
  if ((serv_sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    fprintf(stderr, "%s: socket error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }
      
  /* Construct local address structure */
  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  serv_addr.sin_port = htons(serv_port);

  /* Bind to the local address */
  if (bind(serv_sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    fprintf(stderr, "%s: bind error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  /* Mark the socket so it will listen for incoming connections */
  if (listen(serv_sock, MAX_PENDING) < 0) {
    fprintf(stderr, "%s: listen error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  /* set up a signal handler for broken pipes */
  if (signal(SIGPIPE, sig_pipe) == SIG_ERR) {
    fprintf(stderr, "%s: signal error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  if (pthread_mutex_init(&SEARCH_MUTEX, NULL) != 0) {
    fprintf(stderr, "%s: mutex_init error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  if (pthread_cond_init(&SEARCH_COND, NULL) != 0) {
    fprintf(stderr, "%s: cond_init error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  /* create the pipes for the parent and child to talk over */
  if (pipe(pfd1) != 0 || pipe(pfd2) != 0) {
    fprintf(stderr, "%s: pipe error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  /* fork the hmmer program and set up the communication channels */
  if ((pid = fork()) == -1) {
    fprintf(stderr, "%s: fork error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  if (pid > 0) {

    /* close the unused portions of the parents pipe */
    close(pfd1[0]);
    close(pfd2[1]);

    for ( ;; ) {

      /* Wait for a client to connect */
      n = sizeof(clnt_addr);
      if ((clnt_sock = accept(serv_sock, (struct sockaddr *)&clnt_addr, (unsigned int *)&n)) < 0) {
        fprintf(stderr, "%s: accept error %d - %s\n", argv[0], errno, strerror(errno));
        exit(1);
      }

      printf("Handling client %s\n", inet_ntoa(clnt_addr.sin_addr));

      args = malloc(sizeof(THREAD_ARGS));
      if (args == NULL) {
        fprintf(stderr, "%s: malloc error\n", argv[0]);
        exit(1);
      }
      args->pfd1 = pfd1[1];
      args->pfd2 = pfd2[0];
      args->sock = clnt_sock;

      if (pthread_create(&threadID, NULL, thread_main, (void *)args) != 0) {
        fprintf(stderr, "%s: pthread_create error %d - %s\n", argv[0], errno, strerror(errno));
        exit(1);
      }
    }

  } else {

    int len;
    char *pgm;

    char **argp;
    char *arglist[MAX_ARGS];
    //char *arglist[] = { "./src/phmmer", "--daemon", "-", "../../sequences/sprot57", NULL };

    /* close the unused portions of the childs pipe */
    close(pfd1[1]);
    close(pfd2[0]);

    /* replace the childs stdin using the pipes input stream */
    if (pfd1[0] != STDIN_FILENO) {
      if (dup2(pfd1[0], STDIN_FILENO) != STDIN_FILENO) {
	fprintf(stderr, "%s: dup2 error on stdin %d - %s\n", argv[0], errno, strerror(errno));
	exit(1);
      }
      close(pfd1[0]);
    }

    /* replace the childs stdout using the pipes output stream */
    if (pfd2[1] != STDOUT_FILENO) {
      if (dup2(pfd2[1], STDOUT_FILENO) != STDOUT_FILENO) {
    	fprintf(stderr, "%s: dup2 error on stdin %d - %s\n", argv[0], errno, strerror(errno));
    	exit(1);
      }
      close(pfd2[1]);
    }

    /* process the argument list */
    argp = arglist;
    i = (strcmp(argv[1], "-p") == 0) ? 3 : 1;
    *argp++ = pgm = argv[i++];

    len = strlen(pgm);
    len -= strlen("hmmscan");
    if (len < 0) len = 0;

    if (strcmp(pgm+len, "hmmpgmd") != 0) *argp++ = "--daemon";
    while (i < argc-1) {
      *argp++ = argv[i++];
    }

    /* if we are running hmmscan, the '-' is the last arguement */
    len = strlen(pgm);
    len -= strlen("hmmscan");
    if (len < 0) len = 0;
    if (strcmp(pgm+len, "hmmpgmd") == 0) {
      *argp++ = argv[i];
      *argp   = NULL;
    } else if (strcmp(pgm+len, "hmmscan") == 0) {
      *argp++ = argv[i];
      *argp++ = "-";
      *argp   = NULL;
    } else {
      *argp++ = "-";
      *argp++ = argv[i];
      *argp   = NULL;
    }
    
    printf("Cmdline: ");
    argp = arglist;
    while (*argp) printf("%s ", *argp++);
    printf("\n");

    /* start the hmmer program */
    if (execv(arglist[0], arglist) == -1) {
      fprintf(stderr, "%s: execv error on %s %d - %s\n", argv[0], arglist[0], errno, strerror(errno));
      exit(1);
    }
  }

  return 0;
}


/*****************************************************************
 * @LICENSE@
 *****************************************************************/

