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

#define MAX_ARGS       512
#define MAX_BUFFER_LEN 2048
#define MAX_SEQ        (64*1024)
#define MAX_PENDING    5

#define SERVER_PORT      41139

static void
sig_pipe(int signo)
{
  fprintf(stderr, "Exiting due to SIGPIPE\n");
  exit(1);
}

int
main(int argc, char **argv)
{
  int           i;
  int           eod;
  int           pfd1[2];
  int           pfd2[2];

  int           total;
  int32_t       n;

  pid_t         pid;

  char          buffer[MAX_BUFFER_LEN];
  char          seq[MAX_SEQ];

  int                serv_sock;
  int                clnt_sock;
  unsigned short     serv_port;
  struct sockaddr_in serv_addr;
  struct sockaddr_in clnt_addr;

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

      /* Receive message from client */
      if ((n = recv(clnt_sock, seq, MAX_SEQ-1, 0)) < 0) {
        fprintf(stderr, "%s: recv error %d - %s\n", argv[0], errno, strerror(errno));
        exit(1);
      }
      seq[n] = 0;

      while (n > 0) {
        printf("Received %d\n", n);
        //printf("%s", seq);

        n = strlen(seq) + 1;
        if (write(pfd1[1], seq, n) != n) {
          fprintf(stderr, "%s: write (size %d) error %d - %s\n", argv[0], n, errno, strerror(errno));
          exit(1);
        }

        eod = 0;
        total = 0;
        while (!eod) {
          if ((n = read(pfd2[0], buffer, MAX_BUFFER_LEN-1)) == -1) {
            fprintf(stderr, "%s: read error %d - %s\n", argv[0], errno, strerror(errno));
            exit(1);
          }

          if (n == 0) {
            fprintf(stderr, "%s: child pipe closed\n", argv[0]);
            break;
          }

          buffer[n] = 0;
          total += ++n;

          //printf("%6d\n", n);
          //printf("%s", buffer);

          /* Send the length of the output to the client */
          if (send(clnt_sock, &n, sizeof(n), 0) != sizeof(n)) {
            fprintf(stderr, "%s: send (size %d) error %d - %s\n", argv[0], (int)sizeof(n), errno, strerror(errno));
            exit(1);
          }

          /* Send the output to the client */
          if (send(clnt_sock, buffer, n, 0) != n) {
            fprintf(stderr, "%s: send (size %d) error %d - %s\n", argv[0], n, errno, strerror(errno));
            exit(1);
          }

          /* scan backwards looking for the end of the line */
          i = n - 2;
          while (i > 1 && isspace(buffer[i])) i--;

          /* scan backwards looking for the start of the line */
          while (i > 1 && buffer[i] != '\n' && buffer[i] != '\r') i--;

          eod = (buffer[i+1] == '/' && buffer[i+2] == '/');
        }

        printf ("Report sent: %d bytes\n", total);

        /* Notify the client that we have seen the end of the report */
        n = 0;
        if (send(clnt_sock, &n, sizeof(n), 0) != sizeof(n)) {
          fprintf(stderr, "%s: send (size %d) error %d - %s\n", argv[0], (int)sizeof(n), errno, strerror(errno));
          exit(1);
        }

        /* Receive message from client */
        if ((n = recv(clnt_sock, seq, MAX_SEQ-1, 0)) < 0) {
          fprintf(stderr, "%s: recv error %d - %s\n", argv[0], errno, strerror(errno));
          exit(1);
        }
        seq[n] = 0;
      }

      printf("Closing client %s\n", inet_ntoa(clnt_addr.sin_addr));

      close(clnt_sock);
    }

  } else {

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
    *argp++ = argv[i++];
    *argp++ = "--daemon";
    while (i < argc-1) {
      *argp++ = argv[i++];
    }
    *argp++ = "-";
    *argp++ = argv[i];
    *argp   = NULL;

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

