/* hmmd: hmmer deamon client.
 * 
 * MSF, Tue Aug 10, 2010 [Janelia]
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

#define SERVER_PORT      41139
#define MAX_BUF_LEN      (64*1024)
#define MAX_READ_LEN     2048

static void
sig_int(int signo)
{
  fprintf(stderr, "Exiting due to ctrl-c\n");
  exit(1);
}

static void 
usage(char *pgm)
{
  fprintf(stderr, "Usage: %s [-i addr] [-p port]\n", pgm);
  fprintf(stderr, "    -i addr : ip address running daemon (default: 127.0.0.1)\n");
  fprintf(stderr, "    -p port : port daemon listens on (default: 41139)\n");
  exit(1);
}

int main(int argc, char *argv[])
{
  int              i;
  int              n;
  int              eod;

  int32_t          len;
  int32_t          total;

  char             seq[MAX_BUF_LEN];
  char             report[MAX_BUF_LEN];
  char             buffer[MAX_READ_LEN];

  int                 sock;
  char                serv_ip[64];
  unsigned short      serv_port;
  struct sockaddr_in  serv_addr;

  /* set up defaults */
  strcpy(serv_ip, "127.0.0.1");
  serv_port = SERVER_PORT;

  i = 1;
  while (i < argc) {
    if (i + 1 >= argc) usage(argv[0]);
    if (argv[i][0] != '-') usage(argv[0]);
    if (argv[i][1] == 0 || argv[i][2] != 0) usage(argv[0]);
    switch (argv[i][1]) {
    case 'p':
      serv_port = atoi(argv[i+1]);
      ++i;
      break;
    case 'i':
      strcpy(serv_ip, argv[i+1]);
      ++i;
      break;
    default:
      usage(argv[0]);
    }
    ++i;
  }

  /* set up a signal handler for broken pipes */
  if (signal(SIGINT, sig_int) == SIG_ERR) {
    fprintf(stderr, "%s: signal error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  /* Create a reliable, stream socket using TCP */
  if ((sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP)) < 0) {
    fprintf(stderr, "%s: socket error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  /* Construct the server address structure */
  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family      = AF_INET;
  serv_addr.sin_addr.s_addr = inet_addr(serv_ip);
  serv_addr.sin_port        = htons(serv_port);

  /* Establish the connection to the echo server */
  if (connect(sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0) {
    fprintf(stderr, "%s: connect error %d - %s\n", argv[0], errno, strerror(errno));
    exit(1);
  }

  seq[0] = 0;
  while (strncmp(seq, "//", 2) != 0) {

    eod = 0;
    seq[0] = 0;
    fprintf(stdout, "\n\nEnter next sequence:\n");
    while (!eod) {
      if (fgets(buffer, MAX_READ_LEN, stdin) != NULL) {
        strcat(seq, buffer);
        eod = (strncmp(buffer, "//", 2) == 0);
      } else {
        eod = 1;
      }
    }

    /* send the sequence to the daemon and read the results */
    if (strncmp(seq, "//", 2) != 0) {
      
      n = strlen(seq) + 1;

      printf ("Sending data %d\n", n);

      /* Send the string to the server */
      if (send(sock, seq, n, 0) != n) {
        fprintf(stderr, "%s: send (size %d) error %d - %s\n", argv[0], n, errno, strerror(errno));
        exit(1);
      }

      if ((n = recv(sock, &len, sizeof(len), 0)) != sizeof(len)) {
        fprintf(stderr, "%s: recv error %d - %s\n", argv[0], errno, strerror(errno));
        exit(1);
      }

      total = 0;
      while (len != 0) {

        total += len;
        //printf ("received: %d\n", len);

        while (len > 0) {
          if ((n = recv(sock, report, len, 0)) <= 0) {
            fprintf(stderr, "%s: recv report error %d - %s\n", argv[0], errno, strerror(errno));
            exit(1);
          }

          report[n] = 0;
          if (fputs(report, stdout) == EOF) {
            fprintf(stderr, "%s: fputs error %d - %s\n", argv[0], errno, strerror(errno));
            exit(1);
          }
          
          len -= n;
        }

        if ((n = recv(sock, &len, sizeof(len), 0)) <= 0) {
          fprintf(stderr, "%s: recv len error %d - %s\n", argv[0], errno, strerror(errno));
          exit(1);
        }
      }

      //printf("Total: %d\n", total);
    }
  }

  close(sock);
  return 0;
}
