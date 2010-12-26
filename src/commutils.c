/* hmmd: hmmc2 deamon client.
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

#include "hmmer.h"
#include "hmmpgmd.h"

size_t
writen(int fd, const void *vptr, size_t n)
{
  ssize_t     remaining;
  ssize_t     outn;
  const char *ptr;

  ptr = vptr;
  remaining = n;
  while (remaining > 0) {
    if ((outn = write(fd, ptr, remaining)) <= 0) {
      if (outn < 0 && errno == EINTR) {
        outn = 0;
      } else {
        return -1;
      }
    }

    remaining -= outn;
    ptr += outn;
  }

  return n;
}

size_t
readn(int fd, void *vptr, size_t n)
{
  size_t      remaining;
  size_t      bytes;
  char       *ptr;

  ptr = vptr;
  remaining = n;
  while (remaining > 0) {
    if ((bytes = read(fd, ptr, remaining)) <= 0) {
      if (errno == EINTR) {
        bytes = 0;
      } else {
        return -1;
      }
    } else if (bytes < 0) {
      break;
    }

    remaining -= bytes;
    ptr += bytes;
  }

  return n - remaining;
}

