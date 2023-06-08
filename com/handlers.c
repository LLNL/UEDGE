/**
 * @file handlers.c
 *
 * Purpose: bring some C methods into Fortran
 *
 * $Id: handlers.c,v 7.1 2021/03/05 21:58:20 meyer8 Exp $
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>
#include <fcntl.h>

#if defined(FC_FUNC)
int FC_FUNC(fcntl, FCNTL)(int *fildes, int *cmd, int *args) {
#else
int fcntl_(int *fildes, int *cmd, int *args) {
#endif
  int i;
  i = fcntl(*fildes, *cmd, *args);
  return i;
}

#if defined(FC_FUNC)
int FC_FUNC(rdfile, RDFILE)(int *fildes, char *buf, int *nbyte) {
#else
int rdfile_(int *fildes, char *buf, int *nbyte) {
#endif
  int i;
  i = read(*fildes, buf, *nbyte);
  return i;
}

