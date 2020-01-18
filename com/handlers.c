/**
 * @file handlers.c
 *
 * Purpose: bring some C methods into Fortran
 *
 * $Id: handlers.c,v 7.0 2018/02/28 18:32:49 meyer8 Exp $
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if defined(FC_FUNC)
FC_FUNC(fcntl, FCNTL)(int *fildes, int *cmd, int *args) {
#else
fcntl_(int *fildes, int *cmd, int *args) {
#endif
  int i;
  i = fcntl(*fildes, *cmd, *args);
  return i;
}

#if defined(FC_FUNC)
FC_FUNC(rdfile, RDFILE)(int *fildes, char *buf, int *nbyte) {
#else
rdfile_(int *fildes, char *buf, int *nbyte) {
#endif
  int i;
  i = read(*fildes, buf, *nbyte);
  return i;
}

