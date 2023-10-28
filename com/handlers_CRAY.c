FCNTL (int *fildes,int *cmd,int *args) {
  int i;
  i = fcntl(*fildes,*cmd,*args);
  return i;
}

RDFILE (int *fildes, char *buf, int *nbyte) {
  int i;
  i = read(*fildes, buf, *nbyte);
  return i;
}
