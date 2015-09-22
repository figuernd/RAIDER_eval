#ifndef UTILS_H
#define UTILS_H
void die(const char *msg){
  printf("Error, aborting program...\n");
  printf("%s", msg);
  printf("\n");
  exit(-1);
}
#endif
