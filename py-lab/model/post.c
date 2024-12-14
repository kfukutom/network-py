#include <stdio.h>
#include "my_c_library.c" // Include the C library

int main() {
  
  int result = my_function(5, 3);
  printf("Result of my_function(5, 3): %d\n", result);

  return 0;
}