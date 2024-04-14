#include "gpuError.h"
#include <stdio.h>
#include <stdlib.h>

void hipError(
      hipError_t err,
      const char *file,
      int line)
{
   if (err != hipSuccess) {
       printf("*** HIP Error in %s at line %i\n", file, line );
       printf("%s\n",hipGetErrorString( err ));
       exit( EXIT_FAILURE );
   }
}

void checkError(
      const char* file,
      int line)
{
   hipError_t error = hipGetLastError();
   if (error != hipSuccess) {
      printf("*** HIP Error in %s, line %i\n", file, line);
      printf("%s\n", hipGetErrorString(error));
      exit (EXIT_FAILURE);
   }
}

void promptError(
      const char* condition,
      const char* file,
      int line)
{
   printf("\n***ERROR in %s, line %i !\n", file, line);
   printf("%s check failed!\n\n", condition);
   exit( EXIT_FAILURE );
}