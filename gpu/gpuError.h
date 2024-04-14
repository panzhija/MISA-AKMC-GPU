#ifndef HIP_ERROR_H
#define HIP_ERROR_H

#include <hip/hip_runtime.h>

extern void hipError(
      hipError_t err,
      const char *file,
      int line);

extern void checkError(
      const char *file,
      int line);

extern void promptError(
      const char* condition,
      const char* file,
      int line);

#define HANDLE_HIP( err ) (hipError( err, __FILE__, __LINE__ ))
#define CHECK_FOR_ERROR() (checkError( __FILE__, __LINE__ ))
#define CHECK(expr) \
         if (!(expr)) {promptError(#expr,__FILE__,__LINE__);}

#endif /*HIP_ERROR_H*/