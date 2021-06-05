#include "defs.h"

!
! These preprocessor macros might become handy, specially when we are
! compiling in debug mode.
!

! Checking error rc
#define CHECK_ERR_CODE(localPet, rc)\
   if (rc .NE. 0) then;\
      write (*, "(A, I4, A, I6, A, A, A, I4)") "Processor ", localPet, \
         " exited line: ", __LINE__, " in ", __FILE__, " with error: ", rc; \
      stop;\
   end if;

! Just throwing error
#define THROW_ERR(localPet, errCode, errStr) \
   if (errCode .NE. 0) then; \
      write (*, "(A, I4, A, I6, A, A, A, A)") "Processor: ", localPet, \
                " exited line: ", __LINE__, " in: ", __FILE__, \
                " with error: ", errStr; \
      stop; \
   end if;

! Throwing error if we are in debug mode
#ifdef DEBUG_MODE
#define THROW_DEBUG_ERR(localPet, errCode, errStr) \
   if (errCode .NE. 0) then; \
      write (*, "(A, I4, A, I6, A, A, A, A)") "Processor: ", localPet, \
                " exited line: ", __LINE__, " in: ", __FILE__, \
                " with error: ", errStr; \
      stop; \
   end if;
#else
#define THROW_DEBUG_ERR(localPet, errCode, lineNumber, fileName, errStr)
#endif

!
! You do not need to change anything below here if you are a user
!
#ifdef DEBUG_MANIAC
#define DEBUG_MODE
#endif

#ifndef DEBUG_MODE
#define RELEASE_MODE
#endif
