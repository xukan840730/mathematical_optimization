#ifndef _ASSERT_H_
#define _ASSERT_H_

#define xassert(f) { if (!(f)) { int* p = 0; *p = 0; } }
#define ASSERT(f) { if (!(f)) { int* p = 0; *p = 0; } }


#endif