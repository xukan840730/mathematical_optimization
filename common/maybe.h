#ifndef _MAYBE_H_
#define _MAYBE_H_

#include "assert.h"

template<typename T1>
struct Maybe
{
	bool valid;
	T1 data;	

	Maybe() { valid = false; }
	Maybe(const T1& d) { data = d; valid = true; }

	T1& get() { ASSERT(valid); return data; }
	const T1& get() const { ASSERT(valid); return data; }
};

#endif