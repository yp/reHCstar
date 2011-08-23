// 
// Copyright (c) 2010, Benjamin Kaufmann
// 
// This file is part of Clasp. See http://www.cs.uni-potsdam.de/clasp/ 
// 
// Clasp is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// Clasp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Clasp; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

#ifndef CLASP_UTIL_ATOMIC_H_INCLUDED
#define CLASP_UTIL_ATOMIC_H_INCLUDED
#ifdef _MSC_VER
#pragma once
#endif

#ifndef DISABLE_MULTI_THREADING
#include <tbb/atomic.h>
namespace std {
	using tbb::atomic;
}
#else
namespace std {
template <class T>
struct atomic {
	typedef T value_type;
	atomic() : value(value_type()) {}
	atomic& operator=(value_type t) { value = t; return *this; }
	operator value_type() const { return value; }
	value_type operator+=(value_type v) { return value += v; }
	value_type operator-=(value_type v) { return value -= v; }
	value_type operator++()             { return ++value;    }
	value_type operator--()             { return --value;    }
	value_type operator->() const       { return value;    }
	value_type compare_and_swap(value_type y, value_type z) {
		if (value == z) {
			value = y;
			return z;
		}
		return value;
	}
	T value;
};
}
#endif

#endif
