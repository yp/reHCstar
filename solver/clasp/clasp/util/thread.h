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

#ifndef CLASP_UTIL_THREAD_H_INCLUDED
#define CLASP_UTIL_THREAD_H_INCLUDED

#ifndef DISABLE_MULTI_THREADING
#if _WIN32||_WIN64
#define WIN32_LEAN_AND_MEAN // exclude APIs such as Cryptography, DDE, RPC, Shell, and Windows Sockets.
#define NOMINMAX            // do not let windows.h define macros min and max
#endif
#include <tbb/compat/thread> // replace with std::thread once available
#include <tbb/compat/condition_variable>
#include <tbb/mutex.h>
#include <tbb/spin_mutex.h>
#else

namespace tbb {
class NullMutex {   
public:   
	class scoped_lock {   
	public:   
		scoped_lock() {}
		scoped_lock( NullMutex&, bool=false ) {}   
		~scoped_lock() {}
		void acquire( NullMutex& ) {}
		bool upgrade_to_writer()   { return true; }
    bool downgrade_to_reader() { return true; }
		bool try_acquire( NullMutex&, bool = true ) { return true; }
		void release() {}
	};
	NullMutex()     {}
	void lock()     {}
	bool try_lock() { return true; }
	void unlock()   {}
	static const bool is_rw_mutex = false;   
	static const bool is_recursive_mutex = true;
	static const bool is_fair_mutex = true;
private:
	NullMutex(const NullMutex&);   
	NullMutex& operator=(const NullMutex&);   
};

typedef NullMutex mutex; 
typedef NullMutex spin_mutex;

}

#endif

#endif
