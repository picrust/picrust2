/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of ivy_mike.
 * 
 *  ivy_mike is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ivy_mike is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ivy_mike.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __ivy_mike__thread_h
#define __ivy_mike__thread_h
#include "compiler_capabilities.h"


#if IVY_MIKE__USE_CPP11
#include <thread>
#include <mutex>

// implementing ivy_mike::thread has become kind of trivial on a conforming c++11 compiler... well done everyone (for once, this includes microsoft)
// TODO: get rid of the ivy_mike namespace alltogether.

namespace ivy_mike {
	typedef std::thread thread;
	typedef std::mutex mutex;

}

#else

// WARNING: these thread implementations are not very exception safe. use boost/just/std::threads if you are interested in correctness...

#ifdef WIN32
// funny win32 api fact: you can deactivate everything and the thread implementation still works.
// This clearly proves that threads are the only part of the win32 api which does not suck. What does the "32" stand for, btw.?
#include "ivymike/disable_shit.h"
#include <windows.h> // oh noes, the dreadful msdn advise from hell: "Header Winbase.h (include Windows.h)"
#include <process.h>    /* _beginthread, _endthread */
// And I heard a voice in the midst of the four beasts,
// And I looked and behold: a pale horse.
// And his name, that sat on him, was Windows.h.
// And Instant Namespace Pollution followed with him.

//#if defined(WIN32) && !defined(_M_X64)
//#define GODAWFULLY_STUPID_32_BIT_WINDOWS_CALLING_CONVENTION __stdcall
//#else
//#define GODAWFULLY_STUPID_32_BIT_WINDOWS_CALLING_CONVENTION
//#endif


namespace ivy_mike {
// TODO: the windows implementation is completely untested, but seems to work for pw_dist


class thread {
public:
	typedef HANDLE native_handle_type;
	

private:
	native_handle_type m_thread;

	thread( const thread &other );
	const thread& operator=(const thread &other );

	template<typename callable>
    static unsigned int __stdcall call( void *f ) {
        std::auto_ptr<callable> c(static_cast<callable *>(f));
        try {
            (*c)();
        } catch( std::runtime_error x ) {
            std::cerr << "uncaught std::runtime_error in ivy_mike::thread:\n" << x.what() << std::endl;
//             std::cerr << x.what() << std::endl;
            
            throw;
        } catch( std::exception x ) {
            std::cerr << "uncaught std::exception in ivy_mike::thread:\n" << x.what() << std::endl;
//             std::cerr << x.what() << std::endl;
            
            throw;
        }
       
		_endthreadex(0);
        return 0;
    }
    
public:

	thread() : m_thread(NULL) {}

	template<typename callable>
    thread( const callable &c ) {
#if 0
		m_thread = CreateThread( 
            NULL,                   // default security attributes
            0,                      // use default stack size  
            call<callable>,       // thread function name
            new callable(c),          // argument to thread function 
            0,                      // use default creation flags 
            NULL);   // returns the thread identifier 

        if( m_thread == NULL) {
            throw std::runtime_error( "could not create thread: don't know why." );
        } 
#else
		// holy crap, the c runtime on windows is a real mess... CreateThread seems to be kind of unsafe for some strange reason
		m_thread = (HANDLE)_beginthreadex(NULL,                   // default security attributes
            0,                      // use default stack size  
            call<callable>,       // thread function name
            new callable(c),          // argument to thread function 
            0,                      // use default creation flags 
            NULL);   // returns the thread identifier 
		
#endif
    }
    
    ~thread() {
        // this (=nothing) is actually the right thing to do (TM) if I read 30.3.1.3 if the iso c++0x standard corrctly
        if( joinable() ) {
            std::cerr << "ivy_mike::thread warning: destructor of joinable thread called. possible ressource leak\n";
        }
    }
    
    void swap( thread &other ) {
        std::swap( m_thread, other.m_thread );
        //std::swap( m_valid_thread, other.m_valid_thread );
    }
    
    bool joinable() {
        return m_thread != NULL;
    }
    
    void join() {
        
        if( joinable() ) {
			std::cerr << "joining ..." << std::endl;
			WaitForSingleObject( m_thread, INFINITE );
			std::cerr << "done" << std::endl;
            m_thread = NULL;
        }
    }
    
    native_handle_type native_handle() {
        if( !joinable() ) {
            throw std::runtime_error( "native_handle: thread not joinable" );
        }
        return m_thread;
    }
};



class mutex {
    CRITICAL_SECTION m_cs;

public:
    
    mutex() {
		InitializeCriticalSection( &m_cs );
    }
    ~mutex() {
		DeleteCriticalSection( &m_cs );
    }
    
    inline void lock() {
		EnterCriticalSection( &m_cs );
    }
    
    inline void unlock() {
		LeaveCriticalSection( &m_cs );
    }
};



}

#else

#include <pthread.h>
#include <vector>
#include <memory>
#include <stdexcept>
#include <cerrno>
namespace ivy_mike {
class thread {
    
public:
    typedef pthread_t native_handle_type;
private:
    native_handle_type m_thread;
    bool m_valid_thread;
    
    thread( const thread &other );
    const thread &operator=( const thread &other );
    
    template<typename callable>
    static void *call( void *f ) {
        std::auto_ptr<callable> c(static_cast<callable *>(f));
        
         try {
            (*c)();
        } catch( std::runtime_error x ) {
            std::cerr << "uncaught std::runtime_error in ivy_mike::thread: " << x.what() << std::endl;
//             std::cerr << x.what() << std::endl;
            
            throw;
        } catch( std::exception x ) {
            std::cerr << "uncaught std::exception in ivy_mike::thread" << std::endl;
//             std::cerr << x.what() << std::endl;
            
            throw;
        }
       
        return 0;
    }
    
public:
    
    thread() : m_valid_thread( false ) {
        
    }
    
    template<typename callable>
    thread( const callable &c ) {
        // dynamically allocating a copy of callable seems to be the only way to get the object into the 
		// thread start function without specializing the whole thread object...
        int ret = pthread_create( &m_thread, 0, call<callable>, new callable(c) );
        
        if( ret == EAGAIN ) {
            throw std::runtime_error( "could not create thread: resource_unavailable_try_again" );
        } else if( ret != 0 ) {
            throw std::runtime_error( "could not create thread: internal error" );
        } else {
            m_valid_thread = true;
        }
    }
    
    ~thread() {
        // this (=nothing) is actually the right thing to do (TM) if I read 30.3.1.3 if the iso c++0x standard corrctly
        if( joinable() ) {
            std::cerr << "ivy_mike::thread warning: destructor of joinable thread called. possible ressource leak\n";
        }
    }
    
    void swap( thread &other ) {
        std::swap( m_thread, other.m_thread );
        std::swap( m_valid_thread, other.m_valid_thread );
    }
    
    bool joinable() {
        return m_valid_thread;
    }
    
    void join() {
        
        if( joinable() ) {
            void *rv;
            pthread_join(m_thread, &rv );
            m_valid_thread = false;
        }
    }
    
    native_handle_type native_handle() {
        if( !joinable() ) {
            throw std::runtime_error( "native_handle: thread not joinable" );
        }
        return m_thread;
    }
};



class mutex {
    pthread_mutex_t m_mtx;

public:
    
    mutex() {
        pthread_mutex_init( &m_mtx, 0);
    }
    ~mutex() {
        pthread_mutex_destroy(&m_mtx);
    }
    
    inline void lock() {
        pthread_mutex_lock(&m_mtx);
    }
    
    inline void unlock() {
        pthread_mutex_unlock(&m_mtx);
    }
};



}
#endif

#endif

namespace ivy_mike {
class thread_group {
    std::vector<thread *> m_threads;
    
public:
    
    ~thread_group() {
         //std::cout << "thread_group destructor: fallback join:\n";
         for( std::vector<thread *>::iterator it = m_threads.begin(); it != m_threads.end(); ++it ) {
             if( (*it)->joinable() ) {
             
                 std::cout << "WARNING: joinable thread in ~thread_group(). This probably means there was an uncaught exception: " << (*it)->joinable() << "\n";
             }
         }
        
        join_all();
        
        for( std::vector<thread *>::iterator it = m_threads.begin(); it != m_threads.end(); ++it ) {
            delete *it;
        }
        
    }
    
    
    
    template<typename callable>
    void create_thread( const callable &c ) {
        m_threads.push_back(0); // may throw. so pre-allocate before the thread is created
        m_threads.back() = new thread( c );
        
        
    }
    
    
    void join_all() throw() {
        
        try {
            for( std::vector<thread *>::iterator it = m_threads.begin(); it != m_threads.end(); ++it ) {
                if( (*it)->joinable() ) {
                    (*it)->join();
                }
                
                //delete (*it);
            }
            
           // m_threads.resize(0);
        } catch(...) {
            std::cerr << "BUG: unexpected exception in thread_group::join_all\n"; // kind of stupid: printing this message might throw...
        }
    }
    
    size_t size() {
        return m_threads.size();
    }
};

inline void swap( thread &t1, thread &t2 ) {
    t1.swap(t2);
}


template <typename mtx_t>
class lock_guard {
    mtx_t &m_mtx;
    

    lock_guard( const lock_guard & );
    lock_guard &operator=(const lock_guard & );

public:
    lock_guard( mtx_t &mtx ) 
     : m_mtx(mtx)
    {
        m_mtx.lock();
    }
    
    ~lock_guard() {
        m_mtx.unlock();
    }
};
}




#endif
