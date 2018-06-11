#include <iostream>
#include <algorithm>
#include <iterator>
#include "ivymike/fasta.h"
#include "ivymike/multiple_alignment.h"
#include "ivymike/algorithm.h"
#include "ivymike/tree_traversal_utils.h"

static bool nogap( char c ) {
    c = toupper(c);
    return c != '-' && c != 'N';
}

#if __cplusplus <= 199711L && !defined(_MSC_VER)
namespace std {
// OH NOES, copy_if is c++11... rip it from gnu
 template<typename _InputIterator, typename _OutputIterator,
           typename _Predicate>
_OutputIterator copy_if(_InputIterator __first, _InputIterator __last,
                        _OutputIterator __result, _Predicate __pred)
{
    
    for (; __first != __last; ++__first)
        if (__pred(*__first))
        {
            *__result = *__first;
            ++__result;
        }
        return __result;
}
}
#endif

int main() {
    ivy_mike::multiple_alignment ma;
    ma.load_phylip( std::cin );
    
    for( size_t i = 0; i < ma.names.size(); ++i ) {
        std::cout << ">" << ma.names[i] << "\n";
        
        
        const std::vector< uint8_t > &d = ma.data.at(i);
        
        
        std::copy_if( d.begin(), d.end(), std::ostream_iterator<char>(std::cout), nogap );
        std::cout << "\n";
        
    }
    
    
    return 1;
}