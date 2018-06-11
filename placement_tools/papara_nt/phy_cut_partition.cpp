#include <stdexcept>
#include <iomanip>
#include <iterator>
#include <fstream>

#include "ivymike/large_phylip.h"
#include "blast_partassign.h"
#include "papara.h"

using ivy_mike::large_phylip;

int main( int argc, char *argv[] ) {
    const char *in_phylip = argv[1];
    const char *partfile = argv[2];
    const char *part_name = argv[3];
    
    
    assert( argc == 4 );
    assert( in_phylip != 0 );
    assert( partfile != 0 );
    assert( part_name != 0 );
    
    
    large_phylip lp( in_phylip );
//     lp.map();
    
    
    size_t col_min = -1;
    size_t col_max = -1;
    {
        std::ifstream is( partfile );
        
        
        while( is.good() ) {
            partassign::partition part = partassign::next_partition(is);
            if( part.start == -1 ) {
                throw std::runtime_error( "partition not found" );
            }
            
            std::cout << "part: " << part.start << " "<< part.gene_name << "\n";
            
            assert( part.start >= 0 );
            assert( part.end > part.start );
            
            if( part.gene_name == part_name ) {
                col_min = part.start;
                col_max = part.end + 1;
                break;
            }
        }
    }

    assert( col_min != size_t(-1));
    assert( col_min != col_max );
    
    std::ostream &os = std::cout;
    
    size_t name_width = lp.max_name_len() + 1;
   
    os << lp.size() << " " << col_max - col_min << "\n"; 
    for( int i = 0, size = lp.size(); i < size; ++i ) {
        std::string name = lp.name_at(i);
        
        ivy_mike::u1_t* seq_start = lp.sequence_begin_at(i);
        ivy_mike::u1_t* seq_end = lp.sequence_end_at(i);
        
        size_t seq_len = std::distance( seq_start, seq_end );
        assert( col_min < seq_len );
        assert( col_max <= seq_len );
        
        
        os << std::setw(name_width) << std::left << name;
        std::copy( seq_start + col_min, seq_start + col_max, std::ostream_iterator<char>(os) );
        os << "\n";
    }
    
    
}
