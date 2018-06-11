
#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <boost/numeric/ublas/matrix.hpp>
#include "raxml_interface.h"



int main( int argc, char * argv[] ) {
    
    if( argc != 2 ) {
        std::cerr << "missing parameters: file name\n";
        return 1;
    }
    char *file_name = argv[1];
    
    std::ifstream is( file_name, std::ios::binary );
    
    auto pvecs = read_binary_anc_probs(is);
    
    for( auto & pv : pvecs )
//     auto &pv = pvecs.at(7);
    {
//         std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
        
        
        for( auto it = pv.begin2(); it != pv.end2(); ++it ) {
            size_t num = std::count_if(it.begin(), it.end(), [](double a){return a > 0.255;} );
            
            if( !true ) {
                if( num > 1 ) {
                    std::cout << "meeeep\n";
                    
                }
            } else {
                
                //std::cout << num << " ";
                std::cout << std::distance( pv.begin2(), it ) << " ";
                std::for_each( it.begin(), it.end(), [&](double v) { std::cout << std::setw(14) << v; } );
                std::cout << "\n";
            }
        }
        
        
        
//         break; // stop after first vector
    }
    
    
}