/*
 * Copyright (C) 2009-2012 Simon A. Berger
 * 
 * This file is part of papara.
 * 
 *  papara is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  papara is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with papara.  If not, see <http://www.gnu.org/licenses/>.
 */

// FIXME: someone tries to sabotage us on windows by sneaking in the awful min/max macros. Strangely it only happens in this source file... 



#include "ivymike/disable_shit.h" 
#include <iostream>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>


#include "ivymike/fasta.h"
#include "stepwise_align.h"


typedef std::vector<uint8_t> primer_sequence;

std::vector<primer_sequence> load_primers( const std::string &primer_name ) {
    std::vector<primer_sequence> pl;
    
    if( primer_name.empty() ) {
        return pl;
    }
    
    std::ifstream is( primer_name.c_str() );
    
    if( !is.good() ) {
        throw std::runtime_error( "cannot open primer file" );
    }
    
    std::vector<std::string> names;
    
    ivy_mike::read_fasta( is, names, pl );
    
    return pl;
    
    
}


std::vector<uint8_t> to_vec( const std::string &s ) {
    return std::vector<uint8_t>( s.begin(), s.end() );
}

std::vector<primer_sequence> builtin_primers() {
    std::vector<primer_sequence> pl;
    
    pl.push_back( to_vec( "CTTGGTCATTTAGAGGAAGT" ) );
    pl.push_back( to_vec( "CGATGAAGAACGCAG" ) );
    
    return pl;
}

void process_fasta( std::istream &is, const std::string &primer_name )
{

    std::vector<primer_sequence> primers = builtin_primers();
//     std::vector<std::string> primers = load_primers( primer_name );
//     primers.push_back("CTTGGTCATTTAGAGGAAGT");
//     primers.push_back("CGATGAAGAACGCAG");
    
    std::vector<std::string> names;
    std::vector<std::vector<uint8_t> > data;
    
    ivy_mike::read_fasta( is, names, data, false );
    
    
    assert( names.size() == data.size());
    
    std::vector<size_t> rind;
    for( size_t i = 0; i < names.size(); ++i ) {
        rind.push_back(i);
    }
    
    std::random_shuffle( rind.begin(), rind.end() );
    
    const size_t target_size = 10;
    
    size_t num_written = 0;
    
    
    int match_score = 5;
    ivy_mike::scoring_matrix sm(5, -3);
    
    
    while( num_written < target_size && !rind.empty()) {
        size_t r = rind.back();
        rind.pop_back();
        
        
        
        const std::vector<uint8_t> &seq = data.at(r);
        
        
        
        size_t num_rejects = 0;
        for( size_t j = 0; j < primers.size(); ++j ) {
//             std::vector<uint8_t> p( primers[j].begin(), primers[j].end() ); // align_freeshift overwrites the input sequences (i.e., in/out parameters)
            const std::vector<uint8_t> &p = primers[j];
            float score = align_freeshift_score( sm, seq, p, -5, -2, false );
            float escore = match_score * p.size();
            if( score < 0.9 * escore ) {
                ++num_rejects;
            }
            
//             std::cout << "score " << j << ": " << score << " " << p.size() * match_score << "\n";
        }
        
        if( num_rejects != 0 ) {
//             std::cout << "reject!\n";
            continue;
        }
        
        std::cout << ">" << names.at(r) << "\n";
        
        std::copy( data.at(r).begin(), data.at(r).end(), std::ostream_iterator<char>(std::cout) );
        std::cout << "\n";        
        
        ++num_written;
    }
    if( num_written == 0 ) {
        throw std::runtime_error( "none of the input sequences were selected. bailing out." );
    }
    
}


void process_fasta_weighted( std::istream &is, const std::string &primer_name )
{

    std::vector<primer_sequence> primers = builtin_primers();
    //std::vector<primer_sequence> primers = load_primers( primer_name );
//     primers.push_back("CTTGGTCATTTAGAGGAAGT");
//     primers.push_back("CGATGAAGAACGCAG");
    
    std::vector<std::string> names;
    std::vector<std::vector<uint8_t> > data;
    
    ivy_mike::read_fasta( is, names, data, false );
    
    int baseweight = 1000;
    std::vector<int> weights( data.size(), baseweight );
    
    const size_t n = names.size();
    assert( n == data.size());
    
    std::vector<size_t> rind;
    for( size_t i = 0; i < names.size(); ++i ) {
        rind.push_back(i);
    }
    
    std::random_shuffle( rind.begin(), rind.end() );
    
    const size_t target_size = 10;
    
    size_t num_written = 0;
    
    
    int match_score = 5;
    ivy_mike::scoring_matrix sm(5, -3);
    
    
    for( size_t i = 0; i < n; ++i ) {
        // calculate sequence weighting based on primer matches
        
        const std::vector<uint8_t> &seq = data.at(i);
        
        for( size_t j = 0; j < primers.size(); ++j ) {
            const std::vector<uint8_t> &p = primers[j];
            float score = align_freeshift_score( sm, seq, p, -5, -2, false );
            float escore = match_score * p.size();
            // use 'fixed point' arithmetics for the weigths
            weights[i] += floor((score / escore) * baseweight * 10);
        }
        
//         std::cout << "weight: " << i << " " << weights[i] << "\n";
    }
    
    
    while( num_written < target_size && !data.empty() ) {
        // random selection: each remaining sequence forms a bin/interval sized according to the weights.
        // The weigths are integers (approximating fixed-point weights), to prevent float/roundoff weirdness.
        
        int weight_sum = std::accumulate( weights.begin(), weights.end(), 0 );
        
        int r = rand() % weight_sum;
        
        int cur_end = 0;
        size_t sel = size_t(-1);
        size_t n_remain = names.size();
        
        
        // select the bin into which the random value r falls
        for( sel = 0; sel < n_remain; ++sel ) {
            cur_end += weights[sel];
            
            if( cur_end > r ) {
                break;
            }
            
        }
        
        if( sel >= n_remain ) {
            throw std::runtime_error( "bad math error (i.e., I'm too stupid to get simple integer math right)" );
        }
        
        std::cerr << "============== select " << names.at(sel) << " weight " << weights.at(sel) << std::endl;
        std::cout << ">" << names.at(sel) << "\n";
        
        std::copy( data.at(sel).begin(), data.at(sel).end(), std::ostream_iterator<char>(std::cout) );
        std::cout << std::endl;        
        
        
        // erase selected sequence (i.e, random sampling without replacement)
        names.erase( names.begin() + sel );
        data.erase( data.begin() + sel );
        weights.erase( weights.begin() + sel );
        
        ++num_written;
        
        
    }
        
        
    
    
//     while( num_written < target_size && !rind.empty()) {
//         size_t r = rind.back();
//         rind.pop_back();
//         
//         
//         
//         std::vector<uint8_t> seq = data.at(r);
//         
//         
//         
//         size_t num_rejects = 0;
//         for( size_t j = 0; j < primers.size(); ++j ) {
//             std::vector<uint8_t> p( primers[j].begin(), primers[j].end() ); // align_freeshift overwrites the input sequences (i.e., in/out parameters)
//             //const std::vector<uint8_t> &p = primers[j];
//             float score = align_freeshift( sm, seq, p, -5, -2, false );
//             float escore = match_score * p.size();
//             if( score < 0.9 * escore ) {
//                 ++num_rejects;
//             }
//             
// //             std::cout << "score " << j << ": " << score << " " << p.size() * match_score << "\n";
//         }
//         
//         if( num_rejects != 0 ) {
// //             std::cout << "reject!\n";
//             continue;
//         }
//         
//         std::cout << ">" << names.at(r) << "\n";
//         
//         std::copy( data.at(r).begin(), data.at(r).end(), std::ostream_iterator<char>(std::cout) );
//         std::cout << "\n";        
//         
//         ++num_written;
//     }
    if( num_written == 0 ) {
        throw std::runtime_error( "none of the input sequences were selected. bailing out." );
    }
    
}


int main( int argc, char *argv[] ) {
 
    
    if( argc < 2 ) {
        std::cerr << "missing random seed/inout files" << std::endl;
		return -1;
    }
    
    unsigned int seed = atoi( argv[1] );
    std::cerr << "rand seed: " << seed << std::endl;
    
    std::string primer_name;
    if( argc > 3 ) {
        primer_name = argv[3];
    }
    
    srand(seed);
    {
        std::ifstream is( argv[2] );
        
        if( !is.good() ) {
            std::cerr << "WARNING: cannot read: " << argv[2] << "\n";
            
        }
        
        //process_fasta_weighted( is, primer_name );
        process_fasta( is, primer_name );
    }    
//     if( argc == 2 ) {
//         process_fasta_weighted( std::cin, primer_name );
//     } else {
//         for( int i = 2; i < argc; ++i ) {
//             std::ifstream is( argv[i] );
//             
//             if( !is.good() ) {
//                 std::cerr << "WARNING: cannot read: " << argv[i] << "\n";
//                 continue;
//             }
//             
//             process_fasta_weighted( is, primer_name );
//         }
//     }
}

