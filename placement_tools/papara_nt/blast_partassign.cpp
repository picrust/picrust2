#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <functional>
#include <algorithm>
#include <cctype>
#include <cassert>
#include <map>

#include "blast_partassign.h"
#include "papara.h"

using partassign::partition;
using partassign::blast_hit;
using papara::queries;
using papara::references;

blast_hit partassign::next_hit( std::istream &is ) {
    std::string line;
    
    assert( is.good() );
    
    std::getline( is, line );
    

    // TODO: ok returning an default constructed blast_hit to signal EOF is not the best idea... I actually planned that EOF 
    // should be checked _before_ calling this function, but it didn't turn out to work...
    
    if( line.empty() ) {
        return blast_hit();
    }
    
//     for( size_t i = 0; i < 12; ++i ) {
//         int s;
//         ss >> s;
//         std::cout << i << " " << s << "\n";
//         
//         
//     }
    
    std::string qs, ref;
    float ident;
    int len, mismatch, gap_open;
    
    std::string evalue;
    blast_hit hit;
    
    //bool valid = ss >> hit.qs_name >> hit.ref_name >> ident >> len >> mismatch >> gap_open >> hit.qs_start >> hit.qs_end >> hit.ref_start >> hit.ref_end >> evalue >> hit.bit_score;
    std::istringstream ss(line);
    bool valid = !(ss >> hit.qs_name >> hit.ref_name >> ident >> len >> mismatch >> gap_open >> hit.qs_start >> hit.qs_end >> hit.ref_start >> hit.ref_end >> evalue >> hit.bit_score).fail();

    // convert indices into proper form
    hit.ref_start -= 1;
    hit.ref_end -= 1;
    hit.qs_start -= 1;
    hit.qs_end -= 1;
    
    if( !valid ) {
        std::cerr << "bad line: " << line << "\n";
        throw std::runtime_error( "could not read line in blast output file.\n" );
    }
    
    return hit;
    
}

template<typename T>
T from_string( const std::string &str ) {
    std::istringstream ss(str);
    T x;
    ss >> x;
    return x;
}
// static bool not_space( char x ) {
//     return !isspace(x); // doing the same with pre-c++11 functional is ridiculus
// }

template<typename iiter>
static void part_check_error( iiter it, iiter end, const std::string &line ) {
    if( it == end ) {
        std::stringstream ss;
        ss << "error parsing line of partition file:\n" << line; 
        
        throw std::runtime_error( ss.str() );
    }
    
}

partition partassign::next_partition( std::istream &is ) {
    // TODO: there is absolutely no (=ZILCH) error checking in this function
    
    std::string line;
    
    assert( is.good() );
    
    std::getline( is, line );
    
    std::istringstream ss(line);

    // TODO: d.t.o.
    if( line.empty() ) {
        return partition();
    }
    
    
    // parse model name:
    // it = start if line (maybe scan for non-space?)
    // it_next = first ','
    std::string::iterator it = line.begin();
    std::string::iterator it_next = std::find( it, line.end(), ',' );
    
    std::string model_name( it, it_next );
//     std::cout << "model: '" << model_name << "'\n";
    
    it = it_next + 1;
    
    
//     std::find_if( it, line.end(), not_space );
    
    // parse gene name:
    // it = first non_space character after the ','
    // it_next = first space after the gene name
    it = std::find_if( it, line.end(), std::not1(std::ptr_fun<int,int>(std::isspace) ) ); // here we go, scott meyers would be proud of me ;)
    part_check_error( it, line.end(), line );
    
    it_next = std::find_if( it, line.end(), std::ptr_fun<int,int>(std::isspace) );
    part_check_error( it_next, line.end(), line );
    
    std::string gene_name( it, it_next );
//     std::cout << "gene: '" << gene_name << "'\n";
    
    it = std::find( it_next, line.end(), '=' );
    part_check_error( it, line.end(), line );
    
    ++it;
    part_check_error( it, line.end(), line );
    
    it = std::find_if( it, line.end(), std::not1(std::ptr_fun<int,int>(std::isspace) ) );
    part_check_error( it, line.end(), line );
    
    it_next = std::find( it, line.end(), '-' );
    part_check_error( it, line.end(), line );
    
    std::string start_str( it, it_next );
    it = it_next + 1;
    std::string end_str( it, line.end() );
//     std::cout << "start: '" << start_str << "'\n";
//     std::cout << "end: '" << end_str << "'\n";
//     
    partition part;
    part.start = from_string<int>(start_str) - 1;
    part.end = from_string<int>(end_str) - 1;
    part.gene_name = gene_name;
    return part;
}

// int main() {
//     
//     {
//         std::ifstream is( "test.model" );
//         assert( is.good() );
//         
//         while( is.good() ) {
//             partition p = partassign::next_partition(is);
//             
//             if( p.start == -1 ) {
//                 break;
//             }
//             
//             std::cout << "part: " << p.start << " " << p.end << "\n";
//             
//             
//         }
//         
//     }
//     
//     {
//         std::ifstream is( "test.out" );
//         assert( is.good() );
//         
//         std::map<std::string,blast_hit> hit_map;
//         
//         while( is.good() ) {
//             
//             blast_hit hit = partassign::next_hit( is );
//             
//             if( hit.qs_name.empty() ) {
//                 break;
//             }
//             
//             //std::cout << "hit: " << hit.qs_start << " " << hit.bit_score << "\n";
//             
//             
//             // check for multiple hits per QS, keep the hit with the highest bit-score
//             std::map< std::string, blast_hit >::iterator it = hit_map.lower_bound( hit.qs_name );
//             if( it != hit_map.end() && it->second.qs_name == hit.qs_name ) {
//                 
//                 if( it->second.bit_score < hit.bit_score ) {
//                     it->second = hit;
//                     //                 std::cout << "replace\n";
//                 }
//             } else {
//                 hit_map.insert( it, std::make_pair( hit.qs_name, hit ));
//             }
//         }
//     }
//     
// }
namespace partassign
{
part_assignment::part_assignment ( std::istream& blast_out, std::istream& part_file )
{
    while ( part_file.good() ) {
        partitions_.push_back ( partassign::next_partition ( part_file ) );
        if ( partitions_.back().start == -1 ) { // returning a partition with negaitve indices is next_partition's way of signalling EOF
            partitions_.pop_back();
            break;
        }
    }


    std::map<std::string,partassign::blast_hit> hit_map;

    while ( blast_out.good() ) {

        blast_hit hit = partassign::next_hit ( blast_out );

        if ( hit.qs_name.empty() ) {
            break;
        }

        //std::cout << "hit: " << hit.qs_start << " " << hit.bit_score << "\n";


        // check for multiple hits per QS, keep the hit with the highest bit-score
        std::map< std::string, blast_hit >::iterator it = hit_map.lower_bound ( hit.qs_name );
        if ( it != hit_map.end() && it->second.qs_name == hit.qs_name ) {

            if ( it->second.bit_score < hit.bit_score ) {
                it->second = hit;
                //                 std::cout << "replace\n";
            }
        } else {
            hit_map.insert ( it, std::make_pair ( hit.qs_name, hit ) );
        }
    }
    hits_.swap(hit_map);

//     for ( std::map<std::string,blast_hit>::iterator it = hit_map.begin(), eit = hit_map.end(); it != eit; ++it ) {
//         int start = it->second.ref_start;
//         int end = it->second.ref_end;
// 
//         int part_idx = -1;
// 
//         for ( size_t i = 0; i < partitions_.size(); ++i ) {
//             const partassign::partition &part = partitions_[i];
// 
//             if ( start >= part.start && end <= part.end ) {
//                 part_idx = int ( i );
//                 break;
//             }
//         }
// 
//         if ( part_idx == -1 ) {
//             std::cerr << "QS cannot be uniqely assigned to a single partition: " << it->first << "[" << start << "-" << end << "\n";
//             throw std::runtime_error ( "partitons incompatible with balst hits" );
//         }

//         a [it->first] = part_idx;
// 
//     }

}
const blast_hit& part_assignment::get_blast_hit ( const std::string& qs_name ) const 
{
    std::map< std::string, partassign::blast_hit >::const_iterator it = hits_.find ( qs_name );
    if ( it == hits_.end() ) {
        throw std::runtime_error ( "qs name not found" );
    }

    return it->second;
}
// const partition& part_assignment::partition ( const std::string& name ) const
// {
//     std::map< std::string, int >::const_iterator it = assignments_.find ( name );
// 
//     if ( it == assignments_.end() ) {
//         std::stringstream ss;
//         ss << "no partition assignment for qs " << name;
//         throw std::runtime_error ( ss.str() );
//     }
// 
//     return partitions_.at ( it->second );
// }


template<typename pvec_t, typename seq_tag>
std::vector<std::pair<size_t,size_t> > resolve_qs_bounds( references<pvec_t,seq_tag> &refs, queries<seq_tag> &qs, const partassign::part_assignment &part_assign ) {
    std::vector<std::pair<size_t,size_t> > bounds;
    
    for( size_t i = 0; i < qs.size(); ++i ) {
        const std::string &qs_name = qs.name_at(i);
        const partassign::blast_hit &hit = part_assign.get_blast_hit( qs_name );
        
         size_t ref_idx = refs.find_name( hit.ref_name );
         
         if( ref_idx == size_t(-1) ) {
             throw std::runtime_error( "ref name of blast hit not found" );
         }
         const std::vector<int> &ng_map = refs.ng_map_at(ref_idx);
         
         if( size_t(hit.ref_start) >= ng_map.size() || size_t(hit.ref_end) >= ng_map.size() ) {
             std::cerr << hit.ref_start << " " << hit.ref_end << " " << ng_map.size() << "\n";
             throw std::runtime_error( "blast hit region outside of reference sequence" );
         }
         
         // map position in (non-gappy) ref sequence onto alignment column
         int col_start = ng_map.at(hit.ref_start);
         int col_end = ng_map.at(hit.ref_end);
         
         int part_idx = -1;
         
         std::vector< partassign::partition > partitions = part_assign.partitions();
         for ( size_t i = 0; i < partitions.size(); ++i ) {
             const partassign::partition &part = partitions[i];
             
             if ( col_start >= part.start && col_end <= part.end ) {
                 part_idx = int ( i );
                 break;
             }
         }
         
         std::cout << "qs part: " << qs_name << " " << part_idx << "\n";
        
        
         if ( part_idx == -1 ) {
             std::cerr << "QS cannot be uniquely assigned to a single partition: " << qs_name << " [" << col_start << "-" << col_end << "]\n";
           //  throw std::runtime_error ( "partitons incompatible with blast hits" );
             
             std::cerr << "falling back to full region\n";
             bounds.push_back( std::make_pair( -1, -1 ));
         } else {
         
             bounds.push_back( std::make_pair( partitions[part_idx].start, partitions[part_idx].end ));
         }
    }

    
    return bounds;
}


std::pair<size_t,size_t> partition_bounds( std::istream &is, const std::string &name ) {
    std::pair<size_t,size_t> bounds(-1,-1);
    
    while ( is.good() ) {
     
        partition p = partassign::next_partition ( is );
        if ( p.start == -1 ) { // returning a partition with negaitve indices is next_partition's way of signalling EOF
            break;
        }
        
//         std::cout << "gene: " << p.gene_name << "\n";
        
        if( p.gene_name == name ) {
            return std::pair<size_t,size_t>(p.start, p.end);
        }
        
    }
    return std::pair<size_t,size_t>(-1,-1);
    
}

// combinatorial explosion hazard ahead... if another function comes along put it into a driver class just like papara::driver.

template std::vector<std::pair<size_t,size_t> > resolve_qs_bounds<pvec_cgap,sequence_model::tag_aa>( references<pvec_cgap,sequence_model::tag_aa> &refs, queries<sequence_model::tag_aa> &qs, const partassign::part_assignment &part_assign );
template std::vector<std::pair<size_t,size_t> > resolve_qs_bounds<pvec_cgap,sequence_model::tag_dna>( references<pvec_cgap,sequence_model::tag_dna> &refs, queries<sequence_model::tag_dna> &qs, const partassign::part_assignment &part_assign );
template std::vector<std::pair<size_t,size_t> > resolve_qs_bounds<pvec_pgap,sequence_model::tag_aa>( references<pvec_pgap,sequence_model::tag_aa> &refs, queries<sequence_model::tag_aa> &qs, const partassign::part_assignment &part_assign );
template std::vector<std::pair<size_t,size_t> > resolve_qs_bounds<pvec_pgap,sequence_model::tag_dna>( references<pvec_pgap,sequence_model::tag_dna> &refs, queries<sequence_model::tag_dna> &qs, const partassign::part_assignment &part_assign );
}
