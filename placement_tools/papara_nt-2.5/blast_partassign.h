#ifndef __blast_partassign_h
#define __blast_partassign_h

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace papara {
template<typename pvec_t, typename seq_tag>
class references;

template<typename seq_tag>
class queries;


}


namespace partassign {

struct blast_hit {
    std::string qs_name;
    std::string ref_name;
    
    int qs_start;
    int qs_end;
    int ref_start;
    int ref_end;
    
    float bit_score;
//     std::string bit_score;
//     double evalue;
    
};

class partition {
public:
    partition() : start(-1), end(-1) {}
    int start;
    int end;
    std::string gene_name;
};


// static bool not_space( char x ) {
//     return !isspace(x); // doing the same with pre-c++11 functional is ridiculus
// }

blast_hit next_hit( std::istream &is );
partition next_partition( std::istream &is ); 


class part_assignment {
public:
    part_assignment( std::istream &blast_out, std::istream &part_file ) ;
    
//     const partassign::partition &partition( const std::string &name ) const ;
    const blast_hit &get_blast_hit( const std::string &qs_name ) const;
    const std::vector<partassign::partition> &partitions() const {
        return partitions_;
    }
    
private:
    std::vector<partassign::partition> partitions_;
    //std::map<std::string,int> assignments_;
    std::map<std::string,partassign::blast_hit> hits_;
    
    
};






template<typename pvec_t, typename seq_tag>
std::vector<std::pair<size_t,size_t> > resolve_qs_bounds( papara::references<pvec_t,seq_tag> &refs, papara::queries<seq_tag> &qs, const partassign::part_assignment &part_assign ); 


std::pair<size_t,size_t> partition_bounds( std::istream &is, const std::string &name );


}
#endif