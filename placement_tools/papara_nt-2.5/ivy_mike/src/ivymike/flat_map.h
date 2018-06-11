#ifndef ivy_mike__flat_map__h
#define ivy_mike__flat_map__h
namespace ivy_mike {
template<typename K, typename V>
class flat_map {
public:
    flat_map() : sorted_(true) {}
    
    void sort() {
        if( !sorted_ ) {
            std::sort( pairs_.begin(), pairs_.end() );
            sorted_ = true;
        }
    }
    
    void put_fast( const K &key, const V &value ) {
        pairs_.emplace_back( key, value );
        
        sorted_ = false;
    }
    
    void put_fast( const K &key, V &&value ) {
        pairs_.emplace_back( key, value );
        
        sorted_ = false;
    }
    template<typename... Args>
    void emplace_fast( const K &key, Args&&... args) {
        pairs_.emplace_back( key, args... );
        
        sorted_ = false;
    
    }
    
    
    void put( const K &key, const V &value ) {
        if( !sorted_ ) {
            throw std::runtime_error( "flat_map::put on unsorted map" );
        }
        
        ipair p{key, value};
        auto lb = std::lower_bound( pairs_.begin(), pairs_.end(), p );
        
        pairs_.insert(lb, std::move(p));
    }
    
    const V * get( const K &key ) const {
        if( !sorted_ ) {
            throw std::runtime_error( "flat_map::get on unsorted map" );
        }
        auto lb = std::lower_bound( pairs_.begin(), pairs_.end(), ipair(key, V()) );
        
        if( lb == pairs_.end() || lb->key_ != key ) {
            return nullptr;
        } else {
            return &lb->value_;
        }
    }
    
    void reserve( size_t s ) {
        pairs_.reserve(s);
    }
        
    
private:
    struct ipair {
        K key_;
        V value_;
        
        ipair( const K &key, const V &value ) : key_(key), value_(value) {}
        
        
        inline bool operator<( const ipair &other ) const {
            return key_ < other.key_;
        }
        
        
    };
    
    
    std::vector<ipair> pairs_;
    bool sorted_;
};   
    
    
}

#endif