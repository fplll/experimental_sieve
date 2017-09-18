//
// HyperplaneLSH.h
//

#ifndef HYPERPLANE_LSH_H
#define HYPERPLANE_LSH_H


#include "Typedefs.h"
#include <vector>
#include "SieveUtility.h"
#include "DebugAll.h"
#include "LatticeBases.h"


#define T 12
#define K 10
#define MaxBucketSize 13


namespace GaussSieve{
    
    
    template <class SieveTraits> class Bucket_Element {
        
        using ListStoredPoint = typename SieveTraits::GaussList_StoredPoint;
        
    public:
        
        Bucket_Element()=delete;
        Bucket_Element(const Bucket_Element &pointer_to_lattice_point) = delete; 
        Bucket_Element(Bucket_Element &&pointer_to_lattice_point) = default;
        
        Bucket_Element& operator=(Bucket_Element const &that) =delete;
        Bucket_Element& operator=(Bucket_Element && that) =default;


        ~Bucket_Element() {}
        
        Bucket_Element (ListStoredPoint const * v)
        {
            this->pointer_to_lattice_point = v;
        }
        
        inline ListStoredPoint const& get_point() const {return *pointer_to_lattice_point;}
        
        
    private:
        ListStoredPoint const* pointer_to_lattice_point;
    
    };
    
    
    template<class SieveTraits, class ET> struct HashTablesClass{
        
        using Bucket = std::vector<Bucket_Element<SieveTraits>>;
        
        
        HashTablesClass() = default;
        HashTablesClass(HashTablesClass const &obj) = delete;
        HashTablesClass(HashTablesClass  &&obj) = delete;
        
        HashTablesClass & operator=(HashTablesClass const & obj) = delete;
        HashTablesClass & operator=(HashTablesClass  && obj) = delete;
        
        //Bucket* &operator[](int idx) { return HashTables[idx]; };
        //Bucket const &operator[](int idx) const { return HashTables[idx]; };
        
        ~HashTablesClass() = default;
        
    public:
        void initialize_hash_tables(typename SieveTraits::DimensionType N);
        
        void add_to_hash_tables (typename SieveTraits::GaussList_StoredPoint const* v);
        void remove_from_hash_tables(typename SieveTraits::GaussList_StoredPoint* v);
        void delete_from_hash_tables (typename SieveTraits::GaussList_StoredPoint* v);
        int hash (typename SieveTraits::GaussList_StoredPoint const& v, int t);
        
        Bucket& candidates (int ind1, int ind2) {return HashTables[ind1][ind2];};
        
        unsigned short get_num_of_tables() const {return T;};
        
        void print();
        
    private:
        Bucket HashTables[T][1 << (K-1)];
        uint_fast16_t A[T][K][2];
        //unsigned short A[T][K][2];
        //unsigned short num_of_tables;
        
        
    };
    
    
    
    template<class SieveTraits, class ET>
    int HashTablesClass<SieveTraits, ET>::hash(typename SieveTraits::GaussList_StoredPoint const& v, int t)
    {
        int res = 0;
        for(int k = 0; k < K; k++){
            res <<= 1;
            uint_fast16_t coord1 = this->A[t][k][0];
            uint_fast16_t coord2 = this->A[t][k][1];
            ET IP = v[coord1] - v[coord2];
            if(IP > 0)
                res++;
        }
	
        // Merge buckets u and 2^K - u - 1
        if(res >= (1 << (K-1)))
            res = 2 *  (1 << (K-1)) - res - 1;
	return res;
    }
    
    // Add a vector v to a hash table bucket b
    
    /*
    template<class SieveTraits>
    void add_to_bucket(typename SieveTraits::Bucket* b, typename SieveTraits::GaussList_ReturnType const* v)
    {
        // TODO
        if(b.size() == MaxBucketSize){
                std::cout<< "a bucket is full " << std::endl;
                assert(false);
        }
	
        // Insert v into the bucket
         b.emplace_back(v);
    }
    */
    template<class SieveTraits, class ET>
    void HashTablesClass<SieveTraits, ET>::add_to_hash_tables (typename SieveTraits::GaussList_StoredPoint const* v)
    {
        for(int i=0; i<T; ++i)
        {
            int hash_value = hash(*v, i);
            std::cout << "h(x) = " << hash_value << std::endl;
            
            // TODO
            if(this->HashTables[i][hash_value].size() == MaxBucketSize){
                std::cout<< "a bucket is full " << std::endl;
                assert(false);
            }
	
            // Insert v into the bucket
            
            this->HashTables[i][hash_value].emplace_back(v);
        }
        
        std::cout <<"one element is hashed" << std::endl;
        
    }
    
    
    template<class SieveTraits, class ET>
    void HashTablesClass<SieveTraits, ET>::remove_from_hash_tables(typename SieveTraits::GaussList_StoredPoint* v)
    {
        /*
        // Find w's position in the hash bucket
        int vPos = 0;
        while(b[vPos] != v && vPos < b.size()){
            vPos++;
        }
        if(vPos >= b.size()){
            perror("Vector not found in bucket...\n");
            exit(-1);
        }		
        //TODO:CHECK
        b.erase(vPos);
         */
    }
    
    
    template<class SieveTraits, class ET>
    void HashTablesClass<SieveTraits, ET>::initialize_hash_tables(typename SieveTraits::DimensionType N)
    {
        //TODO: DIMENSION IS NEEDED
        std::cout << "N = " << N << std::endl;
        
        // Initialize hash tables as empty
        
        for(int t = 0; t < T; t++){
            // Initialize random sparse hash vectors by choosing two non-zero entries
            for(int k = 0; k < K; k++){
                this->A[t][k][0] = (rand() % N);
                this->A[t][k][1] = (rand() % N);
                while(this->A[t][k][1] == this->A[t][k][0])
                {
                    this->A[t][k][1] = (rand() % N);
                }
            }
            
            /*
            for(int b = 0; b <  (1 << (K-1)); b++){
                this->HashTables[t][b].reserve(MaxBucketSize);
            }
             */
        }
    }
    
    template<class SieveTraits, class ET>
    void HashTablesClass<SieveTraits, ET>::print()
    {
        for (int t = 0; t<T; ++t)
        {
            std::cout <<"Hash-table #" << t+1 << std::endl;
            for (int k =0; k< (1 << (K-1)); ++k )
            {
                //std::cout << this->HashTables[t][k].size() << std::endl;
                if (this->HashTables[t][k].size() > 0 )
                {
                    std::cout<< k <<"-th bucket has " << this->HashTables[t][k].size() << " elements" << std::endl; 
                    //for (auto it=this->HashTables[t][k].cbegin(); it!=this->HashTables[t][k].cend(); ++it)
                    //    std::cout<<(*it).get_point().get_norm2()<< " ";
                }
                
            }
        }
    }
    
}


#endif