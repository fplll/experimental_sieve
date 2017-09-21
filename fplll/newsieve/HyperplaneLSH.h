//
// HyperplaneLSH.h
//

#ifndef HYPERPLANE_LSH_H
#define HYPERPLANE_LSH_H


#include "DebugAll.h"
#include "Typedefs.h"
#include <vector>
#include "SieveUtility.h"
#include "LatticeBases.h"

#include <random>





namespace GaussSieve{


    template <class SieveTraits> class Bucket_Element {

        using ListStoredPoint = typename SieveTraits::GaussList_StoredPoint;

    public:

        Bucket_Element()=default;
        Bucket_Element(const Bucket_Element & new_pointer_to_lattice_point) = delete;
        Bucket_Element(Bucket_Element && new_pointer_to_lattice_point) = default;

        Bucket_Element& operator=(Bucket_Element const &that) =delete;
        Bucket_Element& operator=(Bucket_Element && that) =default;
        
        inline bool operator==(Bucket_Element const &another) const
        {
            if(another == this->pointer_to_lattice_point) return true;
            return false;
        }


        ~Bucket_Element() {}

        Bucket_Element (ListStoredPoint const * const v)
        {
            this->pointer_to_lattice_point = v;
        }

        inline ListStoredPoint const& get_point() const {return *pointer_to_lattice_point;}
        //inline ListStoredPoint const* get_pointer() const {return pointer_to_lattice_point;}

    private:
        ListStoredPoint const* pointer_to_lattice_point;

    };


    template<class SieveTraits, class ET> struct HashTable{

        using Bucket = std::list<Bucket_Element<SieveTraits>>;


        HashTable() = default;
        HashTable(HashTable const &obj) = delete;
        HashTable(HashTable  &&obj) = delete;

        HashTable & operator=(HashTable const & obj) = delete;
        HashTable & operator=(HashTable  && obj) = delete;

        //Bucket* &operator[](int idx) { return hash_tables[idx]; };
        //Bucket const &operator[](int idx) const { return hash_tables[idx]; };

        ~HashTable() = default;

    public:
        
        void print_ith_table(unsigned short i);
        

        Bucket& candidates (int hash_val) {return hash_table[hash_val];};
        
        uint_fast16_t get_hash_fnct (int ind, int coord) {return hash_fnct[ind][coord];};

        constexpr unsigned short get_num_of_tables() const {return SieveTraits::number_of_hash_tables;};
        
        int hash_to_t (typename SieveTraits::GaussList_StoredPoint const& v, int t);

    private:
        Bucket hash_table[1 << SieveTraits::number_of_hash_functions];
        uint_fast16_t hash_fnct[SieveTraits::number_of_hash_functions][2];
        //unsigned short A[number_of_hash_tables][number_of_hash_functions][2];
        //unsigned short num_of_tables;

    };
    
    template<class SieveTraits, class ET> struct HashTableS{

        HashTableS() = default;
        HashTableS(HashTableS const &obj) = delete;
        HashTableS(HashTableS &&obj) = delete;
        
        HashTableS & operator=(HashTableS const & obj) = delete;
        HashTableS & operator=(HashTableS  && obj) = delete;
        
        ~HashTableS() = default;
        
    public:
        
        void add_to_hash_tables (typename SieveTraits::GaussList_StoredPoint const* v);
        void remove_from_hash_tables(int const &hash_value, typename SieveTraits::GaussList_StoredPoint const* v, int table_index);
        
        
        void initialize_hash_tables(typename SieveTraits::DimensionType N);
        void print_all_tables();
        
         HashTable<SieveTraits, ET>& get_ith_hash_table (int i) {return *all_hash_tables[i];};
        
    private:
    
        HashTable<SieveTraits, ET> all_hash_tables[SieveTraits::number_of_hash_tables];

        
    };



    template<class SieveTraits, class ET>
    int HashTable<SieveTraits, ET>::hash_to_t(typename SieveTraits::GaussList_StoredPoint const& v, int t)
    {
        int res = 0;
        for(int k = 0; k < SieveTraits::number_of_hash_functions; k++){
            res <<= 1;
            uint_fast16_t coord1 = this->get_hash_fnct(k,0);
            uint_fast16_t coord2 = this->get_hash_fnct(k,0);
            //uint_fast16_t coord1 = (all_hash_tables.get_ith_hash_table(t)).get_hash_fnct(k,0);
            //uint_fast16_t coord2 = (all_hash_tables.get_ith_hash_table(t)).get_hash_fnct(k,1);
            ET IP = v[coord1] - v[coord2];
            if(IP > 0)
                res++;
        }

        // Merge buckets u and 2^number_of_hash_functions - u - 1
        if(res >= (1 << (SieveTraits::number_of_hash_functions-1)))
            res = 2 *  (1 << (SieveTraits::number_of_hash_functions-1)) - res - 1;
	return res;
    }

    template<class SieveTraits, class ET>
    void HashTableS<SieveTraits, ET>::add_to_hash_tables (typename SieveTraits::GaussList_StoredPoint const* v)
    {

        // the received (via a pointer) lattice point is copied into a new lattice-point
        // a pointer to the newly created lattice-pointed is stored in hash
        for(int i=0; i<SieveTraits::number_of_hash_tables; ++i)
        {
            
            int hash_value = all_hash_tables[i].hash_to_t(*v, i);

//            bool erased = false;
            /*
            if(this->hash_tables[i][hash_value].size() == SieveTraits::max_bucket_size){
                std::cout<< "a bucket is full " << std::endl;

                //find a vector longer than v and throw it away

                for (auto it = hash_tables[i][hash_value].cbegin(); it !=hash_tables[i][hash_value].cend(); ++it)
                {
                        if ((*it).get_point().get_norm2() > (*v).get_norm2() )
                        {
                            remove_from_hash_tables(&(*it).get_point(), SieveTraits::number_of_hash_tables); //we do not have number_of_hash_tables-th hash-table
                            erased = true;
                            break;
                        }
                }

                //this->print_all_tables();
                //assert(false);
            }
            */

            // Insert v into the bucket
//            if (erased)
                //(this->all_hash_tables[i]).hash_table[hash_value].emplace_back(v);
//            print_all_tables();
        }

        //std::cout <<"one element is hashed" << std::endl;

    }

    //TODO
    template<class SieveTraits, class ET>
    void HashTableS<SieveTraits, ET>::remove_from_hash_tables(int const &hash_value, 
            typename SieveTraits::GaussList_StoredPoint const* address, int table_index)
    {
        
        for(int t=0; t<SieveTraits::number_of_hash_tables; ++t)
        {
            if (t == table_index)
                continue;
          
            /*
            for (auto it = hash_tables[t][hash_value].cbegin(); it !=hash_tables[t][hash_value].cend(); ++it)
            {
                if (&(*it).get_point() == address)
                {
                    std::cout << &(*it).get_point() <<std::endl<<std::flush;
                    this->hash_tables[t][hash_value].erase(it);
                    break;
                }
            }
             */
             this->all_hash_tables[t][hash_value].remove(address);
            
        }
        //std::cout << "removed from all hash-tables" <<std::endl;
    }


    template<class SieveTraits, class ET>
    void HashTableS<SieveTraits, ET>::initialize_hash_tables(typename SieveTraits::DimensionType N)
    {
        //TODO: DIMENSION IS NEEDED
        std::cout << "N = " << N << std::endl;

        for(int t = 0; t < SieveTraits::number_of_hash_tables; t++){
            
            // Initialize random sparse hash vectors by choosing two non-zero entries
            for(int k = 0; k < SieveTraits::number_of_hash_functions; k++){
                
                //TODO: REDO
                //all_hash_tables[t].hash_fnct[k][0] = (rand() % N);
                //this->A[t][k][0] = (rand() % N);
                //this->A[t][k][1] = (rand() % N);
                /*
                while(this->A[t][k][1] == this->A[t][k][0])
                {
                    this->A[t][k][1] = (rand() % N);
                }
                 */
            }

            /*
            for(int b = 0; b <  (1 << (number_of_hash_functions-1)); b++){
                this->hash_tables[t][b].reserve(max_bucket_size);
            }
            */
        }
    }

    template<class SieveTraits, class ET>
    void HashTableS<SieveTraits, ET>::print_all_tables()
    {
    /*
    long number = 0;
    for (int t= 0; t < SieveTraits::number_of_hash_tables;++t)
    {
      for(int k = 0; k < (1<<SieveTraits::number_of_hash_functions-1);++k )
      {
        number+=hash_tables[t][k].size();
      }
    }
    std::cout << "Total number of elements in the hash tables:" << number << std::endl;
    */
        for (int t = 0; t<SieveTraits::number_of_hash_tables; ++t)
        {
            std::cout <<"Hash-table #" << t+1 << std::endl;
            for (int k =0; k< (1 << (SieveTraits::number_of_hash_functions-1)); ++k )
            {
                //std::cout << this->hash_tables[t][k].size() << std::endl;
                if (this->all_hash_tables[t][k].size() >0 )
                {
                    std::cout<< k <<"-th bucket has " << this->all_hash_tables[t][k].size() << " elements" << std::endl;
                    for (auto it=this->all_hash_tables[t][k].cbegin(); it!=this->all_hash_tables[t][k].cend(); ++it)
                        std::cout<<(*it).get_point().get_norm2()<< " ";
                    std::cout << std::endl;
                }

            }
        }
 
    }
    
    //TODO
    template<class SieveTraits, class ET>
    void HashTable<SieveTraits, ET>::print_ith_table(unsigned short i)
    {
        /*
        for (int k =0; k< (1 << (SieveTraits::number_of_hash_functions-1)); ++k )
        {
            //std::cout << this->hash_tables[t][k].size() << std::endl;
            if (this->hash_tables[i][k].size() >0 )
            {
                std::cout<< k <<"-th bucket has " << this->hash_tables[i][k].size() << " elements" << std::endl;
                for (auto it=this->hash_tables[i][k].cbegin(); it!=this->hash_tables[i][k].cend(); ++it)
                            std::cout<<(*it).get_point().get_norm2()<< " ";
                std::cout << std::endl;
            }

        }
         */ 
    }

}


#endif
