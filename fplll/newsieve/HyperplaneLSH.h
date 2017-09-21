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
      
      Bucket_Element (ListStoredPoint const * const v)
      {
        this->pointer_to_lattice_point = v;
      }

      Bucket_Element& operator=(Bucket_Element const &that) =delete;
      Bucket_Element& operator=(Bucket_Element && that) =default;
        
      inline bool operator==(Bucket_Element const &another) const
      {
          if(another == this->pointer_to_lattice_point) return true;
          return false;
      }

      ~Bucket_Element() {}

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
        
      void print_ith_table();
      

      Bucket& candidates (int hash_val) {return hash_table[hash_val];};
      
      uint_fast16_t get_hash_fnct (int ind, int coord) {return hash_fnct[ind][coord];};

      constexpr unsigned short get_num_of_tables() const {return SieveTraits::number_of_hash_tables;};
      
      int hash_to_t (typename SieveTraits::GaussList_StoredPoint const& v);
      
      void add_to_hash_table (typename SieveTraits::GaussList_StoredPoint const* v);
      
      void initialize_hash_table(typename SieveTraits::DimensionType N);

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
        
      
      void remove_from_hash_tables(int const &hash_value, typename SieveTraits::GaussList_StoredPoint const* v, int table_index);
      
      void add_to_all_hash_tables (typename SieveTraits::GaussList_StoredPoint const& v);
        
      
      void initialize_hash_tables(typename SieveTraits::DimensionType N);
      void print_all_tables();
        
         HashTable<SieveTraits, ET>& get_ith_hash_table (int i) {return *all_hash_tables[i];};
        
    private:
    
        HashTable<SieveTraits, ET> all_hash_tables[SieveTraits::number_of_hash_tables];

        
    };


    template<class SieveTraits, class ET>
    int HashTable<SieveTraits, ET>::hash_to_t(typename SieveTraits::GaussList_StoredPoint const& v)
    {
        int res = 0;
        for(int k = 0; k < SieveTraits::number_of_hash_functions; k++){
            res <<= 1;
            uint_fast16_t coord1 = this->get_hash_fnct(k,0);
            uint_fast16_t coord2 = this->get_hash_fnct(k,1);
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
    void HashTable<SieveTraits, ET>::add_to_hash_table(typename SieveTraits::GaussList_StoredPoint const* v)
    {
      int hash_value = hash_to_t(*v);
      this->hash_table[hash_value].emplace_back(v);
      
    }

    template<class SieveTraits, class ET>
    void HashTableS<SieveTraits, ET>::add_to_all_hash_tables (typename SieveTraits::GaussList_StoredPoint const& v)
    {

      // the received (via a pointer) lattice point is copied into a new lattice-point
      // a pointer to the newly created lattice-pointed is stored in hash
      
      //typename SieveTraits::GaussList_StoredPoint v_copy = v.make_copy();
      typename SieveTraits::GaussList_StoredPoint* v_copy  = new typename SieveTraits::GaussList_StoredPoint(v.make_copy());
      
      for(int i=0; i<SieveTraits::number_of_hash_tables; ++i)
      {
            
          this->all_hash_tables[i].add_to_hash_table(v_copy);

      }

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
  void HashTable<SieveTraits, ET>::initialize_hash_table(typename SieveTraits::DimensionType N)
  {
    for(int k = 0; k < SieveTraits::number_of_hash_functions; k++){
      this->hash_fnct[k][0] = (rand() % N);
      this->hash_fnct[k][1] = (rand() % N);
      
      while (this->hash_fnct[k][0] == this->hash_fnct[k][1])
      {
        this->hash_fnct[k][1] = (rand() % N);
      }
      
    }
  }


  template<class SieveTraits, class ET>
  void HashTableS<SieveTraits, ET>::initialize_hash_tables(typename SieveTraits::DimensionType N)
  {
      
    std::cout << "N = " << N << std::endl;
    for(int t = 0; t < SieveTraits::number_of_hash_tables; t++)
    {
      all_hash_tables[t].initialize_hash_table(N);
    }
  }

  template<class SieveTraits, class ET>
  void HashTableS<SieveTraits, ET>::print_all_tables()
  {
    //TODO: REDO
    /*
    long total_hashes = 0;
    for (int t= 0; t < SieveTraits::number_of_hash_tables;++t)
    {
      for(int k = 0; k < (1<<SieveTraits::number_of_hash_functions-1);++k )
      {
        total_hashes+=all_hash_tables[t][k].size();
      }
    }
  std::cout << "Total number of elements in the hash tables:" << total_hashes << std::endl;
  */
      for (int t = 0; t<SieveTraits::number_of_hash_tables; ++t)
      {
        std::cout <<"Hash-table #" << t+1 << std::endl;
        this->all_hash_tables[t].print_ith_table();
      }

  }
    
    //TODO
  template<class SieveTraits, class ET>
  void HashTable<SieveTraits, ET>::print_ith_table()
  {
    
    for (int k =0; k< (1 << (SieveTraits::number_of_hash_functions-1)); ++k )
    {
      if (this->hash_table[k].size() >0 )
      {
        std::cout<< "H[" <<k <<"]:  ";
        for (auto it=this->hash_table[k].cbegin(); it!=this->hash_table[k].cend(); ++it)
                          std::cout<<(*it).get_point().get_norm2()<< " ";
        std::cout << std::endl;
        }
      }
    }
  }


#endif
