/**
  DefaultTermConds.h

  This file provides default and example termination conditions.
  These inherit from TerminationCondition, declared in TerminationConditions.h
  See that file for an explanation of the interface.
*/

// clang-format: currently adheres to clang-format

#ifndef DEFAULT_TERM_CONDS_H
#define DEFAULT_TERM_CONDS_H

#include "DefaultIncludes.h"
#include "TerminationConditions.h"  // consider merging with that file

namespace GaussSieve
{

// default Termination Conditions here:

/**
  This Termination condition never triggers.
*/

template <class SieveTraits, bool MT>
class NeverTerminationCondition final : public TerminationCondition<SieveTraits, MT>
{
public:
  virtual int check(typename SieveTraits::TermCond_QueryType const &lattice_point) override
  {
    return 0;
  }
  virtual ~NeverTerminationCondition(){};
  virtual TerminationConditionType termination_condition_type() const override
  {
    return TerminationConditionType::never_terminate;
  }
};

/**
  LengthTerminationCondition stores a norm^2 as internal data.
  It checks whether the vectors are shorter than this.
  The length bound is set during construction and cannot be modified currently.
*/

template <class SieveTraits, bool MT>
class LengthTerminationCondition : public TerminationCondition<SieveTraits, MT>
{

private:  // shorthands
  using LengthType         = typename SieveTraits::LengthType;
  using TermCond_QueryType = typename SieveTraits::TermCond_QueryType;
  using TerminationCondition<SieveTraits, MT>::sieveptr;

public:
  explicit LengthTerminationCondition(LengthType const &init_target_norm2)
      : target_norm2(init_target_norm2)
  {
  }
  virtual inline int check_vec(TermCond_QueryType const &lattice_point) override;
  virtual inline int check() override;

  // run-time type information:
  virtual TerminationConditionType termination_condition_type() const override
  {
    return TerminationConditionType::length_condition;
  }
  virtual TermCondDependencyType dependency_type() const override
  {
    // only needs to be called if a new shortest vector is found.
    return TermCondDependencyType::new_shortest_vector;
  }

  virtual ~LengthTerminationCondition() {}

  // exception in reading dumps -- unused
  class bad_dumpread_LengthTermCond : public std::runtime_error
  {
  public:
    bad_dumpread_LengthTermCond() : runtime_error("Dump read failed for LengthTerminationCondition")
    {
    }
  };  // exception indicating that read from dump failed.

  // dumping to stream
  virtual std::ostream &dump_to_stream(std::ostream &os) const override
  {
    os << "Target Norm^2 =" << target_norm2 << std::endl;
    return os;
  }

  // reading back in:
  virtual std::istream &read_from_stream(std::istream &is) override
  {
    if (!string_consume(is, "Target Norm^2 ="))
    {
      throw bad_dumpread_LengthTermCond();
    }
    is >> target_norm2;
    return is;
  }

private:
  LengthType target_norm2;
};

/**
  MinkowskiTerminationCondition checks by length, as above.
  The length is computed from the basis by a reasonable version of Minkowski's bound.

  Currently only works for certain data types of bases.
  Consider refactoring template dependencies...
*/

template <class SieveTraits, bool MT>
class MinkowskiTerminationCondition final : public TerminationCondition<SieveTraits, MT>
{
private:  // shorthands
  using LengthType         = typename SieveTraits::LengthType;
  using TermCond_QueryType = typename SieveTraits::TermCond_QueryType;
  using TerminationCondition<SieveTraits, MT>::sieveptr;

public:
  // target_length may be unitialised. We are guaranteed that init() is run before use.
  MinkowskiTerminationCondition() : target_norm2() {}

  virtual inline int check() override;
  virtual inline int check_vec(TermCond_QueryType const &lattice_point) override;
  virtual void custom_init() override;

  // serves as run-time type information:
  virtual TerminationConditionType termination_condition_type() const override
  {
    return TerminationConditionType::minkowski_condition;
  }
  virtual TermCondDependencyType dependency_type() const override
  {
    // only needs to be called if a new shortest vector is found.
    return TermCondDependencyType::new_shortest_vector;
  }
  virtual ~MinkowskiTerminationCondition() {}

private:
  LengthType target_norm2;
};

}  // namespace GaussSieve

#endif
