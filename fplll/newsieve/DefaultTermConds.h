/**
  DefaultTermConds.h

  This file provides default and example termination conditions.
  These inherit from TerminationCondition, declared in TermCondNew.h
  See that file for an explanation of the interface.
*/

// clang-format: currently adheres to clang-format

#ifndef DEFAULT_TERM_CONDS_H
#define DEFAULT_TERM_CONDS_H

#include "DefaultIncludes.h"
#include "TermCondNew.h"  // consider merging with that file

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
  using LengthType = typename SieveTraits::LengthType;
  virtual int check(Sieve<SieveTraits, MT> *const sieve) override { return 0; };
  virtual int check_vec(Sieve<SieveTraits, MT> *const sieve, LengthType const &length2) override
  {
    return 0;
  };
  virtual bool is_simple() const override { return true; };
  virtual ~NeverTerminationCondition(){};
  virtual TerminationConditionType termination_condition_type() const override
  {
    return TerminationConditionType::never_terminate;
  };  // run-time type information.
};

/**
  LengthTerminationCondition stores a length as internal data.
  It checks whether the vectors are shorter than this length.
  The length bound is set during construction and cannot be modified currently.
*/

template <class SieveTraits, bool MT>
class LengthTerminationCondition : public TerminationCondition<SieveTraits, MT>
{
public:
  using LengthType = typename SieveTraits::LengthType;
  explicit LengthTerminationCondition(LengthType const &init_target_length)
      : target_length(init_target_length){};
  virtual inline int check(Sieve<SieveTraits, MT> *const sieve) override;
  virtual int check_vec(Sieve<SieveTraits, MT> *const sieve, LengthType const &length2) override
  {
    return (length2 <= target_length) ? 1 : 0;
  };
  virtual bool is_simple() const override { return true; };
  virtual TerminationConditionType termination_condition_type() const override
  {
    return TerminationConditionType::length_condition;
  };  // run-time type information.
  virtual ~LengthTerminationCondition(){};
  class bad_dumpread_LengthTermCond : public std::runtime_error
  {
  public:
    bad_dumpread_LengthTermCond() : runtime_error("Dump read failed for LengthTerminationCondition")
    {
    }
  };  // exception indicating that read from dump failed.
  virtual std::ostream &dump_to_stream(std::ostream &os) override
  {
    os << "Target Length=" << target_length << std::endl;
    return os;
  }
  virtual std::istream &read_from_stream(std::istream &is) override
  {
    if (string_consume(is, "Target Length="))
      throw bad_dumpread_LengthTermCond();
    is >> target_length;
    return is;
  };

private:
  LengthType target_length;  // TODO: Make const (requires change to read_from_stream)
};

/**
  MinkowskiTerminationCondition checks by length, as above.
  The length is computed from the basis by a reasonable version of Minkowski's bound.

  Currently only works for certain data types of bases.
  Consider refactoring template dependencies...
*/

template <class SieveTraits, bool MT>
class MinkowskiTerminationCondition final
    : public TerminationCondition<SieveTraits, MT>  // Length Termination Condition
{
public:
  using LengthType = typename SieveTraits::LengthType;
  // target_length may be unitialised. We are guaranteed that init() is run before use.
  MinkowskiTerminationCondition() : target_length(){};

  // returns (sieve -> get_best_length2() <= target_length ) ? 1 : 0;
  virtual inline int check(Sieve<SieveTraits, MT> *const sieve) override;

  virtual int check_vec(Sieve<SieveTraits, MT> *const sieve, LengthType const &length2) override
  {
    return (length2 <= target_length) ? 1 : 0;
  }
  virtual bool is_simple() const override { return true; };

  // serves as run-time type information:
  virtual TerminationConditionType termination_condition_type() const override
  {
    return TerminationConditionType::minkowski_condition;
  }
  virtual ~MinkowskiTerminationCondition(){};
  inline virtual void init(Sieve<SieveTraits, MT> *const sieve) override;

private:
  LengthType target_length;
};

}  // namespace GaussSieve

#endif
