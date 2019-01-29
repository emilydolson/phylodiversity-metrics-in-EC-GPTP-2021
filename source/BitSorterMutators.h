// This file was written by Alex Lalejini

#ifndef ALEX_BITSORTER_MUTATORS_H
#define ALEX_BITSORTER_MUTATORS_H

#include "base/vector.h"
#include "tools/Random.h"
#include "tools/random_utils.h"

#include "org_types.h"

/// Bit sorter organism mutator
struct BitSorterMutator {

  size_t MAX_NETWORK_SIZE;  ///< Maximum size network can grow
  size_t MIN_NETWORK_SIZE;  ///< Minimum size network can shrink
  size_t SORT_SEQ_SIZE;     ///< Sort input size (defines range for i,j values)

  double PER_INDEX_SUB;
  double PER_PAIR_DUP;
  double PER_PAIR_INS;
  double PER_PAIR_DEL;
  double PER_PAIR_SWAP;

  BitSorterMutator() 
    : MAX_NETWORK_SIZE(64),
      MIN_NETWORK_SIZE(1),
      SORT_SEQ_SIZE(4),
      PER_INDEX_SUB(0.001),
      PER_PAIR_DUP(0.001),
      PER_PAIR_INS(0.001),
      PER_PAIR_DEL(0.001),
      PER_PAIR_SWAP(0.001)
  { ; }

  size_t Mutate(emp::Random & rnd, sorting_t & genome) {
    size_t muts = 0;
    size_t expected_size = genome.GetSize();

    // Doing this in reverse preserves indices
    for (int i = genome.GetSize() - 1; i >= 0; --i) {
      // Deletions!
      if (rnd.P(PER_PAIR_DEL) && (expected_size > MIN_NETWORK_SIZE)) {
        ++muts;
        --expected_size;
        genome.RemoveCompare(i);
        continue;
      }

      // Do we insert?
      if (rnd.P(PER_PAIR_INS) && (expected_size < MAX_NETWORK_SIZE)) {
        ++muts;
        ++expected_size;
        // Insert randomly
        genome.InsertCompare(i, rnd.GetUInt(SORT_SEQ_SIZE), rnd.GetUInt(SORT_SEQ_SIZE));
      }
      
      auto comparator = genome.GetComparator(i);
      // Do we duplicate?
      if (rnd.P(PER_PAIR_DUP) && (expected_size < MAX_NETWORK_SIZE)) {
        ++muts;
        ++expected_size;
        // Duplicate!
        genome.InsertCompare(i+1, comparator.first, comparator.second);

      }

      // Per-index substitutions?
      if (rnd.P(PER_INDEX_SUB)) {
        genome.EditCompare(i, rnd.GetUInt(SORT_SEQ_SIZE), comparator.second);
        ++muts;
      }
      if (rnd.P(PER_INDEX_SUB)) {
        genome.EditCompare(i, comparator.first, rnd.GetUInt(SORT_SEQ_SIZE));
        ++muts;
      }
    }

    // How about swaps?
    if (PER_PAIR_SWAP > 0) {
      for (size_t i = 0; i < genome.GetSize(); ++i) {
        if (rnd.P(PER_PAIR_SWAP)) {
          // Select two random positions
          auto comparator = genome.GetComparator(i);
          genome.InsertCompare(i, comparator.second, comparator.first);
          ++muts;
        }
      }
    }

    return muts;
  }

  sorting_t GenRandomBitSorter(emp::Random & rnd) {
    sorting_t rando_sorter;
    // How big?
    const size_t sorter_size = rnd.GetUInt(MIN_NETWORK_SIZE, MAX_NETWORK_SIZE);
    // Build!
    for (size_t i = 0; i < sorter_size; ++i) {
      rando_sorter.AddCompare(rnd.GetUInt(SORT_SEQ_SIZE),rnd.GetUInt(SORT_SEQ_SIZE));
    }
    return rando_sorter;
  }
};



#endif