#ifndef ECO_EC_ORG_TYPES_H
#define ECO_EC_ORG_TYPES_H

#include "hardware/AvidaGP.h"
#include "hardware/BitSorter.h"
#include "tools/BitVector.h"
#include "tools/string_utils.h"
#include "base/vector.h"
#include "base/Ptr.h"

class BitSorterOrg : public emp::BitSorter {
public:
    using emp::BitSorter::BitSorter;
    using gene_t = emp::BitSorter::bits_t;

    emp::vector<double> phenotype;
    int pos;
    double fitness;

};

class AvidaGPOrg : public emp::AvidaGP {
public:
    // AvidaGPOrg(const genome_t & in_genome) : emp::AvidaGP(in_genome) { ; }
    // AvidaGPOrg(emp::Ptr<const inst_lib_t> inst_lib) : emp::AvidaGP(Genome(inst_lib)) { ; }
    // AvidaGPOrg(const inst_lib_t & inst_lib) : emp::AvidaGP(Genome(&inst_lib)) { ; }

    using emp::AvidaGP::AvidaGP;

    AvidaGPOrg() = default;
    AvidaGPOrg(const AvidaGPOrg &) = default;
    AvidaGPOrg(AvidaGPOrg &&) = default;

    emp::vector<double> phenotype;
    int pos;
    double fitness;

    using gene_t = emp::AvidaGP::inst_t;
};

class BitVectorOrg : public emp::BitVector {
public:
    using emp::BitVector::BitVector;
    using gene_t = bool;

    emp::vector<double> phenotype;
    int pos;
    double fitness;
};

class RealValueVectorOrg : public emp::vector<double> {
public:
    RealValueVectorOrg(size_t N) : emp::vector<double>(N) {;}
    using gene_t = double;

    emp::vector<double> phenotype;
    int pos;
    double fitness;
};


using gp_t = AvidaGPOrg;
using sorting_t = BitSorterOrg;
using bit_t = BitVectorOrg;
using rv_t = RealValueVectorOrg;

namespace std {
  /// Hash function to allow BitVector to be used with maps and sets (must be in std).
  template <>
  struct hash<BitVectorOrg> {
    std::size_t operator()(const BitVectorOrg & b) const {
      return b.Hash();
    }
  };

  /// operator<< to work with ostream (must be in std to work)
  inline std::ostream & operator<<(std::ostream & out, const BitVectorOrg & bit_v) {
    bit_v.Print(out);
    return out;
  }

}

#endif