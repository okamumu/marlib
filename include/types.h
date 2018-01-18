#pragma once

#include <iostream>

namespace marlib {

  typedef int size_type;

  enum trans_t : char {
    Trans = 'T',
    NoTrans = 'N',
  };

  inline
  trans_t trans_c(const trans_t& trans) {
    switch (trans) {
      case Trans:
      return NoTrans;
      case NoTrans:
      return Trans;
    }
  }

}
