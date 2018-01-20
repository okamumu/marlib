/*
  phase-type distribution
*/


namespace marlib {

  template <typename VectorT, typename MatrixT> struct phparams;
  template <typename VectorT, typename MatrixT> struct phparams_unif;
  template <typename VectorT, typename MatrixT> struct phase_estep_results;

  // algorithms

  enum dist_t : char {
    PDF = 'P',
    CDF = 'D',
    CCDF = 'C',
  };

  /*
    prototype for algoritms
    main codes are located in phase.cpp
  */

  template <typename VectorT, typename MatrixT>
  phparams_unif<VectorT,MatrixT>& phase_unif(phparams_unif<VectorT,MatrixT>& ph);

  template <typename VectorT, typename MatrixT>
  VectorT& phase_cdf_unif(const dist_t& dtype,
    const phparams_unif<VectorT,MatrixT>& ph,
    const VectorT& t, VectorT& x);

  template <typename VectorT, typename MatrixT>
  phase_estep_results<VectorT,MatrixT>& phase_estep_wtime_unif(
    const phparams_unif<VectorT,MatrixT>& ph,
    VectorT tdat,
    VectorT wdat,
    phase_estep_results<VectorT,MatrixT>& eres);

  template <typename VectorT, typename MatrixT, typename iVectorT>
  phase_estep_results<VectorT,MatrixT>& phase_estep_group_trunc(
    const phparams_unif<VectorT,MatrixT>& ph,
    VectorT tdat,
		iVectorT gdat,
		const typename iVectorT::ValueType& gdatlast,
		iVectorT idat,
    phase_estep_results<VectorT,MatrixT>& eres);

  template <typename VectorT, typename MatrixT, typename iVectorT>
  phase_estep_results<VectorT,MatrixT>& phase_estep_group_trunc_poi(
    const phparams_unif<VectorT,MatrixT>& ph,
    const typename VectorT::ValueType& omega,
    VectorT tdat,
		iVectorT gdat,
		const typename iVectorT::ValueType& gdatlast,
		iVectorT idat,
    phase_estep_results<VectorT,MatrixT>& eres);

  template <typename VectorT, typename MatrixT>
  phparams<VectorT,MatrixT>& phase_mstep_gen(
    const phase_estep_results<VectorT,MatrixT>& eres,
    phparams<VectorT,MatrixT>& ph);

  template <typename VectorT, typename MatrixT>
  void bidiag2cf1(phparams<VectorT,MatrixT>& ph);

  /*
    phase params
  */

  template <typename VectorT, typename MatrixT>
  struct phparams {
    phparams(const VectorT& _alpha, const VectorT& _xi, const MatrixT& _T)
    : alpha(_alpha), baralpha(_alpha.clone()), xi(_xi), T(_T) {}

    size_type dim() const {
      return alpha.size();
    }

    ////// print
    virtual std::ostream& print(std::ostream& os) const {
      os << "alpha=(" << alpha << ")" << std::endl;
      os << "xi=(" << xi << ")" << std::endl;
      os << "T=" << std::endl;
      os << T << std::endl;
      return os;
    }

    template <typename VectorTT, typename MatrixTT>
    friend std::ostream& operator<< (std::ostream& os, const phparams<VectorTT,MatrixTT>& v);

    VectorT alpha;
    VectorT baralpha;
    VectorT xi;
    MatrixT T;

    virtual phparams<VectorT,MatrixT>& mstep(
      const phase_estep_results<VectorT,MatrixT>& eres) = 0;

    virtual VectorT& cdf(const dist_t& dtype, const VectorT& t, VectorT& x) const = 0;

  };

  template <typename VectorT, typename MatrixT>
  std::ostream& operator<<(std::ostream& os, const phparams<VectorT,MatrixT>& v) {
    return v.print(os);
  }

  /*
    phase params with uniformization
  */

  template <typename VectorT, typename MatrixT>
  struct phparams_unif : public phparams<VectorT,MatrixT> {

    const typename VectorT::ValueType unif_factor = 1.01;
    const typename VectorT::ValueType pdf_eps = 1.0e-8;

    typedef phparams<VectorT,MatrixT> super;

    phparams_unif(const VectorT& _alpha, const VectorT& _xi, const MatrixT& _T)
    : super(_alpha, _xi, _T), P(_T.clone()), qv(0) { }

    virtual std::ostream& print(std::ostream& os) const {
      super::print(os);
      os << "P=" << std::endl;
      os << P << std::endl;
      os << "qv=" << qv << std::endl;
      return os;
    }

    MatrixT P;
    typename VectorT::ValueType qv;

    // methods
    phparams_unif<VectorT,MatrixT>& unif() {
      return phase_unif(*this);
    }
  };

  /*
    estep results
  */

  template <typename VectorT, typename MatrixT>
  struct phase_estep_results {
    typename VectorT::ValueType llf;
    typename VectorT::ValueType etotal;
    VectorT eb;
    VectorT ey;
    VectorT ez;
    MatrixT en;

    phase_estep_results(const phparams<VectorT,MatrixT>& ph)
    : llf(0), etotal(0), eb(ph.dim()), ey(ph.dim()), ez(ph.dim()), en(ph.T.clone()) {}

    // for test
    virtual std::ostream& print(std::ostream& os) const {
      os << "total=" << etotal << std::endl;
      os << "eb=" << eb << " " << dasum(eb) << std::endl;
      os << "ey=" << ey << " " << dasum(ey) << std::endl;
      os << "ez=" << ez << " " << dasum(ez) << std::endl;
      os << "en=" << std::endl;
      os << en << std::endl;
      os << "llf=" << llf << std::endl;
      return os;
    }
  };

  template <typename VectorT, typename MatrixT>
  std::ostream& operator<<(std::ostream& os, const phase_estep_results<VectorT,MatrixT>& v) {
    return v.print(os);
  }

  /*
    gph unif
  */

  template <typename VectorT, typename MatrixT, typename iVectorT>
  struct gph_unif : public phparams_unif<VectorT,MatrixT> {

    typedef phparams_unif<VectorT,MatrixT> super;

    gph_unif(const VectorT& _alpha, const VectorT& _xi, const MatrixT& _T)
    : super(_alpha, _xi, _T) { }

    // methods
    VectorT& cdf(const dist_t& dtype,
      const VectorT& t, VectorT& x) const {
      return phase_cdf_unif(dtype, *this, t, x);
    }

    phase_estep_results<VectorT,MatrixT>& estep(
      const VectorT& tdat,
      const VectorT& wdat,
      phase_estep_results<VectorT,MatrixT>& eres) const {
      return phase_estep_wtime_unif(*this, tdat, wdat, eres);
    }

    phase_estep_results<VectorT,MatrixT>& estep(
      const VectorT& tdat,
      const iVectorT& gdat,
      const typename iVectorT::ValueType& gdatlast,
      const iVectorT& idat,
      phase_estep_results<VectorT,MatrixT>& eres) const {
      return phase_estep_group_trunc(*this, tdat, gdat, gdatlast, idat, eres);
    }

    phase_estep_results<VectorT,MatrixT>& estep(
      const typename VectorT::ValueType& omega,
      const VectorT& tdat,
      const iVectorT& gdat,
      const typename iVectorT::ValueType& gdatlast,
      const iVectorT& idat,
      phase_estep_results<VectorT,MatrixT>& eres) const {
      return phase_estep_group_trunc_poi(*this, omega, tdat, gdat, gdatlast, idat, eres);
    }

    phparams<VectorT,MatrixT>& mstep(
      const phase_estep_results<VectorT,MatrixT>& eres) {
      return phase_mstep_gen(eres, *this);
    }
  };

  /*
    cf1 unif
  */

  template <typename VectorT, typename MatrixT, typename iVectorT>
  struct cf1_unif : public phparams_unif<VectorT,MatrixT> {

    typedef phparams_unif<VectorT,MatrixT> super;

    cf1_unif(const VectorT& _alpha, const VectorT& _xi, const MatrixT& _T)
    : super(_alpha, _xi, _T) { }

    // methods
    VectorT& cdf(const dist_t& dtype,
      const VectorT& t, VectorT& x) const {
      return phase_cdf_unif(dtype, *this, t, x);
    }

    phase_estep_results<VectorT,MatrixT>& estep(
      const VectorT& tdat,
      const VectorT& wdat,
      phase_estep_results<VectorT,MatrixT>& eres) const {
      return phase_estep_wtime_unif(*this, tdat, wdat, eres);
    }

    phase_estep_results<VectorT,MatrixT>& estep(
      const VectorT& tdat,
      const iVectorT& gdat,
      const typename iVectorT::ValueType& gdatlast,
      const iVectorT& idat,
      phase_estep_results<VectorT,MatrixT>& eres) const {
      return phase_estep_group_trunc(*this, tdat, gdat, gdatlast, idat, eres);
    }

    phase_estep_results<VectorT,MatrixT>& estep(
      const typename VectorT::ValueType& omega,
      const VectorT& tdat,
      const iVectorT& gdat,
      const typename iVectorT::ValueType& gdatlast,
      const iVectorT& idat,
      phase_estep_results<VectorT,MatrixT>& eres) const {
      return phase_estep_group_trunc_poi(*this, omega, tdat, gdat, gdatlast, idat, eres);
    }

    phparams<VectorT,MatrixT>& mstep(
      const phase_estep_results<VectorT,MatrixT>& eres) {
      phase_mstep_gen(eres, *this);
      bidiag2cf1(*this);
      return *this;
    }
  };
}
