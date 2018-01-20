/*
  phase-type distribution
*/

#include "marlib.hpp"

namespace marlib {

  template <typename VectorT, typename MatrixT>
  phparams_unif<VectorT,MatrixT>& phase_unif(phparams_unif<VectorT,MatrixT>& ph) {
    ph.qv = damax(ph.T.diag()) * ph.unif_factor;
    ph.P = ph.T;
    ph.P /= ph.qv;
    ph.P.diag() += 1;
    return ph;
  }

  template <typename VectorT, typename MatrixT>
  VectorT& phase_cdf_unif(const dist_t& dtype,
    const phparams_unif<VectorT,MatrixT>& ph,
    const VectorT& t, VectorT& x) {

    using ValueT = typename VectorT::ValueType;
    using RangeT = typename VectorT::RangeType;

    VectorT y = ph.alpha.clone();
    poisson<ValueT,RangeT> pois(ph.qv * damax(t), ph.pdf_eps);
    for (RangeT i=t.begin(), j=x.begin(); i<=t.end(); i++, j++) {
      pois.set(ph.qv * t(i), ph.pdf_eps);
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, y, y, 0);
      switch (dtype) {
        case PDF:
          x(j) = ddot(y, ph.xi);
          break;
        case CDF:
          x(j) = 1 - dasum(y);
          break;
        case CCDF:
          x(j) = dasum(y);
          break;
      }
    }
    return x;
  }

  template <typename VectorT, typename MatrixT>
  phase_estep_results<VectorT,MatrixT>& phase_estep_wtime_unif(
    const phparams_unif<VectorT,MatrixT>& ph,
    VectorT tdat,
    VectorT wdat,
    phase_estep_results<VectorT,MatrixT>& eres) {

    using ValueT = typename VectorT::ValueType;
    using RangeT = typename VectorT::RangeType;

    size_type dim = ph.alpha.size();
    size_type m = tdat.size();
    VectorT blf(m);
    tdat.set_range(range<RangeT>(1,m));
    wdat.set_range(range<RangeT>(1,m));
    blf.set_range(range<RangeT>(1,m));

    array<VectorT*> vf(m+1);
    array<VectorT*> vb(m+1);
    array<VectorT*> vc(m+1);
    ValueT* pool = new ValueT [dim * 3 * (m + 1)];
    ValueT* p = pool;
    for (RangeT k=0; k<=m; k++) {
      vf.ptr(k) = new VectorT(dim, p);
      p += dim;
      vb.ptr(k) = new VectorT(dim, p);
      p += dim;
      vc.ptr(k) = new VectorT(dim, p);
      p += dim;
    }

    eres.llf = 0;
    ValueT tllf = 0;

    eres.etotal = 0;
    eres.eb = 0;
    eres.ey = 0;
    eres.ez = 0;
    eres.en = 0;

    poisson<ValueT,RangeT> pois(0, poisson<ValueT,RangeT>::rightbound(ph.qv * damax(tdat), ph.pdf_eps) + 1);

		// forward & backward
		vf[0] = ph.alpha;
		vb[0] = ph.xi;
		for (RangeT k=1; k<=m; k++) {
      pois.set(ph.qv * tdat(k), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k), ph.pdf_eps) + 1);

      // forward
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, vf[k-1], vf[k], 0);
      ValueT scalef = ddot(vf[k], ph.xi);
      vf[k] /= scalef;
      // dscal(1/scalef, vf[k]);
      daxpy(wdat(k), vf[k], eres.ey);

      // backward
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(NoTrans, ph.P, ph.qv, pois, vb[k-1], vb[k], 0);
      ValueT scaleb = ddot(ph.alpha, vb[k]);
      vb[k] /= scaleb;
      // dscal(1/scaleb, vb[k]);
			daxpy(wdat(k), vb[k], eres.eb);

      blf(k) = scalef;
      tllf += log(blf(k));
			eres.llf += wdat(k) * tllf;
			eres.etotal += wdat(k);
		}

    vc[m] = 0;
		daxpy(wdat(m)/blf(m), ph.alpha, vc[m]);
		for (RangeT k=m-1; k>=1; k--) {
      pois.set(ph.qv * tdat(k+1), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k+1), ph.pdf_eps) + 1);
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, vc[k+1], vc[k], 0);
      vc[k] /= blf(k);
			// dscal(1/blf(k), vc[k]);
			daxpy(wdat(k)/blf(k), ph.alpha, vc[k]);
		}
		for (RangeT k=1; k<=m; k++) {
      pois.set(ph.qv * tdat(k), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k), ph.pdf_eps) + 1);
      mexpconv_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, NoTrans, ph.P, ph.qv, pois, vc[k], vb[k-1], vb[k-1], eres.en);
		}

		eres.eb *= ph.alpha;
		eres.ez = eres.en.diag();
		eres.en *= ph.T;
    eres.en.diag() = 0;
		eres.ey *= ph.xi;

    for (RangeT k=0; k<=m; k++) {
      delete vf.ptr(k);
      delete vb.ptr(k);
      delete vc.ptr(k);
    }
    delete [] pool;

    return eres;
  }

  template <typename VectorT, typename MatrixT, typename iVectorT>
  phase_estep_results<VectorT,MatrixT>& phase_estep_group_trunc(
    const phparams_unif<VectorT,MatrixT>& ph,
    VectorT tdat,
		iVectorT gdat,
		const typename iVectorT::ValueType& gdatlast,
		iVectorT idat,
    phase_estep_results<VectorT,MatrixT>& eres) {

    using ValueT = typename VectorT::ValueType;
    using RangeT = typename VectorT::RangeType;
    using IntegerT = typename iVectorT::ValueType;

    size_type n = ph.dim();
    size_type m = tdat.size();
    tdat.set_range(range<RangeT>(1,m));
    gdat.set_range(range<RangeT>(1,m));
    idat.set_range(range<RangeT>(1,m));

    VectorT vf(n);
    VectorT tildevf(n);
    VectorT tildevb(n);

    VectorT wg(m+1);
    VectorT wp(m+1);
    wg.set_range(range<RangeT>(1,m+1));
    wp.set_range(range<RangeT>(1,m+1));

    array<VectorT*> barvf(m+1);
    array<VectorT*> barvb(m+1);
    array<VectorT*> vb(m+1);
    array<VectorT*> vc(m+1);
    ValueT* pool = new ValueT [n * 4 * (m + 1)];
    ValueT* p = pool;
    for (size_type k=0; k<=m; k++) {
      barvf.ptr(k) = new VectorT(n, p); p += n;
      barvb.ptr(k) = new VectorT(n, p); p += n;
      vb.ptr(k) = new VectorT(n, p); p += n;
      vc.ptr(k) = new VectorT(n, p); p += n;
    }

    eres.llf = 0;
    eres.eb = 0;
    eres.ey = 0;
    eres.ez = 0;
    eres.en = 0;

    poisson<ValueT,RangeT> pois(0, poisson<ValueT,RangeT>::rightbound(ph.qv * damax(tdat), ph.pdf_eps) + 1);

    barvf[0] = ph.baralpha;
		barvb[0] = 1;
		vb[0] = ph.xi;
		ValueT nn = 0;
		ValueT uu = 0;

    for (size_type k=1; k<=m; k++) {
      pois.set(ph.qv * tdat(k), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k), ph.pdf_eps) + 1);

      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, barvf[k-1], barvf[k], 0);
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(NoTrans, ph.P, ph.qv, pois, barvb[k-1], barvb[k], 0);
      dgemv(NoTrans, ValueT(-1), ph.T, barvb[k], ValueT(0), vb[k]);

      tildevf = barvf[k-1];
      tildevf -= barvf[k];
      tildevb = barvb[k-1];
      tildevb -= barvb[k];

      if (gdat(k) >= 0 && tdat(k) != 0) {
        ValueT tmp = ddot(ph.alpha, tildevb);
        eres.llf += gdat(k) * log(tmp) - lfact<ValueT,IntegerT>(gdat(k));
        nn += gdat(k);
        uu += tmp;
        wg(k) = gdat(k) / tmp;
        daxpy(wg(k), tildevb, eres.eb);
        daxpy(wg(k), tildevf, eres.ey);
      }
      if (idat(k) == 1) {
				dgemv(Trans, ValueT(-1), ph.T, barvf[k], ValueT(0), vf);
				ValueT tmp = ddot(ph.alpha, vb[k]);
				eres.llf += log(tmp);
				nn += 1;
        wp(k) = 1 / tmp;
				daxpy(wp(k), vb[k], eres.eb);
				daxpy(wp(k), vf, eres.ey);
			}
		}

    // for the interval [t_m, infinity)
		if (gdatlast >= 0) {
			ValueT tmp = ddot(ph.alpha, barvb[m]);
			eres.llf += gdatlast * log(tmp) - lfact<ValueT,IntegerT>(gdatlast);
			nn += gdatlast;
			uu += tmp;
      wg(m+1) = gdatlast / tmp;
			daxpy(wg(m+1), barvb[m], eres.eb);
			daxpy(wg(m+1), barvf[m], eres.ey);
		}

    // compupte weights for unobserved periods
		for (size_type k=1; k<=m; k++) {
			if (gdat(k) == -1) {
        tildevf = barvf[k-1];
        tildevf -= barvf[k];
        tildevb = barvb[k-1];
        tildevb -= barvb[k];
        wg(k) = nn / uu;
				daxpy(wg(k), tildevb, eres.eb);
				daxpy(wg(k), tildevf, eres.ey);
			}
		}
		if (gdatlast == -1) {
      wg(m+1) = nn / uu;
			daxpy(wg(m+1), barvb[m], eres.eb);
			daxpy(wg(m+1), barvf[m], eres.ey);
		}
		eres.llf += lgamma<ValueT,IntegerT>(nn + 1.0) - nn * log(uu);

    // compute vectors for convolution
		vc[m] = 0;
		daxpy(wg(m+1)-wg(m), ph.baralpha, vc[m]);
		if (idat(m) == 1) {
			daxpy(wp(m), ph.alpha, vc[m]);
		}
		for (size_type k=m-1; k>=1; k--) {
      pois.set(ph.qv * tdat(k+1), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k+1), ph.pdf_eps) + 1);
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, vc[k+1], vc[k], 0);
			daxpy(wg(k+1)-wg(k), ph.baralpha, vc[k]);
			if (idat(k) == 1) {
				daxpy(wp(k), ph.alpha, vc[k]);
			}
		}
		for (size_type k=1; k<=m; k++) {
      pois.set(ph.qv * tdat(k), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k), ph.pdf_eps) + 1);
			dger(NoTrans, wg(k+1)-wg(k), ph.baralpha, barvb[k], eres.en);
      mexpconv_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, NoTrans, ph.P, ph.qv, pois, vc[k], vb[k-1], vb[k-1], eres.en);
		}
		dger(NoTrans, wg(1), ph.baralpha, barvb[0], eres.en);

    eres.etotal = nn / uu;
		eres.eb *= ph.alpha;
		eres.ez = eres.en.diag();
		eres.en *= ph.T;
    eres.en.diag() = 0;
		eres.ey *= ph.xi;

    for (size_type k=0; k<=m; k++) {
      delete barvf.ptr(k);
      delete barvb.ptr(k);
      delete vb.ptr(k);
      delete vc.ptr(k);
    }
    delete [] pool;

    return eres;
	}

  template <typename VectorT, typename MatrixT, typename iVectorT>
  phase_estep_results<VectorT,MatrixT>& phase_estep_group_trunc_poi(
    const phparams_unif<VectorT,MatrixT>& ph,
    const typename VectorT::ValueType& omega,
    VectorT tdat,
		iVectorT gdat,
		const typename iVectorT::ValueType& gdatlast,
		iVectorT idat,
    phase_estep_results<VectorT,MatrixT>& eres) {

    using ValueT = typename VectorT::ValueType;
    using RangeT = typename VectorT::RangeType;
    using IntegerT = typename iVectorT::ValueType;

    size_type n = ph.dim();
    size_type m = tdat.size();
    tdat.set_range(range<RangeT>(1,m));
    gdat.set_range(range<RangeT>(1,m));
    idat.set_range(range<RangeT>(1,m));

    VectorT vf(n);
    VectorT tildevf(n);
    VectorT tildevb(n);

    VectorT wg(m+1);
    VectorT wp(m+1);
    wg.set_range(range<RangeT>(1,m+1));
    wp.set_range(range<RangeT>(1,m+1));

    array<VectorT*> barvf(m+1);
    array<VectorT*> barvb(m+1);
    array<VectorT*> vb(m+1);
    array<VectorT*> vc(m+1);
    ValueT* pool = new ValueT [n * 4 * (m + 1)];
    ValueT* p = pool;
    for (size_type k=0; k<=m; k++) {
      barvf.ptr(k) = new VectorT(n, p); p += n;
      barvb.ptr(k) = new VectorT(n, p); p += n;
      vb.ptr(k) = new VectorT(n, p); p += n;
      vc.ptr(k) = new VectorT(n, p); p += n;
    }

    eres.llf = 0;
    eres.eb = 0;
    eres.ey = 0;
    eres.ez = 0;
    eres.en = 0;

    poisson<ValueT,RangeT> pois(0, poisson<ValueT,RangeT>::rightbound(ph.qv * damax(tdat), ph.pdf_eps) + 1);

    barvf[0] = ph.baralpha;
		barvb[0] = 1;
		vb[0] = ph.xi;
		ValueT nn = 0;
		ValueT uu = 0;

    for (size_type k=1; k<=m; k++) {
      pois.set(ph.qv * tdat(k), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k), ph.pdf_eps) + 1);

      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, barvf[k-1], barvf[k], 0);
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(NoTrans, ph.P, ph.qv, pois, barvb[k-1], barvb[k], 0);
      dgemv(NoTrans, ValueT(-1), ph.T, barvb[k], ValueT(0), vb[k]);

      tildevf = barvf[k-1];
      tildevf -= barvf[k];
      tildevb = barvb[k-1];
      tildevb -= barvb[k];

      if (gdat(k) >= 0 && tdat(k) != 0) {
        ValueT tmp = ddot(ph.alpha, tildevb);
        eres.llf += gdat(k) * log(tmp) - lfact<ValueT,IntegerT>(gdat(k));
        nn += gdat(k);
        uu += tmp;
        wg(k) = gdat(k) / tmp;
        daxpy(wg(k), tildevb, eres.eb);
        daxpy(wg(k), tildevf, eres.ey);
      }
      if (idat(k) == 1) {
				dgemv(Trans, ValueT(-1), ph.T, barvf[k], ValueT(0), vf);
				ValueT tmp = ddot(ph.alpha, vb[k]);
				eres.llf += log(tmp);
				nn += 1;
        wp(k) = 1 / tmp;
				daxpy(wp(k), vb[k], eres.eb);
				daxpy(wp(k), vf, eres.ey);
			}
		}

    // for the interval [t_m, infinity)
		if (gdatlast >= 0) {
			ValueT tmp = ddot(ph.alpha, barvb[m]);
			eres.llf += gdatlast * log(tmp) - lfact<ValueT,IntegerT>(gdatlast);
			nn += gdatlast;
			uu += tmp;
      wg(m+1) = gdatlast / tmp;
			daxpy(wg(m+1), barvb[m], eres.eb);
			daxpy(wg(m+1), barvf[m], eres.ey);
		}

    // compupte weights for unobserved periods
		for (size_type k=1; k<=m; k++) {
			if (gdat(k) == -1) {
        tildevf = barvf[k-1];
        tildevf -= barvf[k];
        tildevb = barvb[k-1];
        tildevb -= barvb[k];
        wg(k) = omega;
				daxpy(wg(k), tildevb, eres.eb);
				daxpy(wg(k), tildevf, eres.ey);
			}
		}
		if (gdatlast == -1) {
      wg(m+1) = omega;
			daxpy(wg(m+1), barvb[m], eres.eb);
			daxpy(wg(m+1), barvf[m], eres.ey);
		}
		eres.llf += nn * log(omega) - omega * uu;

    // compute vectors for convolution
		vc[m] = 0;
		daxpy(wg(m+1)-wg(m), ph.baralpha, vc[m]);
		if (idat(m) == 1) {
			daxpy(wp(m), ph.alpha, vc[m]);
		}
		for (size_type k=m-1; k>=1; k--) {
      pois.set(ph.qv * tdat(k+1), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k+1), ph.pdf_eps) + 1);
      mexp_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, ph.P, ph.qv, pois, vc[k+1], vc[k], 0);
			daxpy(wg(k+1)-wg(k), ph.baralpha, vc[k]);
			if (idat(k) == 1) {
				daxpy(wp(k), ph.alpha, vc[k]);
			}
		}
		for (size_type k=1; k<=m; k++) {
      pois.set(ph.qv * tdat(k), 0, poisson<ValueT,RangeT>::rightbound(ph.qv * tdat(k), ph.pdf_eps) + 1);
			dger(NoTrans, wg(k+1)-wg(k), ph.baralpha, barvb[k], eres.en);
      mexpconv_unif<ValueT,RangeT,MatrixT,VectorT>(Trans, NoTrans, ph.P, ph.qv, pois, vc[k], vb[k-1], vb[k-1], eres.en);
		}
		dger(NoTrans, wg(1), ph.baralpha, barvb[0], eres.en);

    eres.etotal = nn + omega * (1 - uu);
		eres.eb *= ph.alpha;
		eres.ez = eres.en.diag();
		eres.en *= ph.T;
    eres.en.diag() = 0;
		eres.ey *= ph.xi;

    for (size_type k=0; k<=m; k++) {
      delete barvf.ptr(k);
      delete barvb.ptr(k);
      delete vb.ptr(k);
      delete vc.ptr(k);
    }
    delete [] pool;

    return eres;
	}

  ////////////////
  // M step
  ////////////////

  template <typename VectorT, typename MatrixT>
  phparams<VectorT,MatrixT>& phase_mstep_gen(
    const phase_estep_results<VectorT,MatrixT>& eres,
    phparams<VectorT,MatrixT>& ph) {

    ph.alpha = eres.eb;
    ph.alpha /= eres.etotal;

    ph.xi = eres.ey;
    ph.xi /= eres.ez;

    ph.T = eres.en;
    ph.T /= eres.ez;
    ph.T.diag() = 0;

    VectorT tmp = ph.T.row_sum();
    tmp += ph.xi;
    ph.T.diag() -= tmp;

    return ph;
  }

  template <typename ValueT, typename RangeT>
  void cf1_swap(const RangeT i, const RangeT j, vector<ValueT,RangeT>& alpha,
    vector<ValueT*,RangeT>& rate, vector<ValueT*,RangeT>& diag) {
    assert(alpha.begin() == rate.begin() && alpha.begin() == diag.begin());
    ValueT tmp;
    ValueT w = rate(j) / rate(i);
    alpha(i) += (1 - w) * alpha(j);
    alpha(j) *= w;
    tmp = rate(j);
    rate(j) = rate(i);
    rate(i) = tmp;
    diag(i) = -rate(i);
    diag(j) = -rate(j);
  }

  template <typename ValueT, typename RangeT>
  void cf1_sort(vector<ValueT,RangeT> alpha,
    vector<ValueT*,RangeT> rate, vector<ValueT*,RangeT> diag) {
    assert(alpha.size() == rate.size() && alpha.size() == diag.size());
    size_type n = alpha.size();
    alpha.set_range(range<RangeT>(1,n));
    rate.set_range(range<RangeT>(1,n));
    diag.set_range(range<RangeT>(1,n));
    for (RangeT i=1; i<=n-1; i++) {
			if (rate(i) > rate(i+1)) {
				cf1_swap(i, i+1, alpha, rate, diag);
				for (RangeT j=i; j>=2; j--) {
					if (rate(j-1) < rate(j)) {
						break;
					}
					cf1_swap(j-1, j, alpha, rate, diag);
				}
			}
		}
  }

  template <typename VectorT, typename MatrixT>
  void bidiag2cf1(phparams<VectorT,MatrixT>& ph) {

    using ValueT = typename VectorT::ValueType;
    using RangeT = typename VectorT::RangeType;

    vector<ValueT*,RangeT> rate(ph.dim());
    vector<ValueT*,RangeT> diag(ph.dim());
    ph.T.diag(diag, 0, diag.begin());
    ph.T.diag(rate, -1, rate.begin());
    rate.ptr(rate.end()) = &ph.xi(ph.xi.end());
    cf1_sort(ph.alpha, rate, diag);
  }

  // dense

  template struct phparams<vector<double,int>,dense_matrix<double,int>>;
  template struct phparams_unif<vector<double,int>,dense_matrix<double,int>>;
  template struct gph_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>>;
  template struct cf1_unif<vector<double,int>,dense_matrix<double,int>,vector<int,int>>;

  // csr

  template struct phparams<vector<double,int>,csr_matrix<double,int>>;
  template struct phparams_unif<vector<double,int>,csr_matrix<double,int>>;
  template struct gph_unif<vector<double,int>,csr_matrix<double,int>,vector<int,int>>;
  template struct cf1_unif<vector<double,int>,csr_matrix<double,int>,vector<int,int>>;

}
