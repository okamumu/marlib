/*
poisson
*/

#include "marlib.hpp"

namespace marlib {

	template <typename ValueT, typename RangeT>
	poisson<ValueT,RangeT>::poisson(const ValueT lambda, const ValueT eps)
	: m_lambda(lambda), m_left(0), m_right(rightbound(lambda,eps)), m_prob(m_right+1), m_weight(set_prob()) {}

	template <typename ValueT, typename RangeT>
	poisson<ValueT,RangeT>::poisson(const RangeT minleft, const RangeT maxright)
	: m_left(minleft), m_right(maxright), m_prob(maxright-minleft+1) {}

	template <typename ValueT, typename RangeT>
	poisson<ValueT,RangeT>::~poisson() {}

	// template <typename ValueT, typename RangeT>
	// const RangeT poisson<ValueT,RangeT>::left() const {
	// 	return m_left;
	// }
  //
	// template <typename ValueT, typename RangeT>
	// const RangeT poisson<ValueT,RangeT>::right() const {
	// 	return m_right;
	// }
  //
	// template <typename ValueT, typename RangeT>
	// array<ValueT>& poisson<ValueT,RangeT>::get_prob() {
	// 	return m_prob;
	// }
  //
	// template <typename ValueT, typename RangeT>
	// const ValueT poisson<ValueT,RangeT>::weight() const {
	// 	return m_weight;
	// }
  //
	// template <typename ValueT, typename RangeT>
	// const ValueT poisson<ValueT,RangeT>::lambda() const {
	// 	return m_lambda;
	// }
  //
	// template <typename ValueT, typename RangeT>
	// const ValueT poisson<ValueT,RangeT>::operator()(const RangeT i) const {
	// 	return m_prob[i - m_left];
	// }

	template <typename ValueT, typename RangeT>
	void poisson<ValueT,RangeT>::set(const ValueT lambda, const ValueT eps) {
		set(lambda, 0, rightbound(lambda, eps));
	}

	template <typename ValueT, typename RangeT>
	void poisson<ValueT,RangeT>::set(const ValueT lambda, const RangeT left, const RangeT right) {
		// if (m_lambda != lambda) {
			m_lambda = lambda;
			m_left = left;
			m_right = right;
			assert(m_prob.size() >= m_right - m_left + 1);
			m_weight = set_prob();
		// }
	}

	/*
	Description
	return a tail probability of standard normal distribution
	Parameters
	IN
	x: input value
	OUT
	return value
	*/

	template <typename ValueT, typename RangeT>
	ValueT poisson<ValueT,RangeT>::normalt(const ValueT x) {
		ValueT x2 = x*x;
		ValueT tmp = x;
		ValueT sum = 1 / tmp;
		tmp = tmp * x2;
		sum = sum - 1 / tmp;
		tmp = tmp * x2;
		sum = sum + 3 / tmp;
		tmp = tmp * x2;
		sum = sum - 15 / tmp;
		tmp = tmp * x2;
		sum = sum + 105 / tmp;
		return (log(sum) - x2/2 - LOG2PIOVER2);
	}

	/*
	Description
	return a quantile of standard normal distribution
	Parameters
	IN
	p: probability for the quantile
	OUT
	return value
	*/

	template <typename ValueT, typename RangeT>
	ValueT poisson<ValueT,RangeT>::normalq(const ValueT p) {
		const ValueT leps = log(p);
		assert(leps <= NORMALQ_UPPER_LOGP && leps >= NORMALQ_LOWER_LOGP);
		ValueT l = NORMALQ_LOWER_Q;
		ValueT u = NORMALQ_UPPER_Q;
		ValueT m = (l + u) / 2;
		ValueT fm = normalt(m) - leps;
		while (fabs(fm) > NORMALQ_EPS) {
			if (fm > 0) {
				l = m;
			} else {
				u = m;
			}
			m = (l + u)/2;
			fm = normalt(m) - leps;
		}
		return m;
	}

	/*
	! Description: compute the right bound of Poisson range
	!              for a given error tolerance
	!
	! Parameters:
	!   IN
	!    lambda: Poisson rate (mean)
	!    epsi: error tolerance
	!   OUT
	!    right bound is a return value
	*/

	template <typename ValueT, typename RangeT>
	RangeT poisson<ValueT,RangeT>::rightbound(const ValueT lambda, const ValueT eps) {
		RangeT right;
		if (lambda == 0) {
			return 0;
		}
		if (lambda < POISSON_LAMBDA_MIN) {
			ValueT tmp = exp(-lambda);
			ValueT total = tmp;
			right = 0;
			for (RangeT k=1; k<=POISSON_RIGHT_MAX; k++) {
				right++;
				tmp *= lambda / right;
				total += tmp;
				if (total + eps >= 1)
				break;
			}
		} else {
			ValueT z = normalq(eps);
			ValueT tmp = z + sqrt(4 * lambda - 1);
			right = static_cast<RangeT>(tmp * tmp / 4 + 1);
		}
		return right;
	}

	/*
	! Description & parameters :
	!  IN
	!    lambda: Poisson parameter (mean)
	!    left, right: left and right bounds
	!  OUT
	!    prob: Poisson probabilities from left to right bounds
	!    weight: weight value, i.e., exact pmf is prob[i]/weight
	*/

	template <typename ValueT, typename RangeT>
	ValueT poisson<ValueT,RangeT>::set_prob() {
		RangeT mode = static_cast<RangeT>(m_lambda);
		if (mode >= 1) {
			m_prob[mode-m_left] = exp(-m_lambda + mode * log(m_lambda) - LOG2PIOVER2 - (mode + 1.0/2.0) * log(mode) + mode);
		} else {
			m_prob[mode-m_left] = exp(-m_lambda);
		}
		// -- down --
		for (RangeT j=mode; j>=m_left+1; j--) {
			m_prob[j-1-m_left] = m_prob[j-m_left] * j / m_lambda;
		}
		// -- up --
		for (RangeT j=mode; j<=m_right-1; j++) {
			m_prob[j+1-m_left] = m_prob[j-m_left] * m_lambda / (j+1);
		}
		// -- compute W --
		ValueT weight = 0;
		RangeT s = m_left;
		RangeT t = m_right;
		while (s < t) {
			if (m_prob[s-m_left] <= m_prob[t-m_left]) {
				weight += m_prob[s-m_left];
				s++;
			} else {
				weight += m_prob[t-m_left];
				t--;
			}
		}
		weight += m_prob[s-m_left];
		return weight;
	}

	// instance
	template class poisson<double,int>;

}
