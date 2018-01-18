/*
  poisson
 */

 #pragma once

namespace marlib {

	template<typename ValueT, typename RangeT>
	class poisson {
		static constexpr ValueT NORMALQ_LOWER_Q = 3.0;
		static constexpr ValueT NORMALQ_UPPER_Q = 37.0;

		static constexpr ValueT NORMALQ_LOWER_LOGP = -689.0;
		static constexpr ValueT NORMALQ_UPPER_LOGP = -6.6;

		static constexpr ValueT NORMALQ_EPS = 1.0e-8;

		static constexpr ValueT LOG2PIOVER2 = 0.9189385332046727417803297364; // log(2pi) / 2

		static constexpr ValueT POISSON_LAMBDA_MIN = 3.0;
		static constexpr RangeT POISSON_RIGHT_MAX = 23;

	public:
		poisson(const ValueT lambda, const ValueT eps);
    poisson(const RangeT minleft, const RangeT maxright);
		~poisson();

		const RangeT left() const;
		const RangeT right() const;
		array<ValueT>& get_prob();
		const ValueT weight() const;
    const ValueT lambda() const;

		const ValueT operator()(const RangeT i) const;

    static RangeT rightbound(const ValueT lambda, const ValueT eps);

		void set(const ValueT lambda, const ValueT eps);
    void set(const ValueT lambda, const RangeT left, const RangeT right);

	private:
		ValueT m_lambda;
		RangeT m_left;
		RangeT m_right;
		array<ValueT> m_prob;
		ValueT m_weight;

		static ValueT normalt(const ValueT x);
		static ValueT normalq(const ValueT p);
		ValueT set_prob();

	};

  template <typename ValueT, typename RangeT>
	inline const RangeT poisson<ValueT,RangeT>::left() const {
		return m_left;
	}

	template <typename ValueT, typename RangeT>
	inline const RangeT poisson<ValueT,RangeT>::right() const {
		return m_right;
	}

	template <typename ValueT, typename RangeT>
	inline array<ValueT>& poisson<ValueT,RangeT>::get_prob() {
		return m_prob;
	}

	template <typename ValueT, typename RangeT>
	inline const ValueT poisson<ValueT,RangeT>::weight() const {
		return m_weight;
	}

	template <typename ValueT, typename RangeT>
	inline const ValueT poisson<ValueT,RangeT>::lambda() const {
		return m_lambda;
	}

	template <typename ValueT, typename RangeT>
	inline const ValueT poisson<ValueT,RangeT>::operator()(const RangeT i) const {
		return m_prob[i - m_left];
	}
}
