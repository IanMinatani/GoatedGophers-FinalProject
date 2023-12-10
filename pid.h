#include <algorithm>
#include <math.h>

class pid {
	private:
		const float kp_;
		const float ki_;
		const float kd_;
		const float f_co_d_Hz_;
		float u_[2];
		float y_[2];

	public:
		pid(const float kp, const float ki, const float kd, const float f_co_d_Hz) :
			kp_(kp),
			ki_(ki),
			kd_(kd),
			f_co_d_Hz_(f_co_d_Hz) {
			u_[0] = 0;
			u_[1] = 0;
			y_[0] = 0;
			y_[1] = 0;
		}
		~pid() {
			// No-op.
		}
		float step(const float u_k, const float T) {
			const float N = 2*3.14159*f_co_d_Hz_;

			const float b0 = kp_*(1+N*T)+ki_*T*(1+N*T)+kd_*N;
			const float b1 = -(kp_*(2+N*T)+ki_*T+2*kd_*N);
			const float b2 = kp_+kd_*N;

			const float a0 = 1+N*T;
			const float a1 = -(2+N*T);
			const float a2 = 1;

			const float y_k = ((b0*u_k + b1*u_[0] + b2*u_[1]) - (a1*y_[0] + a2*y_[1]))/a0;

			u_[1] = u_[0];
			u_[0] = u_k;

			y_[1] = y_[0];
			y_[0] = y_k;

			return y_k;
		}
};
