
#include "utility.hpp"
#include "particle_systems.hpp"
#include "integrators.hpp"

void eulerStep(ParticleSystem& ps, float step) {
	// YOUR CODE HERE (R1)
	// Implement an Euler integrator.
	State init;
	State final;
	State f;

	init = ps.state();
	f = ps.evalF(init);
	for (int i = 0; i < init.size(); i++)
	{
		final.push_back(init[i] + step * f[i]);
	}
	ps.set_state(final);
};

void trapezoidStep(ParticleSystem& ps, float step) {
	// YOUR CODE HERE (R3)
	// Implement a trapezoid integrator.
	State initial;
	State final;
	State final2;
	State f;
	State f2;

	initial = ps.state();
	f = ps.evalF(initial);
	for (int i = 0; i < initial.size(); i++)
	{
		final.push_back(initial[i] + step * f[i]);

	}

	ps.set_state(final);
	f2 = ps.evalF(ps.state());
	for (int i = 0; i < initial.size(); i++)
	{
		final2.push_back(initial[i] + (step / 2) * (f[i] + f2[i]));

	}

	ps.set_state(final2);
}

void midpointStep(ParticleSystem& ps, float step) {
	const auto& x0 = ps.state();
	auto n = x0.size();
	auto f0 = ps.evalF(x0);
	auto xm = State(n), x1 = State(n);
	for (auto i = 0u; i < n; ++i) {
		xm[i] = x0[i] + (0.5f * step) * f0[i];
	}
	auto fm = ps.evalF(xm);
	for (auto i = 0u; i < n; ++i) {
		x1[i] = x0[i] + step * fm[i];
	}
	ps.set_state(x1);
}

void rk4Step(ParticleSystem& ps, float step) {
	// EXTRA: Implement the RK4 Runge-Kutta integrator.
	State initial,s1,s2,s3;
	State final;
	State f;
	State k1,k2,k3,k4;
	initial = ps.state();
	s1 = ps.state();
	s2 = ps.state();
	s3 = ps.state();
	final = ps.state();
	f = ps.state();
	k1 = ps.state();
	k2 = ps.state();
	k3 = ps.state();
	k4 = ps.state();
	f = ps.evalF(initial);
	auto n = initial.size();
	for (auto i = 0u; i < n; ++i) {
		k1[i] = step * f[i];
		s1[i] = initial[i] + k1[i] * 0.5;
	}
	auto s1f = ps.evalF(s1);
	for (auto i = 0u; i < n; ++i) {
		k2[i] = step * s1f[i];
		s2[i] = initial[i] + k2[i] * 0.5;
	}
	auto s2f = ps.evalF(s2);
	for (auto i = 0u; i < n; ++i) {
		k3[i] = step * s2f[i];
		s3[i] = initial[i] + k3[i];
	}
	auto s3f = ps.evalF(s3);
	for (auto i = 0u; i < n; ++i) {
		k4[i] = step * s3f[i];
		final[i] = initial[i] + (k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) / 6;
	}
	ps.set_state(final);
}

#ifdef EIGEN_SPARSECORE_MODULE_H

void implicit_euler_step(ParticleSystem& ps, float step, SparseMatrix& J, SparseLU& solver, bool initial) {
	// EXTRA: Implement the implicit Euler integrator. (Note that the related formula on page 134 on the lecture slides is missing a 'h'; the formula should be (I-h*Jf(Yi))DY=-F(Yi))
}

void implicit_midpoint_step(ParticleSystem& ps, float step, SparseMatrix& J, SparseLU& solver, bool initial) {
	// EXTRA: Implement the implicit midpoint integrator.
}
#endif
