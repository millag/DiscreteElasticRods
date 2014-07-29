#include "AbstractRK4Integrator.h"

//----------------------------------------------------------------------------------------------------------------------
State AbstractRK4Integrator::evaluate(const State &_initial, float _t	)
{

	State output;
	output.m_position=_initial.m_velocity;
	output.m_velocity=motionFunction(_initial, _t);
	return output;
}
//----------------------------------------------------------------------------------------------------------------------
State AbstractRK4Integrator::evaluate( const State &_initial, float _t,float _dt, const State &_d	)
{

	State state;
	state.m_position = _initial.m_position + _d.m_position*_dt;
	state.m_velocity = _initial.m_velocity + _d.m_velocity*_dt;
	State output;
	output.m_position = state.m_velocity;
	output.m_velocity = motionFunction(state, _t+_dt);
	return output;
}
//----------------------------------------------------------------------------------------------------------------------
void AbstractRK4Integrator::integrate(float _t, float _dt	)
{
	// basic RK4 integration we calculate the initial value then two slopes
	// and combine together see
	// http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods for more details
	State a = evaluate(m_state, _t);
	State b = evaluate(m_state, _t, _dt*0.5f, a);
	State c = evaluate(m_state, _t, _dt*0.5f, b);
	State d = evaluate(m_state, _t, _dt, c);
	const ngl::Vec3 dxdt = 1.0f/6.0f * (a.m_position + 2.0f*(b.m_position + c.m_position) + d.m_position);
	const ngl::Vec3 dvdt = 1.0f/6.0f * (a.m_velocity + 2.0f*(b.m_velocity + c.m_velocity) + d.m_velocity);

	m_state.m_position = m_state.m_position + dxdt*_dt;
	m_state.m_velocity = m_state.m_velocity + dvdt*_dt;
}
