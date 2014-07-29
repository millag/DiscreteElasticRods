#ifndef __ABSTRACTRK4INTEGRATOR_H__
#define __ABSTRACTRK4INTEGRATOR_H__

#include <ngl/Vec3.h>
/// @brief a simple class to hold State values for the integrator, this is used for
/// the derivitives as well as the initial state
/// @note this class is based on the tutorials here
/// http://gafferongames.com/game-physics/
class State
{
	public :
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief the position of the state
		//----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 m_position;
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief the velocity of the state
		//----------------------------------------------------------------------------------------------------------------------
    ngl::Vec3 m_velocity;
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief ctor
		/// @brief[in] _pos initial position
		/// @brief[in] _vel initial velocity
		//----------------------------------------------------------------------------------------------------------------------
		inline State( ngl::Vec3 _pos, ngl::Vec3 _vel) :
										m_position(_pos),
										m_velocity(_vel){;}
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief default ctor
		//----------------------------------------------------------------------------------------------------------------------
		inline State(){;}
};

/// @brief a simple base class to do RK4 integration this class then needs the motionFunction to be
/// implemented so we can integrate the values

class AbstractRK4Integrator
{
	public :

	//----------------------------------------------------------------------------------------------------------------------
	/// @brief ctor passing in inital state
	/// @param[in] _state the state to set for the initial value
	//----------------------------------------------------------------------------------------------------------------------
	AbstractRK4Integrator(State _state ) :
												  m_state(_state){;}
	//----------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------
	AbstractRK4Integrator(){;}
	//----------------------------------------------------------------------------------------------------------------------
		virtual ~AbstractRK4Integrator(){;}

		//----------------------------------------------------------------------------------------------------------------------
		/// @brief the actual function we need to implement for the integration this can be anything
		/// this is a pure virtual method an needs to be implemented in the class inheriting this
	  /// @param[in] _state the state we wish to calculate the motion function for
		/// @param[in] _t the current time step to integrate for
		//----------------------------------------------------------------------------------------------------------------------
		virtual ngl::Vec3 motionFunction(const State &_state,float _t)=0;
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief method to evaluate the state given initial value and time step
		/// @param[in] _initial the initial state
		/// @param[in] _t the current time step value
		//----------------------------------------------------------------------------------------------------------------------
		State evaluate(const State &_initial, float _t);
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief method to evaluate the state given initial value and time step
		/// @param[in] _initial the initial state
		/// @param[in] _t the current time step value
		/// @param[in] _dt the delta time step for integration
		/// @param[in] _d the derivative for the integration step
		//----------------------------------------------------------------------------------------------------------------------
		State evaluate(const State &_initial,float _t, float _dt, const State &_d);
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief method to evaluate the state given initial value and time step
		/// @param[in] _initial the initial state
		/// @param[in] _t the current time step value
		/// @param[in] _dt the delta time step for integration
		//----------------------------------------------------------------------------------------------------------------------
		void integrate(float _t,	 float _dt);

		//----------------------------------------------------------------------------------------------------------------------
		/// @brief accesor to get the current state value
		/// @returns m_state the current state value
		//----------------------------------------------------------------------------------------------------------------------
		inline State getState() const {return m_state;}
	protected :
		//----------------------------------------------------------------------------------------------------------------------
		/// @brief the state value
		//----------------------------------------------------------------------------------------------------------------------
		State m_state;


};

#endif
