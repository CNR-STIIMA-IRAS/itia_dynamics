# include <itia_dynamics_core/itia_kin_solver_moveit.h>



namespace itia
{
  namespace dynamics 
{

KinSolverMoveit::KinSolverMoveit(moveit::core::RobotModelPtr& robot_model)
: 
m_robot_model(robot_model), 
KinSolverBase()
{
//   m_jmg =  m_robot_model->getJointModelGroup();
//   auto m_kb = m_jmg->getSolverInstance();
//   m_kb.getPositionFK();
//   
//   m_kinematic_state.reset(new robot_state::RobotState(m_robot_model));
//   m_kinematic_state->getJacobian(m_jmg);
//   m_kinematic_state->getJacobian();
}



}
}
