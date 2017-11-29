#ifndef __ITIA_DYNAMICS_KIN_SOLVER_MOVEIT__
#define __ITIA_DYNAMICS_KIN_SOLVER_MOVEIT__
# include <itia_dynamics_core/itia_primitives.h>
# include <ros/ros.h>

# include <moveit/robot_model/robot_model.h>
# include <moveit/robot_state/robot_state.h>
# include <itia_dynamics_core/itia_kin_solver_base.h>


namespace itia
{
namespace dynamics 
{
  class KinSolverMoveit : public KinSolverBase
  {
  public:
    KinSolverMoveit(moveit::core::RobotModelPtr& robot_model);
    virtual itia::CMotion fkine(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame) {};
    virtual itia::JMotion ikine(const itia::CMotion& cmotion, const std::string& tool_frame, const std::string& base_frame,  itia::JMotion jmotion_0, int* result) {};
    
  protected:
    moveit::core::RobotModelPtr& m_robot_model;
    
    
    moveit::core::JointModelGroup* m_jmg;
//     kinematics::KinematicsBase m_kb;
    robot_state::RobotStatePtr m_kinematic_state;
  };
  
}
}






# endif