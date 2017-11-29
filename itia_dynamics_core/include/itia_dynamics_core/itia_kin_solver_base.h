#ifndef __ITIA_DYNAMICS_KIN_SOLVER_BASE__
#define __ITIA_DYNAMICS_KIN_SOLVER_BASE__
# include <itia_dynamics_core/itia_primitives.h>
# include <ros/ros.h>

namespace itia
{
namespace dynamics 
{
  class KinSolverBase
  {
  public:
    KinSolverBase() {};
    virtual itia::CMotion fkine(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame) = 0;
    virtual itia::JMotion ikine(const itia::CMotion& cmotion, const std::string& tool_frame, const std::string& base_frame,  itia::JMotion jmotion_0, int* result) = 0;
    virtual Eigen::Matrix<double, 6, -1> jac(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame) = 0;
    virtual unsigned int getJointNr() {return m_joint_nr;};
  protected:
    unsigned int m_joint_nr;
  };
  
//   typedef boost::shared_ptr<itia::dynamics::KinSolverBase> KinSolverBasePtr;
  using KinSolverBasePtr = boost::shared_ptr<itia::dynamics::KinSolverBase> ;
}
}






# endif