#ifndef __ITIA_DYNAMICS_KIN_SOLVER_KDL__
#define __ITIA_DYNAMICS_KIN_SOLVER_KDL__

# include <itia_dynamics_core/itia_kin_solver_base.h>
# include <itia_dynamics_core/kdl_conversions.h>
# include <kdl/chainfksolvervel_recursive.hpp>
# include <kdl/chainfksolverpos_recursive.hpp>
# include <kdl/chainiksolverpos_lma.hpp>
# include <kdl/chainjnttojacsolver.hpp>
# include <kdl/chainiksolvervel_pinv_nso.hpp>
# include <kdl/chainidsolver_recursive_newton_euler.hpp>
# include <boost/shared_array.hpp>
# include <map>
# include <kdl_parser/kdl_parser.hpp>

//  TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
namespace itia
{
namespace dynamics 
{
  class KinSolverEigen: public KinSolverBase
  {
  public:
    KinSolverKDL(const std::string& robot_description);
    itia::CMotion fkine(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame);
    itia::JMotion ikine(const itia::CMotion& cmotion, const std::string& tool_frame, const std::string& base_frame,  itia::JMotion jmotion_0, int* result);
    virtual Eigen::Matrix<double, 6, -1> jac(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame);
    virtual Eigen::VectorXd getTorques(const itia::JMotion& jmotion, const itia::Wrench& ext_wrench, const std::string& tool_frame, const std::string& base_frame);
    virtual Eigen::VectorXd getTorques(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame);
    
  protected:
    KDL::Tree m_tree;
    std::string m_robot_description;
    std::map<std::vector<std::string>, boost::shared_ptr<KDL::ChainFkSolverPos_recursive>> m_pos_fkine;
    std::map<std::vector<std::string>, boost::shared_ptr<KDL::ChainFkSolverVel_recursive>> m_vel_fkine;

    std::map<std::vector<std::string>, boost::shared_ptr<KDL::ChainIdSolver_RNE>> m_dyn_solver;
    
    std::map<std::vector<std::string>, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>> m_pos_ikine;
    std::map<std::vector<std::string>, boost::shared_ptr<KDL::ChainJntToJacSolver>> m_jac_solver;
    std::map<std::vector<std::string>, boost::shared_ptr<KDL::Chain>> m_chain;
    

  };
  
}
}






# endif