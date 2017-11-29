
 // -------------------------------------------------------------------------------- 
 // Copyright (c) 2017 CNR-ITIA <iras@itia.cnr.it>
 // All rights reserved.
 //
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are met:
 //
 // 1. Redistributions of source code must retain the above copyright notice,
 // this list of conditions and the following disclaimer.
 // 2. Redistributions in binary form must reproduce the above copyright
 // notice, this list of conditions and the following disclaimer in the
 // documentation and/or other materials provided with the distribution.
 // 3. Neither the name of mosquitto nor the names of its
 // contributors may be used to endorse or promote products derived from
 // this software without specific prior written permission.
 //
 //
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 // ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 // LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 // CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 // SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 // INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 // CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 // ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 // POSSIBILITY OF SUCH DAMAGE.
 // -------------------------------------------------------------------------------- 

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