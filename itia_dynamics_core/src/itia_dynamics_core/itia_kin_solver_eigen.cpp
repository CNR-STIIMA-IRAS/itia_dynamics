
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

# include <itia_dynamics_core/itia_kin_solver_kdl.h>
#include <kdl/framevel_io.hpp>

// TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
#include <itia_dynamics_core/itia_urdf_parser.h>
namespace itia
{
namespace dynamics 
{
  
itia::dynamics::KinSolverKDL::KinSolverKDL(const std::string& robot_description):
KinSolverBase()
{
  m_robot_description = robot_description;
  if ( !kdl_parser::treeFromParam( robot_description, m_tree ) )
  {
    ROS_ERROR("Error on extracting tree from param");
    throw ("Error on extracting tree from param");
  }
}

CMotion KinSolverKDL::fkine(const JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame)
{
  std::vector<std::string> frames_couple(2);
  
  KDL::Chain chain;
  if (!m_tree.getChain(base_frame,tool_frame,chain))
  {
    ROS_ERROR("Error on extracting chain from tree");
    throw ("Error on extracting chain from tree");
  }
  
  frames_couple.at(0) = base_frame;
  frames_couple.at(1) = tool_frame;
  
  itia::CMotion cmotion;
  cmotion.resize(7, jmotion.cols());
  cmotion.setZero();
  
  if (jmotion.cols() == 1)
  {
    KDL::Frame frame;
    KDL::JntArray joint;
    itia::dynamics::JMotionToJoint(jmotion, &joint);
    try
    {
      m_pos_fkine.at(frames_couple)->JntToCart(joint, frame);
    }
    catch (std::out_of_range)
    {
      m_pos_fkine.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainFkSolverPos_recursive>>(  frames_couple, boost::shared_ptr<KDL::ChainFkSolverPos_recursive>(new KDL::ChainFkSolverPos_recursive(chain))) );
      m_pos_fkine.at(frames_couple)->JntToCart(joint, frame);
    }
    itia::dynamics::frameToCMotion(frame, &cmotion);
  }
  else
  {
    KDL::FrameVel frame_vel;
    KDL::JntArrayVel joint_vel;
    itia::dynamics::JMotionToJoint(jmotion, &joint_vel);
    try
    {
      m_vel_fkine.at(frames_couple)->JntToCart(joint_vel, frame_vel);
    }
    catch (std::out_of_range)
    {
      m_vel_fkine.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainFkSolverVel_recursive>>(  frames_couple, boost::shared_ptr<KDL::ChainFkSolverVel_recursive>(new KDL::ChainFkSolverVel_recursive(chain))) );
      m_vel_fkine.at(frames_couple)->JntToCart(joint_vel, frame_vel);
    }
    itia::dynamics::frameVelToCMotion(frame_vel, &cmotion);
  }
  return cmotion;
}

Eigen::Matrix<double, 6, -1> KinSolverKDL::jac(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame)
{
  std::vector<std::string> frames_couple(2);
  frames_couple.at(0) = base_frame;
  frames_couple.at(1) = tool_frame;
  
  if (m_chain.count(frames_couple) == 0)
  {
    m_chain.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::Chain>>(  frames_couple, boost::shared_ptr<KDL::Chain>(new KDL::Chain)));
    m_tree.getChain(base_frame,tool_frame,*m_chain.at(frames_couple));
    m_pos_ikine.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>>(  frames_couple, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>(new KDL::ChainIkSolverPos_LMA(*m_chain.at(frames_couple)))) );
    m_jac_solver.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainJntToJacSolver>>(  frames_couple, boost::shared_ptr<KDL::ChainJntToJacSolver>(new KDL::ChainJntToJacSolver(*m_chain.at(frames_couple)))) );
  }
  int n_jnt = m_chain.at(frames_couple)->getNrOfJoints();
  KDL::JntArray joint(n_jnt);
  itia::dynamics::JMotionToJoint(jmotion, &joint);
  KDL::Jacobian jac(n_jnt);
  m_jac_solver.at(frames_couple)->JntToJac(joint, jac);
  return jac.data;
}


JMotion KinSolverKDL::ikine(const CMotion& cmotion, const std::string& tool_frame, const std::string& base_frame,  itia::JMotion jmotion_0, int* result)
{
  std::vector<std::string> frames_couple(2);
  frames_couple.at(0) = base_frame;
  frames_couple.at(1) = tool_frame;
  
  if (m_chain.count(frames_couple) == 0)
  {
    m_chain.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::Chain>>(  frames_couple, boost::shared_ptr<KDL::Chain>(new KDL::Chain)));
    m_tree.getChain(base_frame,tool_frame,*m_chain.at(frames_couple));
    m_pos_ikine.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>>(  frames_couple, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>(new KDL::ChainIkSolverPos_LMA(*m_chain.at(frames_couple)))) );
    m_jac_solver.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainJntToJacSolver>>(  frames_couple, boost::shared_ptr<KDL::ChainJntToJacSolver>(new KDL::ChainJntToJacSolver(*m_chain.at(frames_couple)))) );
  }
  
  int n_jnt = m_chain.at(frames_couple)->getNrOfJoints();
  itia::JMotion jmotion;
    
  KDL::Frame frame;
  KDL::JntArray joint(n_jnt);
  KDL::JntArray joint_0(n_jnt);
  
  itia::dynamics::CMotionToFrame(cmotion, &frame);
  itia::dynamics::JMotionToJoint(jmotion_0, &joint_0);
  *result = m_pos_ikine.at(frames_couple)->CartToJnt(joint_0, frame, joint);
  
  if (cmotion.cols() == 1)
  {
    itia::dynamics::JointToJMotion(joint, &jmotion);
  }
  else
  {
    KDL::JntArrayVel joint_vel(n_jnt);
    joint_vel.q = joint;
    
    itia::Twist twist;
    itia::dynamics::CMotionToAffine(cmotion, &twist);
    
    KDL::Jacobian jac(n_jnt);
    m_jac_solver.at(frames_couple)->JntToJac(joint, jac);
    
    Eigen::FullPivLU<Eigen::MatrixXd> jac_decomp(jac.data);
    Eigen::VectorXd dq = jac_decomp.solve(twist.col(0));
    if (jac_decomp.rank()<dq.rows())
    {
      Eigen::MatrixXd null=jac_decomp.kernel();
      ROS_WARN_THROTTLE(0.1,"Singolarity point!");
      
      for (int iC=0;iC<null.cols();iC++)
      {
        Eigen::VectorXd null_versor=null.col(iC);
        null_versor.normalize();
        dq=dq-(dq.dot(null_versor))*null_versor;
      }
    }
    joint_vel.qdot.data = dq;
    itia::dynamics::JointToJMotion(joint_vel, &jmotion);
  }
  return jmotion;
}

Eigen::VectorXd KinSolverKDL::getTorques(const itia::JMotion& jmotion, const Wrench& ext_wrench, const std::string& tool_frame, const std::string& base_frame)
{
  std::vector<std::string> frames_couple(2);
  
  frames_couple.at(0) = base_frame;
  frames_couple.at(1) = tool_frame;
  int nj, ns;
  
  if (m_chain.count(frames_couple) == 0)
  {
    m_chain.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::Chain>>(  frames_couple, boost::shared_ptr<KDL::Chain>(new KDL::Chain)));
    m_tree.getChain(base_frame,tool_frame,*m_chain.at(frames_couple));
    KDL::Vector grav(0, 0, -9.806);
    m_dyn_solver.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainIdSolver_RNE>>(  frames_couple, boost::shared_ptr<KDL::ChainIdSolver_RNE>(new KDL::ChainIdSolver_RNE(*m_chain.at(frames_couple), grav))) );
    m_pos_ikine.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>>(  frames_couple, boost::shared_ptr<KDL::ChainIkSolverPos_LMA>(new KDL::ChainIkSolverPos_LMA(*m_chain.at(frames_couple)))) );
    m_jac_solver.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainJntToJacSolver>>(  frames_couple, boost::shared_ptr<KDL::ChainJntToJacSolver>(new KDL::ChainJntToJacSolver(*m_chain.at(frames_couple)))) );
  }
  
  if (m_dyn_solver.count(frames_couple) == 0)
  {
    KDL::Vector grav(0, 0, -9.806);
    m_dyn_solver.insert(std::pair<std::vector<std::string>, boost::shared_ptr<KDL::ChainIdSolver_RNE>>(  frames_couple, boost::shared_ptr<KDL::ChainIdSolver_RNE>(new KDL::ChainIdSolver_RNE(*m_chain.at(frames_couple), grav))) );
    
  }
  
  ns = m_chain.at( frames_couple )->getNrOfSegments();
  nj = m_chain.at( frames_couple )->getNrOfJoints();
  
  KDL::JntArray q;
  KDL::JntArray Dq;
  KDL::JntArray DDq;
  KDL::JntArray tau(nj);
  
  
  q.data = jmotion.col(0);
  if (jmotion.cols()>1)
    Dq.data = jmotion.col(1);
  else
    Dq.data =  Eigen::VectorXd::Constant(q.data.rows(), 1, 0.0);
  if (jmotion.cols()>2)
    DDq.data = jmotion.col(2);
  else
    DDq.data = Eigen::VectorXd::Constant(q.data.rows(), 1, 0.0);
    
  KDL::Wrenches kdl_ext_wrenches(ns);
  
  itia::dynamics::WrenchToKDLWrench(ext_wrench, &(kdl_ext_wrenches.at(ns-1)));
  if (m_dyn_solver.at(frames_couple)->CartToJnt(q, Dq, DDq, kdl_ext_wrenches, tau))
    ROS_ERROR_THROTTLE(0.1, "Error computing forces!");
  
  return tau.data;  
}

Eigen::VectorXd KinSolverKDL::getTorques(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame)
{
  return getTorques(jmotion, Eigen::VectorXd::Constant(6, 1, 0), tool_frame, base_frame);
}

 
}
}


