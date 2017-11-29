
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

#ifndef __ITIA_MUTILS_ITIA_PRIVITIES__
#define __ITIA_MUTILS_ITIA_PRIVITIES__

#include <assert.h>

# include <Eigen/Geometry>
# include<Eigen/StdVector>

# include <boost/concept_check.hpp>
# include <boost/enable_shared_from_this.hpp>
# include <boost/graph/graph_concepts.hpp>

# include <urdf/model.h>
# include <ros/console.h>

# include <itia_mutils/frame_distance.h>
# include <itia_dynamics_core/ideal_spring.h>
# include <itia_dynamics_core/friction_polynomial2.h>
# include <itia_dynamics_core/friction_polynomial1.h>
# include <itia_dynamics_core/itia_urdf_parser.h>
# include <itia_dynamics_core/itia_spacevect_algebra.h>

namespace itia 
{
  
using Twist = Eigen::Matrix<double, 6, -1> ;
using Wrench = Eigen::Matrix<double, 6, -1>;
using CMotion = Eigen::Matrix<double, 7, -1>;
using JMotion = Eigen::MatrixXd;


namespace dynamics 
{
  class Link;
  
  class Joint: public boost::enable_shared_from_this<itia::dynamics::Joint>
  {
  protected:
    enum
    {
      REVOLUTE, PRISMATIC, FIXED
    } m_type;
    
    
    
    Eigen::Affine3d m_T_pj; // parent <- joint
    Eigen::Affine3d m_last_T_pc;                            // parent <- child
    Eigen::Matrix3d m_last_R_jc;
    
    Eigen::Vector6d m_screw_of_c_in_p;                      // skew of child origin in parent
    
    Eigen::Vector3d m_axis_in_j;
    Eigen::Matrix3d m_skew_axis_in_j;
    Eigen::Matrix3d m_square_skew_axis_in_j;
    Eigen::Matrix3d m_identity;
    
    Eigen::Vector3d m_axis_in_p;
    Eigen::Matrix3d m_skew_axis_in_p;
    Eigen::Matrix3d m_square_skew_axis_in_p;
    Eigen::Matrix3d m_R_pj;
    
    
    double m_last_q;                                        // last value of q
    
    std::string m_name;
    boost::shared_ptr<itia::dynamics::Link> m_parent_link;
    boost::shared_ptr<itia::dynamics::Link> m_child_link;
  
    void computedTpc();
    void computeJacobian();
  public:
    Joint();
    void fromUrdf(const boost::shared_ptr<urdf::Joint>& urdf_joint, const boost::shared_ptr<itia::dynamics::Link>& parent_link, const boost::shared_ptr<urdf::Link>& child_link);
    boost::shared_ptr<itia::dynamics::Joint> pointer();
    std::string getName(){return m_name;};
    boost::shared_ptr <itia::dynamics::Link> getChildLink()  {return m_child_link;}; 
    boost::shared_ptr <itia::dynamics::Link> getParentLink() {return m_parent_link;}; 
    const Eigen::Affine3d& getTransformation(   const double& q = 0);
    const Eigen::Vector6d& getScrew_of_child_in_parent();                      
    
    const bool isFixed() {return (m_type == FIXED);};
  };
  

  
  class Link: public boost::enable_shared_from_this<itia::dynamics::Link>
  {
  protected:
    boost::shared_ptr<itia::dynamics::Joint> m_parent_joint;
    std::vector<boost::shared_ptr<itia::dynamics::Joint>> m_child_joints;                     
    std::vector<boost::shared_ptr<itia::dynamics::Link>> m_child_links;                     
    std::string m_name;
    
    double m_mass;
    Eigen::Vector3d m_cog_in_c;
    Eigen::Matrix<double, 6, 6> m_Inertia_cc;
    std::vector< Eigen::Matrix<double, 6, 6>, Eigen::aligned_allocator<Eigen::Matrix<double, 6, 6>> > m_Inertia_cc_single_term;
    
    
  public:
    
    Link() {};
    void fromUrdf(const boost::shared_ptr<urdf::Link>& urdf_link, 
                  const boost::shared_ptr<itia::dynamics::Joint>& parent_joint = 0);
    boost::shared_ptr<itia::dynamics::Link> pointer();
    
    std::string getName(){return m_name;};
    boost::shared_ptr<itia::dynamics::Joint> getParentJoint() {return m_parent_joint;};
    std::vector<boost::shared_ptr<itia::dynamics::Joint>> getChildrenJoints() {return m_child_joints;};
    boost::shared_ptr<itia::dynamics::Link> findChild(const std::string& name);
    boost::shared_ptr<itia::dynamics::Joint> findChildJoint(const std::string& name);
    const Eigen::Matrix66d& getSpatialInertia() {return m_Inertia_cc;};
    const std::vector< Eigen::Matrix<double, 6, 6>, Eigen::aligned_allocator<Eigen::Matrix<double, 6, 6>> >& getSpatialInertiaTerms() {return m_Inertia_cc_single_term;};
    const double& getMass() {return m_mass;};
    const Eigen::Vector3d& getCog() {return m_cog_in_c;};
    
  };
    
    
  class Chain
  {
  protected:
    std::vector<boost::shared_ptr<itia::dynamics::Link>> m_links;
    std::vector<boost::shared_ptr<itia::dynamics::Joint>> m_joints;
    unsigned int m_joints_number;
    unsigned int m_active_joints_number;
    unsigned int m_links_number;
    
    std::vector<std::string> m_links_name;
    std::map<std::string, unsigned int> m_joints_name;
    
    Eigen::Matrix6Xd m_jacobian;
    
    Eigen::Affine3d m_T_bt;                               // base <- tool
    
    Eigen::VectorXd m_last_q;
    Eigen::VectorXd m_sorted_q;

    Eigen::VectorXd m_last_Dq;
    Eigen::VectorXd m_sorted_Dq;
    
    Eigen::VectorXd m_last_DDq;
    Eigen::VectorXd m_sorted_DDq;

    Eigen::VectorXd m_last_DDDq;
    Eigen::VectorXd m_sorted_DDDq;

    
    Eigen::VectorXd m_joint_torques;
    Eigen::VectorXd m_active_joint_torques;
    Eigen::MatrixXd m_regressor_extended;
    
    Eigen::MatrixXd m_input_to_chain_joint;
    Eigen::MatrixXd m_chain_to_input_joint;
    std::vector<unsigned int> m_active_joints;

    std::vector<Eigen::Affine3d,Eigen::aligned_allocator<Eigen::Affine3d> > m_T_bl;  //
    bool m_is_screws_computed;
    bool m_is_jac_computed;
    bool m_is_vel_computed;
    bool m_is_acc_computed;
    bool m_is_linacc_computed;
    bool m_is_nonlinacc_computed;
    bool m_is_jerk_computed;
    bool m_is_linjerk_computed;
    bool m_is_nonlinjerk_computed;
    bool m_is_wrench_computed;
    bool m_is_regressor_computed;
    
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_screws_of_c_in_b;
    
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_Ds;
    
    // twists of c in b
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_twists;
    // Dtwists of c in b
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_Dtwists;
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_Dtwists_linear_part;
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_Dtwists_nonlinear_part;
    
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_DDtwists;
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_DDtwists_linear_part;
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_DDtwists_nonlinear_part;
    
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_wrenches;
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_inertial_wrenches;
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> m_gravity_wrenches;
    std::vector<Eigen::Matrix<double, 6, 10>,Eigen::aligned_allocator<Eigen::Matrix<double, 6, 10>>> m_wrenches_regressor;
    
    Eigen::Vector3d m_gravity;
    Eigen::MatrixXd m_joint_inertia;
    Eigen::MatrixXd m_joint_inertia_extended;
  
    void computeScrews();
    void computeFrames();
    
  public:
    Chain(const boost::shared_ptr<itia::dynamics::Link>& root_link, const std::string& base_link_name, const std::string& ee_link_name, const Eigen::Vector3d gravity = Eigen::Vector3d::Zero());
    void setInputJointsName(const std::vector<std::string> joints_name);
    unsigned int getLinksNumber() {return m_links_number;};
    unsigned int getJointsNumber() {return m_joints_number;};
    unsigned int getActiveJointsNumber() {return m_active_joints_number;};
    std::vector<std::string> getLinksName() {return m_links_name;};
    
    /*
     * Kinematics methods
     */
    Eigen::Affine3d getTransformation(const Eigen::VectorXd& q);
    Eigen::Matrix6Xd getJacobian(const Eigen::VectorXd& q);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getTwist(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getDTwist(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getDTwistLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& DDq);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getDTwistNonLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getDDTwistLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& DDDq);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getDDTwistNonLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq);
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getDDTwist(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq, const Eigen::VectorXd& DDDq);
    
    /*
     * Dynamics methods
     */
    std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> getWrench(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq, std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > ext_wrenches_in_link_frame);
    Eigen::VectorXd getJointTorque(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq, std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > ext_wrenches_in_link_frame);
    Eigen::VectorXd getJointTorque(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq);
    
    Eigen::MatrixXd getRegressor( const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& Dq, 
                                  const Eigen::VectorXd& DDq);
    
    Eigen::MatrixXd getJointInertia(const Eigen::VectorXd& q);
    
  };
}
}
# endif



/*
 * W_ii=[S_ii, vi;0 0]
 * T_ji=[R_ji, p_ji;0 1]
 * W_jj=T_ji*W_ii*T_ij
 * 
 * T_ij=[R_ji', -R_ji'*p_ji;0 1]
 * W_ii*T_ij = [S_ii*R_ji', -S_ii*R_ji'*p_ji+vi;0 0]
 * Rji*W_ii*T_ij = [R_ji*S_ii*R_ji', -R_ji*S_ii*R_ji'*p_ji+R_ji*vi;0 0]
 * Rji*W_ii*T_ij = [S_jj, -S_jj*p_ji+R_ji*vi;0 0]
 */