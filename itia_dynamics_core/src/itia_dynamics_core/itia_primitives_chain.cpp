
#include <exception>
#include <itia_dynamics_core/itia_primitives.h>

namespace itia 
{
namespace dynamics 
{

Chain::Chain(const boost::shared_ptr<itia::dynamics::Link>& root_link, 
             const std::string& base_link_name, 
             const std::string& ee_link_name, 
             Eigen::Vector3d gravity)
{
  m_is_screws_computed = 
  m_is_jac_computed = 
  m_is_vel_computed = 
  m_is_acc_computed = 
  m_is_nonlinacc_computed = 
  m_is_linacc_computed =
  m_is_jerk_computed = 
  m_is_nonlinjerk_computed = 
  m_is_linjerk_computed =
  m_is_wrench_computed =
  m_is_regressor_computed =
    false;

  m_gravity = gravity;
  boost::shared_ptr<itia::dynamics::Link> base_link = root_link->findChild(base_link_name);
  if (!base_link)
  {
    ROS_ERROR("Base link not found");
    return;
  }
  boost::shared_ptr<itia::dynamics::Link> ee_link = base_link->findChild(ee_link_name);
  if (!ee_link)
  {
    ROS_ERROR("Tool link not found");
    return;
  }
  
  boost::shared_ptr<itia::dynamics::Link> act_link(ee_link);
  while (1)
  {
    m_links.insert(m_links.begin(), act_link);
    if (act_link->getName().compare(base_link_name))
    {
      m_joints.insert(m_joints.begin(), act_link->getParentJoint());
      act_link = act_link->getParentJoint()->getParentLink();
    }
    else
      break;
  }
  
  for (unsigned int idx = 0;idx<m_links.size();idx++)
    m_links_name.push_back( m_links.at(idx)->getName() );
   
  for (unsigned int idx = 0;idx<m_joints.size();idx++)
    m_joints_name.insert( std::pair<std::string, unsigned int>(m_joints.at(idx)->getName(), idx) );
  
  m_active_joints_number = m_joints_number = m_joints.size();
  m_links_number  = m_links.size();
  m_joint_inertia_extended.resize(m_joints_number, m_joints_number);
  m_joint_inertia_extended.setZero();
  
  
  m_input_to_chain_joint.resize(m_joints_number, m_joints_number);
  m_input_to_chain_joint.setIdentity();
  m_chain_to_input_joint = m_input_to_chain_joint;
  
  m_T_bt.setIdentity();
  m_last_q.resize(m_joints_number);
  m_last_q.setZero();
  
  m_last_Dq.resize(m_joints_number);
  m_last_Dq.setZero();
  m_last_DDq.resize(m_joints_number);
  m_last_DDq.setZero();
  m_last_DDDq.resize(m_joints_number);
  m_last_DDDq.setZero();
  
  m_active_joints.resize(m_joints_number);
  
  m_wrenches_regressor.resize(m_links_number);
  
  
  Eigen::Vector6d zero_vector;
  zero_vector.setZero();
  for (unsigned int idx = 0;idx<m_links_number;idx++)
  {
    m_screws_of_c_in_b.push_back(zero_vector);
    m_twists.push_back(zero_vector);
    m_Dtwists.push_back(zero_vector);
    m_Dtwists_linear_part.push_back(zero_vector);
    m_Dtwists_nonlinear_part.push_back(zero_vector);
    m_DDtwists.push_back(zero_vector);
    m_DDtwists_linear_part.push_back(zero_vector);
    m_DDtwists_nonlinear_part.push_back(zero_vector);
    m_wrenches.push_back(zero_vector);
    m_inertial_wrenches.push_back(zero_vector);
    m_gravity_wrenches.push_back(zero_vector);
    m_wrenches_regressor.at(idx).setZero();
  }
  
  
  m_sorted_q.resize(m_joints_number);
  m_sorted_Dq.resize(m_joints_number);
  m_sorted_DDq.resize(m_joints_number);
  m_sorted_DDDq.resize(m_joints_number);
  m_joint_torques.resize(m_joints_number);
  m_active_joint_torques.resize(m_joints_number);
  m_sorted_q.setZero();
  m_sorted_Dq.setZero();
  m_sorted_DDq.setZero();
  m_sorted_DDDq.setZero();
  m_joint_torques.setZero();
  
  m_regressor_extended.resize(m_joints_number, m_joints_number*10);
  m_regressor_extended.setZero();
  
  m_T_bl.resize(m_links_number);
  m_T_bl.at(0).setIdentity();
  computeFrames();
  
};

void Chain::setInputJointsName(const std::vector< std::string > joints_name)
{
  m_input_to_chain_joint.resize(m_joints_number, joints_name.size());
  
  m_active_joints_number = joints_name.size();
  m_input_to_chain_joint.setZero();
  m_active_joints.resize(joints_name.size());
//   m_regressor_extended.resize(m_active_joints_number, m_joints_number*10);
  
  m_regressor_extended.setZero();
  m_joint_inertia.resize(m_active_joints_number, m_active_joints_number);
  m_joint_inertia.setZero();
  
  for (unsigned int idx = 0;idx<joints_name.size();idx++)
    if (m_joints_name.find(joints_name.at(idx)) != m_joints_name.end())
    {
      m_input_to_chain_joint(m_joints_name.find(joints_name.at(idx))->second, idx) = 1;
      m_active_joints.at(idx) = m_joints_name.find(joints_name.at(idx))->second;
    }
    else
      ROS_WARN("Joint named '%s' not found", joints_name.at(idx).c_str());
      
  m_active_joint_torques.resize(m_active_joints_number);
  m_chain_to_input_joint = m_input_to_chain_joint.transpose();
  
  m_last_q.resize(joints_name.size());
  m_last_q.setZero();
  
  m_sorted_q.setZero();
  m_sorted_Dq.setZero();
  m_sorted_DDq.setZero();
  m_sorted_DDDq.setZero();
  m_jacobian.resize(6, joints_name.size());
  m_jacobian.setZero();

  computeFrames();
  
  m_joint_inertia.resize(m_active_joints_number, m_active_joints_number);
  m_joint_inertia.setZero();
  
  Eigen::Matrix4d eye4;
  eye4.setIdentity();
  

}

void Chain::computeFrames()
{
  m_sorted_q = m_input_to_chain_joint*m_last_q;
  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    m_T_bl.at(nl).matrix() = m_T_bl.at(nl-1).matrix() * m_joints.at(nj)->getTransformation(m_sorted_q(nj)).matrix();
  }
  m_T_bt = m_T_bl.at(m_links_number-1);

}

void Chain::computeScrews()
{
  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    spatialRotation(m_joints.at(nj)->getScrew_of_child_in_parent(), m_T_bl.at(nj).linear().matrix(), &(m_screws_of_c_in_b.at(nl)));
  }
  m_is_screws_computed = true;
}

Eigen::Affine3d Chain::getTransformation(const Eigen::VectorXd& q)
{  
  if ( (q == m_last_q) || (m_joints_number == 0) )  
    return m_T_bt;
 
  m_last_q = q;
  m_is_screws_computed = 
  m_is_jac_computed = 
  m_is_vel_computed = 
  m_is_acc_computed = 
  m_is_nonlinacc_computed = 
  m_is_linacc_computed =
  m_is_jerk_computed = 
  m_is_nonlinjerk_computed = 
  m_is_linjerk_computed =
  m_is_wrench_computed = 
  m_is_regressor_computed =
    false;
  
  computeFrames();
  
  return m_T_bt;
}

Eigen::Matrix6Xd Chain::getJacobian(const Eigen::VectorXd& q)
{
  getTransformation(q);
  if (m_joints_number == 0)
    return m_jacobian;
  
  if (m_is_jac_computed)
    return m_jacobian;
    
  if (!m_is_screws_computed)
    computeScrews();
    
  for (unsigned int idx = 0;idx<m_active_joints.size();idx++)
  {
    unsigned int nj = m_active_joints.at(idx);
    unsigned int nl = nj+1;
    if (!m_joints.at(nj)->isFixed())
      m_jacobian.col(idx) = spatialTranslation(m_screws_of_c_in_b.at(nl), m_T_bt.translation()-m_T_bl.at(nl).translation());
  }
  
  m_is_jac_computed = true;
  return m_jacobian;
}


std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> Chain::getTwist(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq)
{
  getTransformation(q);
  m_sorted_Dq = m_input_to_chain_joint*Dq;
  if ( (m_sorted_Dq - m_last_Dq).norm()>1e-12)
  {
    m_is_vel_computed = false;
    m_is_acc_computed = false;
    m_is_jerk_computed = false;
    m_is_wrench_computed = false;
    m_is_regressor_computed = false;
    m_is_linacc_computed = false;
    m_is_nonlinacc_computed = false;
    m_is_nonlinjerk_computed = false;
    m_is_linjerk_computed = false;
  }
  m_last_Dq = m_sorted_Dq;
  if (m_is_vel_computed)
    return m_twists;
  
  if (!m_is_screws_computed)
    computeScrews();

  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    m_twists.at(nl) = spatialTranslation(m_twists.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation()) + 
                      m_screws_of_c_in_b.at(nl) * m_sorted_Dq(nj);
  }
  m_is_vel_computed = true;

  return m_twists;
}

std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > Chain::getDTwistLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& DDq)
{
  getTransformation(q);
  m_sorted_DDq = m_input_to_chain_joint*DDq;
  if ( (m_sorted_DDq - m_last_DDq).norm()>1e-12)
  {
    m_is_acc_computed = false;
    m_is_jerk_computed = false;
    m_is_wrench_computed = false;
    m_is_regressor_computed = false;
    m_is_linacc_computed = false;
    m_is_nonlinacc_computed = false;
    m_is_nonlinjerk_computed = false;
    m_is_linjerk_computed = false;
  }
  m_last_DDq = m_sorted_DDq;
  
  if (m_is_linacc_computed)
    return m_Dtwists_linear_part;
  
  if (!m_is_screws_computed)
    computeScrews();

  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    m_Dtwists_linear_part.at(nl) = spatialTranslation(m_Dtwists_linear_part.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation()) + 
                      m_screws_of_c_in_b.at(nl) * m_sorted_DDq(nj);
  }
  m_is_linacc_computed= true;

  return m_Dtwists_linear_part;
};

std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> Chain::getDTwistNonLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq)
{
  getTransformation(q);
  if (!m_is_vel_computed)
    getTwist(q, Dq);
  if (m_is_nonlinacc_computed)
    return m_Dtwists_nonlinear_part;
  
  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    m_Dtwists_nonlinear_part.at(nl) = spatialTranslation(m_Dtwists_nonlinear_part.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation()) + 
                                      spatialCrossProduct(m_twists.at(nl), m_screws_of_c_in_b.at(nl)) * m_sorted_Dq(nj);
  }
  m_is_nonlinacc_computed = true;

  return m_Dtwists_nonlinear_part;
}

std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > Chain::getDTwist(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq)
{
  
  getTransformation(q);
  
  m_sorted_DDq = m_input_to_chain_joint*DDq;
  
  if ( (m_sorted_DDq - m_last_DDq).norm()>1e-12)
  {
    m_is_acc_computed = false;
    m_is_jerk_computed = false;
    m_is_wrench_computed = false;
    m_is_regressor_computed = false;
    m_is_linacc_computed = false;
    m_is_nonlinacc_computed = false;
    m_is_nonlinjerk_computed = false;
    m_is_linjerk_computed = false;
  }
  m_last_DDq = m_sorted_DDq;
  
  if (m_is_acc_computed)
    return m_Dtwists;
  if (m_is_linacc_computed && m_is_nonlinacc_computed)
  {
    for (unsigned int nl = 1;nl<m_links_number;nl++)
      m_Dtwists.at(nl) =  m_Dtwists_nonlinear_part.at(nl) + m_Dtwists_linear_part.at(nl);
  }
  else 
  {
    m_sorted_DDq = m_input_to_chain_joint*DDq;
    if (!m_is_vel_computed)
      getTwist(q, Dq);
    for (unsigned int nl = 1;nl<m_links_number;nl++)
    {
      unsigned int nj = nl-1;
      m_Dtwists.at(nl) = spatialTranslation(m_Dtwists.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation()) + 
                         spatialCrossProduct(m_twists.at(nl), m_screws_of_c_in_b.at(nl)) * m_sorted_Dq(nj) +  m_screws_of_c_in_b.at(nl) * m_sorted_DDq(nj);
    }
  }
  
  m_is_acc_computed = true;
  
  return m_Dtwists;
}

std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > Chain::getDDTwistLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& DDDq)
{
  getTransformation(q);
  
  m_sorted_DDDq = m_input_to_chain_joint*DDDq;
  if ( (m_sorted_DDDq - m_last_DDDq).norm()>1e-12)
  {
    m_is_jerk_computed = false;
    m_is_nonlinjerk_computed = false;
    m_is_linjerk_computed = false;
  }
  m_last_DDDq = m_sorted_DDDq;
  
  if (m_is_linjerk_computed)
    return m_Dtwists_linear_part;
  
  if (!m_is_screws_computed)
    computeScrews();

  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    m_DDtwists_linear_part.at(nl) = spatialTranslation(m_DDtwists_linear_part.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation()) + 
                      m_screws_of_c_in_b.at(nl) * m_sorted_DDDq(nj);
  }
  m_is_linjerk_computed = true;

  return m_DDtwists_linear_part;
};

std::vector<Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d>> Chain::getDDTwistNonLinearPart(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq)
{
  getTransformation(q);
  if (!m_is_acc_computed)
    getDTwist(q, Dq, DDq);
  
  if (m_is_nonlinjerk_computed)
    return m_DDtwists_nonlinear_part;
  
  /*
    * a(i)= a(i-1)+S(i)*DDq(i)  + v(i) X S(i)*Dq(i) -> DS(i) = v(i) X S(i)
    * j(i)= j(i-1)+S(i)*DDDq(i) + a(i) X S(i)*Dq(i)  + v(i) X S(i)*DDq(i) + v(i) X DS(i)*Dq(i)
    * j(i)= j(i-1)+S(i)*DDDq(i) + v(i) X S(i)*DDq(i) + a(i) X S(i)*Dq(i) + v(i) X v(i) X S(i) * Dq(i)
    */
  Eigen::Vector6d v_cross_s;
  for (unsigned int nl = 1;nl<m_links_number;nl++)
  {
    unsigned int nj = nl-1;
    spatialCrossProduct(m_twists.at(nl), m_screws_of_c_in_b.at(nl), &v_cross_s);
    m_DDtwists_nonlinear_part.at(nl) =  spatialTranslation(m_DDtwists_nonlinear_part.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation()) + 
                                        v_cross_s * m_sorted_DDq(nj) +   /* v(i) X S(i)*DDq(i) */
                                        ( spatialCrossProduct( m_Dtwists.at(nl), m_screws_of_c_in_b.at(nl) ) + /* a(i) X S(i)*Dq(i) */
                                          spatialCrossProduct( m_twists.at(nl), v_cross_s ) ) * m_sorted_Dq(nj); /* v(i) X v(i) X S(i) * Dq(i) */

  }
  m_is_nonlinjerk_computed = true;

  return m_DDtwists_nonlinear_part;
}

std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > Chain::getDDTwist(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq, const Eigen::VectorXd& DDDq)
{
  getTransformation(q);
  m_sorted_DDDq = m_input_to_chain_joint*DDDq;
  if ( (m_sorted_DDDq - m_last_DDDq).norm()>1e-12)
  {
    m_is_jerk_computed = false;
    m_is_nonlinjerk_computed = false;
    m_is_linjerk_computed = false;
  }
  m_last_DDDq = m_sorted_DDDq;
  
  if (m_is_jerk_computed)
    return m_DDtwists;
  if (m_is_linjerk_computed && m_is_nonlinjerk_computed)
  {
    for (unsigned int nl = 1;nl<m_links_number;nl++)
      m_DDtwists.at(nl) =  m_DDtwists_nonlinear_part.at(nl) + m_DDtwists_linear_part.at(nl);
  }
  else 
  {
    m_sorted_DDDq = m_input_to_chain_joint*DDDq;
    if (!m_is_acc_computed)
      getDTwist(q, Dq, DDq);
    Eigen::Vector6d v_cross_s;
    for (unsigned int nl = 1;nl<m_links_number;nl++)
    {
      unsigned int nj = nl-1;
      spatialCrossProduct(m_twists.at(nl), m_screws_of_c_in_b.at(nl), &v_cross_s);
      m_DDtwists.at(nl) = spatialTranslation( m_DDtwists.at(nl-1), m_T_bl.at(nl).translation()-m_T_bl.at(nl-1).translation() ) + 
                          m_screws_of_c_in_b.at(nl) * m_sorted_DDDq(nj) +
                          v_cross_s * m_sorted_DDq(nj) +   /* v(i) X S(i)*DDq(i) */
                          ( spatialCrossProduct( m_Dtwists.at(nl), m_screws_of_c_in_b.at(nl) ) + /* a(i) X S(i)*Dq(i) */
                            spatialCrossProduct( m_twists.at(nl), v_cross_s ) ) * m_sorted_Dq(nj); /* v(i) X v(i) X S(i) * Dq(i) */
    }
  }
  m_is_jerk_computed = true;
  return m_DDtwists;
}




std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > Chain::getWrench(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq, std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > ext_wrenches_in_link_frame)
{
  getDTwist(q, Dq, DDq);
  if (m_is_wrench_computed)
    return m_wrenches;
  
  for (int nl = (m_links_number-1);nl >= 0;nl--)
  { 
    if (nl == 0)
    {
      m_inertial_wrenches.at(nl).setZero();
      m_gravity_wrenches.at(nl).setZero();
    } 
    else
    {
      m_inertial_wrenches.at(nl) = spatialRotation(
                                   m_links.at(nl)->getSpatialInertia() * 
                                  spatialRotation(m_Dtwists.at(nl), m_T_bl.at(nl).linear().transpose())
                                  + 
                                  spatialDualCrossProduct(
                                    spatialRotation(m_twists.at(nl), m_T_bl.at(nl).linear().transpose()), 
                                    m_links.at(nl)->getSpatialInertia() * 
                                    spatialRotation(m_twists.at(nl), m_T_bl.at(nl).linear().transpose())
                                    ), m_T_bl.at(nl).linear());
      m_gravity_wrenches.at(nl).block(0, 0, 3, 1) = -m_links.at(nl)->getMass() *m_gravity;
      m_gravity_wrenches.at(nl).block(3, 0, 3, 1) = -(m_T_bl.at(nl).linear() * m_links.at(nl)->getCog() ).cross(m_links.at(nl)->getMass() * m_gravity);
      
    
    }
    
    if (nl < (m_links_number-1))
      m_wrenches.at(nl) = spatialTranformation(ext_wrenches_in_link_frame.at(nl), m_T_bl.at(nl)) + m_inertial_wrenches.at(nl) + m_gravity_wrenches.at(nl) + spatialDualTranslation( m_wrenches.at(nl+1) , m_T_bl.at(nl).translation()-m_T_bl.at(nl+1).translation() );
    else
      m_wrenches.at(nl) = spatialTranformation(ext_wrenches_in_link_frame.at(nl), m_T_bl.at(nl)) + m_inertial_wrenches.at(nl) + m_gravity_wrenches.at(nl);
      
  }
  
  m_is_wrench_computed = true;
  return m_wrenches;
}

Eigen::VectorXd Chain::getJointTorque(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq, std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > ext_wrenches_in_link_frame)
{
  getWrench(q, Dq, DDq, ext_wrenches_in_link_frame);
  for (unsigned int nj = 0; nj < m_joints_number; nj++)
  {
    unsigned int nl = nj+1;
    m_joint_torques(nj) = m_wrenches.at(nl).dot(m_screws_of_c_in_b.at(nl));
  }
  m_active_joint_torques = m_chain_to_input_joint * m_joint_torques;
  return m_active_joint_torques;
}


Eigen::VectorXd Chain::getJointTorque(const Eigen::VectorXd& q, const Eigen::VectorXd& Dq, const Eigen::VectorXd& DDq)
{
  std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > ext_wrenches_in_link_frame(m_links_number);
  for (unsigned int iL = 0;iL < m_links_number; iL++)
    ext_wrenches_in_link_frame.at(iL).setZero();
  return getJointTorque(q, Dq, DDq, ext_wrenches_in_link_frame);
}

Eigen::MatrixXd Chain::getRegressor(const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& Dq, 
                                    const Eigen::VectorXd& DDq)
{
  if (q.rows() != Dq.rows() )
  {
    ROS_ERROR("Input data dimensions mismatch");
    throw std::invalid_argument("Input data dimensions mismatch");
  } 
    
  if (Dq.rows() != DDq.rows() )
  {
    ROS_ERROR("Input data dimensions mismatch");
    throw std::invalid_argument("Input data dimensions mismatch");
  }  
  
  getDTwist(q, Dq, DDq);
  
  if (m_is_regressor_computed)
  {
    static bool verbose_ = true;
    if (verbose_)
      ROS_WARN("Regressor input element equals to the previous call. Returned the same dynamics Regressor");
    
    return m_chain_to_input_joint * m_regressor_extended;
  }
    
  for (int nl = (m_links_number-1);nl > 0;nl--)
  {
    // m, mcx, mcy, mcz, Ixx, Ixy, Ixz, Iyy, Iyz, Izz
    for (unsigned int iPar = 0;iPar<10;iPar++)
    {
      m_wrenches_regressor.at(nl).col(iPar) = spatialRotation(m_links.at(nl)->getSpatialInertiaTerms().at(iPar) * 
                                                              spatialRotation(m_Dtwists.at(nl), m_T_bl.at(nl).linear().transpose()) + 
                                                              spatialDualCrossProduct(
                                                              spatialRotation(m_twists.at(nl), m_T_bl.at(nl).linear().transpose()), 
                                                              m_links.at(nl)->getSpatialInertiaTerms().at(iPar) * 
                                                              spatialRotation(m_twists.at(nl), m_T_bl.at(nl).linear().transpose())
                                                              ), m_T_bl.at(nl).linear());
    }
    

    m_wrenches_regressor.at(nl).col(0).block(0, 0, 3, 1) -=  m_gravity;
    m_wrenches_regressor.at(nl).col(1).block(3, 0, 3, 1) -=  (m_T_bl.at(nl).linear() * Eigen::Vector3d::UnitX() ).cross(m_gravity);
    m_wrenches_regressor.at(nl).col(2).block(3, 0, 3, 1) -=  (m_T_bl.at(nl).linear() * Eigen::Vector3d::UnitY() ).cross(m_gravity);
    m_wrenches_regressor.at(nl).col(3).block(3, 0, 3, 1) -=  (m_T_bl.at(nl).linear() * Eigen::Vector3d::UnitZ() ).cross(m_gravity);
    
    m_regressor_extended.block(nl-1, (nl-1)*10, 1, 10) = m_screws_of_c_in_b.at(nl).transpose() * m_wrenches_regressor.at(nl);
   
    for (unsigned int nl_following = nl+1; nl_following<m_links_number; nl_following++)
    {
      for (unsigned int iPar = 0;iPar<10;iPar++)
        m_regressor_extended(nl-1, (nl_following-1)*10+iPar) = m_screws_of_c_in_b.at(nl).transpose() * spatialDualTranslation( m_wrenches_regressor.at(nl_following).col(iPar) , m_T_bl.at(nl).translation() - m_T_bl.at(nl_following).translation() );
    }      
  }
  
  m_is_regressor_computed = true;
    
  return m_chain_to_input_joint * m_regressor_extended;
}

Eigen::MatrixXd Chain::getJointInertia(const Eigen::VectorXd& q)
{
  getTransformation(q);
  computeScrews();
  m_joint_inertia_extended.setZero();
  for (unsigned int nj = 0; nj<m_joints_number; nj++)
  {
    Eigen::Matrix6Xd jacobian_nl(6, m_joints_number);
    jacobian_nl.setZero();
    for (unsigned int ij = 0;ij <= nj/*m_active_joints.size()*/;ij++)
    {
      unsigned int il = ij+1;
      if (!m_joints.at(ij)->isFixed())
      {
        jacobian_nl.col(ij) = spatialTranslation(m_screws_of_c_in_b.at(il), m_T_bl.at(nj+1).translation()-m_T_bl.at(il).translation());
        jacobian_nl.col(ij) = spatialRotation(jacobian_nl.col(ij), m_T_bl.at(nj+1).linear().transpose());
      }
    }
    m_joint_inertia_extended += jacobian_nl.transpose() * m_links.at(nj+1)->getSpatialInertia() * jacobian_nl;
  }
  m_joint_inertia = m_chain_to_input_joint * m_joint_inertia_extended * m_input_to_chain_joint;
  return m_joint_inertia;
};


// Eigen::VectorXd Chain::getJoint(const Eigen::Affine3d& T_bt_goal, const Eigen::VectorXd& q0)
// {
//   Eigen::VectorXd q = q0;
//   getJacobian(q);
//   Eigen::Vector6d m_cartesian_error;
//   Eigen::MatrixXd m_error_jacobian;
//   itia::mutils::getFrameDistanceJac(m_T_bt, T_bt_goal, &m_cartesian_error, &m_error_jacobian);
//   
//   
//   /* NOTE
//    * Controller update 
//    * dq(k) = (J+(K*Je*J))\(K*e)
//    * q(k+1)= q(k)+dq(k)
//    */
//   bool  flag = true;
//   int ntrial = 5;
//   int itrial = 0;
//   double m_gain = 1;
//   
//   
//   Eigen::Matrix<double, -1, 6> den_control_matrix;
//   Eigen::VectorXd m_dq(m_active_joints.size());
//   
//   while (flag && itrial<ntrial)
//   {
//     itia::mutils::getFrameDistanceJac(m_T_bt, T_bt_goal, &m_cartesian_error, &m_error_jacobian);
//     den_control_matrix = m_jacobian + m_jacobian*m_error_jacobian*m_gain;
//     Eigen::FullPivLU<Eigen::MatrixXd> den_control_matrix_decomp(den_control_matrix);
//     den_control_matrix_decomp.setThreshold(1e-2);
//     Eigen::Matrix<double, 6, 1> num_control_vector = m_gain*m_cartesian_error;
//     m_dq = den_control_matrix_decomp.solve(num_control_vector);
//     if (den_control_matrix_decomp.rank()<6)
//     {
//       Eigen::MatrixXd null=den_control_matrix_decomp.kernel();
//       ROS_WARN_THROTTLE(0.1,"Singolarity point!");
//       m_dq = m_dq - (null.transpose()-m_dq) *null;
//     }
//     
//     q.noalias() += m_dq;
//     getJacobian(q);
//     itrial++;
//   }
// 
// }


}
}
