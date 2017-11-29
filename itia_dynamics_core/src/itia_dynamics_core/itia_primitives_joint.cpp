#include <itia_dynamics_core/itia_primitives.h>

namespace itia 
{
namespace dynamics 
{



Joint::Joint()
{
  m_last_q = 0;
  
  m_last_T_pc.setIdentity();
  m_last_R_jc.setIdentity();
  
  m_identity.setIdentity();
  m_screw_of_c_in_p.setZero();
  
  computedTpc();
}

/*
 * Rpc=Rpj*Rjc
 * DRpc=Rpj*DRjc
 * S=Rpj*DRjc*Rjc'*Rpj'
 * S=Rpj*Sjc*Rpj'
 * 
 * W_pp = Rpj*W_jj*T_jp = [S_pp', -S_pp*p_pc+R_pj*vj;0 0]
 * 
*/
void Joint::computeJacobian()
{
  if (m_type == REVOLUTE)
  {
    m_screw_of_c_in_p << Eigen::MatrixXd::Constant(3, 1, 0), m_axis_in_p;
  }
  else if (m_type == PRISMATIC)
  {
    m_screw_of_c_in_p << m_axis_in_p, Eigen::MatrixXd::Constant(3, 1, 0);
  }
}


void Joint::computedTpc()
{
  if (m_type == REVOLUTE)
  {
    m_last_R_jc = m_identity+sin(m_last_q) *m_skew_axis_in_j+(1-cos(m_last_q)) *m_square_skew_axis_in_j;
    m_last_T_pc.linear().matrix()= m_R_pj*m_last_R_jc;     //m_T_pj+sin(m_last_q) *m_skew_axis_in_p+(1-cos(m_last_q)) *m_square_skew_axis_in_p;
  }
  else if (m_type == PRISMATIC)
    m_last_T_pc.translation() = m_T_pj.translation()+m_axis_in_p*m_last_q;
}


void Joint::fromUrdf(const boost::shared_ptr< urdf::Joint >& urdf_joint, const boost::shared_ptr< Link >& parent_link, const boost::shared_ptr< urdf::Link >& child_link)
{
  m_parent_link = parent_link;  
//   urdf_joint->name;
  m_T_pj = urdfPoseToAffine(urdf_joint->parent_to_joint_origin_transform);
  m_axis_in_j = urdfVectorToEigen(urdf_joint->axis);
  m_name = urdf_joint->name;
  if (m_axis_in_j.norm()>0)
    m_axis_in_j /= m_axis_in_j.norm();
  
  
  m_skew_axis_in_j = skew(m_axis_in_j);
  m_square_skew_axis_in_j = m_skew_axis_in_j*m_skew_axis_in_j;
  m_identity.setIdentity();
  
  
  // Transform in p frame
  m_R_pj = m_T_pj.linear();
  m_axis_in_p = m_R_pj*m_axis_in_j;
  m_skew_axis_in_p = m_T_pj.linear()*m_skew_axis_in_j;
  m_square_skew_axis_in_p =m_T_pj.linear() *m_square_skew_axis_in_j;
  m_last_T_pc = m_T_pj;
  
  
  if (urdf_joint->type == urdf::Joint::REVOLUTE)
    m_type = itia::dynamics::Joint::REVOLUTE;
  else if (urdf_joint->type == urdf::Joint::CONTINUOUS)
    m_type = itia::dynamics::Joint::REVOLUTE;
  else if (urdf_joint->type == urdf::Joint::PRISMATIC)
  {
    m_type = itia::dynamics::Joint::PRISMATIC;
  }
  else
    m_type = itia::dynamics::Joint::FIXED;
  
  m_child_link.reset(new itia::dynamics::Link() );
  m_child_link->fromUrdf(child_link, pointer());
  computedTpc();
  computeJacobian();
}

boost::shared_ptr<itia::dynamics::Joint> Joint::pointer()
{
  return shared_from_this();
}

const Eigen::Affine3d& Joint::getTransformation(const double& q)
{
  if (q != m_last_q)
  {
    m_last_q = q;
    computedTpc();
  }
  return m_last_T_pc;
}


const Eigen::Vector6d& Joint::getScrew_of_child_in_parent()
{
  return m_screw_of_c_in_p;
}

}
}