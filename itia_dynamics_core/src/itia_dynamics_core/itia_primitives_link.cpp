
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

#include <itia_dynamics_core/itia_primitives.h>

namespace itia 
{
namespace dynamics 
{


  
void Link::fromUrdf(const boost::shared_ptr<urdf::Link>& urdf_link, const boost::shared_ptr<itia::dynamics::Joint>& parent_joint)
{
  m_parent_joint = parent_joint;
  m_name = urdf_link->name;
  for (unsigned int idx = 0; idx<urdf_link->child_joints.size();idx++)
  {
    m_child_joints.push_back( boost::shared_ptr<itia::dynamics::Joint>(new itia::dynamics::Joint()));
    m_child_joints.back()->fromUrdf(urdf_link->child_joints.at(idx),pointer(), urdf_link->child_links.at(idx));
    m_child_links.push_back(m_child_joints.at(idx)->getChildLink());
  }
  
  m_mass = 0;
  Eigen::Matrix3d inertia;
  m_Inertia_cc_single_term.resize(10);
  if (urdf_link->inertial != NULL)
  {
    m_mass = urdf_link->inertial->mass;
    
    inertia(0, 0) = urdf_link->inertial->ixx;
    inertia(0, 1) = urdf_link->inertial->ixy;
    inertia(0, 2) = urdf_link->inertial->ixz;
    inertia(1, 0) = urdf_link->inertial->ixy;
    inertia(1, 1) = urdf_link->inertial->iyy;
    inertia(1, 2) = urdf_link->inertial->iyz;
    inertia(2, 0) = urdf_link->inertial->ixz;
    inertia(2, 1) = urdf_link->inertial->iyz;
    inertia(2, 2) = urdf_link->inertial->izz;
    
    m_cog_in_c(0) = urdf_link->inertial->origin.position.x;
    m_cog_in_c(1) = urdf_link->inertial->origin.position.y;
    m_cog_in_c(2) = urdf_link->inertial->origin.position.z;
    
    Eigen::Quaterniond q_p_cog;
    
    q_p_cog.x() = urdf_link->inertial->origin.rotation.x;
    q_p_cog.y() = urdf_link->inertial->origin.rotation.y;
    q_p_cog.z() = urdf_link->inertial->origin.rotation.z;
    q_p_cog.w() = urdf_link->inertial->origin.rotation.w;
    
    // inertia (p) = R_p_cog * inertia_cog * R_cog_p
    inertia = q_p_cog.toRotationMatrix() * inertia * q_p_cog.toRotationMatrix().transpose();
    computeSpatialInertiaMatrix(inertia, m_cog_in_c, m_mass, m_Inertia_cc);
    
    
  }
  else
  {
    inertia.setZero();
    m_cog_in_c.setZero();
    
    computeSpatialInertiaMatrix(inertia, m_cog_in_c, m_mass, m_Inertia_cc);
    for (int idx = 0;idx<10;idx++)
      m_Inertia_cc_single_term.at(idx).setZero();
  }
  
  
  
  Eigen::Vector3d cog;
  // m, mcx, mcy, mcz, Ixx, Ixy, Ixz, Iyy, Iyz, Izz
  inertia.setZero();
  cog.setZero();
  double mass = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(0));
  mass = 0;
  
  // mcx
  cog(0) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(1));
  cog(0) = 0;
  
  // mcy
  cog(1) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(2));
  cog(1) = 0;
  
  // mcz
  cog(2) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(3));
  cog(2) = 0;
  
  // Ixx
  inertia(0, 0) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(4));
  inertia(0, 0) = 0;
  
  // Ixy
  inertia(0, 1) = 1;
  inertia(1, 0) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(5));
  inertia(0, 1) = 0;
  inertia(1, 0) = 0;
  
  // Ixz
  inertia(0, 2) = 1;
  inertia(2, 0) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(6));
  inertia(0, 2) = 0;
  inertia(2, 0) = 0;
  
  // Iyy
  inertia(1, 1) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(7));
  inertia(1, 1) = 0;
  
  // Iyz
  inertia(1, 2) = 1;
  inertia(2, 1) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(8));
  inertia(1, 2) = 0;
  inertia(2, 1) = 0;
  
  // Izz
  inertia(2, 2) = 1;
  computeSpatialInertiaMatrix(inertia, cog, mass, m_Inertia_cc_single_term.at(9));
  inertia(2, 2) = 0;
  
}
boost::shared_ptr<itia::dynamics::Link> Link::pointer()
{
  return shared_from_this();
};

boost::shared_ptr<itia::dynamics::Link> Link::findChild(const std::string& name)
{
  boost::shared_ptr<itia::dynamics::Link> ptr;
  if (!m_name.compare(name))
    return pointer();
  if (m_child_links.size() == 0)
    return ptr;
  for (unsigned int idx = 0;idx<m_child_links.size();idx++)
  {
    if (!m_child_links.at(idx)->getName().compare(name))
      return m_child_links.at(idx);
    ptr = m_child_links.at(idx)->findChild(name);
    if (ptr)
      return ptr;
  }
  return ptr;
};

boost::shared_ptr< Joint > Link::findChildJoint(const std::string& name)
{
  boost::shared_ptr<itia::dynamics::Joint> ptr;
  if (m_child_joints.size() == 0)
    return ptr;
  for (unsigned int idx = 0;idx<m_child_joints.size();idx++)
  {
    if (!m_child_joints.at(idx)->getName().compare(name))
      return m_child_joints.at(idx);
    ptr = m_child_links.at(idx)->findChildJoint(name);
    if (ptr)
      return ptr;
  }
  return ptr;
}




}
}