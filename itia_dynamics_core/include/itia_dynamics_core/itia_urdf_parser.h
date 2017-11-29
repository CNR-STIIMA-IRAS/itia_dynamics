#ifndef __ITIA_KIN_DERIVATIVES__
#define __ITIA_KIN_DERIVATIVES__


#include "itia_dynamics_core/itia_primitives.h"

#include <urdf/model.h>
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>

#include <kdl_parser/kdl_parser.hpp>


namespace itia
{
namespace dynamics
{
  
  inline Eigen::Affine3d urdfPoseToAffine(urdf::Pose pose)
  {
    Eigen::Affine3d affine;
    affine = Eigen::Quaterniond(pose.rotation.w, pose.rotation.x, pose.rotation.y, pose.rotation.z);
    affine.translation() << pose.position.x, pose.position.y, pose.position.z;
    return affine;  
  }
  
  inline Eigen::Vector3d urdfVectorToEigen(urdf::Vector3 vector)
  {
    Eigen::Vector3d eigen_vector;
    eigen_vector  << vector.x, vector.y, vector.z;
    return eigen_vector;
  }
  
//   class Link;
//   
//   class Joint: public boost::enable_shared_from_this<itia::dynamics::Joint>
//   {
//   protected:
//     enum
//     {
//       REVOLUTE, PRISMATIC, FIXED
//     } m_type;
//     
//     Eigen::Affine3d m_T_pj; // parent <- joint
//     Eigen::Vector3d m_axis;
//     
//     std::string m_name;
//     boost::shared_ptr<itia::dynamics::Link> m_parent_link;
//     boost::shared_ptr<itia::dynamics::Link> m_child_link;
// 
//   public:
//     Joint(const boost::shared_ptr<urdf::Joint>& urdf_joint, const boost::shared_ptr<itia::dynamics::Link>& parent_link, const boost::shared_ptr<urdf::Link>& child_link);
//     boost::shared_ptr<itia::dynamics::Joint> pointer()
//     {
//       return shared_from_this();
//     }
//     std::string getName(){return m_name;};
//     boost::shared_ptr <itia::dynamics::Link> getChildLink() {return m_child_link;}; 
//     boost::shared_ptr <itia::dynamics::Link> getParentLink() {return m_parent_link;}; 
//     
//   };
//   
/*
  
  class Link: public boost::enable_shared_from_this<itia::dynamics::Link>
  {
  protected:
    boost::shared_ptr<itia::dynamics::Joint> m_parent_joint;
    std::vector<boost::shared_ptr<itia::dynamics::Joint>> m_child_joints;                     
    std::vector<boost::shared_ptr<itia::dynamics::Link>> m_child_links;                     
    std::string m_name;
  public:
    Link(const boost::shared_ptr<urdf::Link>& urdf_link, const boost::shared_ptr<itia::dynamics::Joint>& parent_joint = 0)
    :
    m_parent_joint(parent_joint)
    {
      m_name = urdf_link->name;
      ROS_DEBUG_STREAM("Adding link named: " <<  urdf_link->name);
      for (int idx = 0; idx<urdf_link->child_joints.size();idx++)
      {
        pointer();
        ROS_DEBUG_STREAM("Adding joint named: " <<  urdf_link->child_joints.at(idx)->name);
        m_child_joints.push_back( boost::shared_ptr<itia::dynamics::Joint>(
          new itia::dynamics::Joint( urdf_link->child_joints.at(idx),
          pointer(), 
          urdf_link->child_links.at(idx)
        )));
        ROS_DEBUG_STREAM("DONE");
        m_child_links.push_back(m_child_joints.at(idx)->getChildLink());
      }
    }
    boost::shared_ptr<itia::dynamics::Link> pointer()
    {
      return shared_from_this();
    };
    std::string getName(){return m_name;};
    
    boost::shared_ptr<itia::dynamics::Joint> getParentJoint() {return m_parent_joint;};
    boost::shared_ptr<itia::dynamics::Link> findChild(const std::string& name)
    {
      boost::shared_ptr<itia::dynamics::Link> ptr;
      if (m_child_links.size() == 0)
        return ptr;
      for (int idx = 0;idx<m_child_links.size();idx++)
      {
        if (!m_child_links.at(idx)->getName().compare(name))
          return m_child_links.at(idx);
        ptr = m_child_links.at(idx)->findChild(name);
        if (!ptr)
          return ptr;
      }
      return ptr;
    };
  };
    
  Joint::Joint(const boost::shared_ptr< urdf::Joint >& urdf_joint, const boost::shared_ptr< Link >& parent_link, const boost::shared_ptr< urdf::Link >& child_link)
  : 
  m_parent_link(parent_link)  
  {
    
    ROS_DEBUG_STREAM("Adding joint named: " <<  urdf_joint->name);
    m_T_pj = urdfPoseToAffine(urdf_joint->parent_to_joint_origin_transform);
    m_axis = urdfVectorToEigen(urdf_joint->axis);
    m_name = urdf_joint->name;
    
    if (urdf_joint->type == urdf::Joint::REVOLUTE)
      m_type = itia::dynamics::Joint::REVOLUTE;
    else if (urdf_joint->type == urdf::Joint::CONTINUOUS)
      m_type = itia::dynamics::Joint::REVOLUTE;
    else if (urdf_joint->type == urdf::Joint::PRISMATIC)
      m_type = itia::dynamics::Joint::PRISMATIC;
    else
      m_type = itia::dynamics::Joint::FIXED;
    
    m_child_link.reset(
      new itia::dynamics::Link(
        child_link, 
        pointer()) 
    );
  }*/

//   class Tree: public boost::enable_shared_from_this<itia::dynamics::Tree>
//   {
//   protected:
//     std::map<std::string, boost::shared_ptr<itia::dynamics::Link>>  m_links_map;
//     std::map<std::string, boost::shared_ptr<itia::dynamics::Joint>> m_joints_map;
//   public:
//     Tree(const urdf::Model)
//     {
//       
//     };
//     boost::shared_ptr<itia::dynamics::Link> pointer()
//     {
//       return shared_from_this();
//     }
//   };


}
}
# endif