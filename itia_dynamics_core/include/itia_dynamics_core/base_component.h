#ifndef __ITIA_DYNAMICS_CORE__COMPONENT__
#define __ITIA_DYNAMICS_CORE__COMPONENT__

# include <boost/graph/graph_concepts.hpp>
# include <Eigen/Core>
# include <ros/ros.h>

namespace itia 
{
namespace dynamics 
{

class ComponentBase
{
protected:
  
  // Name of the robot
  std::string m_robot_name;
  
  // Name of the Joint where the component is applied
  std::string m_component_joint_name;
  
  // Name of the joints of the tree
  std::vector<std::string>  m_joint_names;
  
  // Number of Joints
  unsigned int m_joints_number;
  
  // Number of the joint component in the list of joint names
  unsigned int m_component_joint_number;
  
  Eigen::VectorXd m_torques;
  Eigen::MatrixXd m_regressor;
  ros::NodeHandle nh;
  
public:
  ComponentBase(const std::string& joint_name, const std::string& robot_name)
  {
    m_component_joint_name= joint_name;
    m_robot_name = robot_name;
    
    if (!nh.getParam(robot_name+"/joint_names", m_joint_names))
    {
      throw std::invalid_argument("PARAMETER '" + robot_name+"/joint_names' NOT FOUND");
    }
    m_joints_number = m_joint_names.size();
    
    unsigned int idx = 0;
    for (; idx<m_joints_number; idx++)
      if ( !m_joint_names.at(idx).compare(m_component_joint_name) )
      {
        m_component_joint_number = idx;
        break;
      }
      
    if (idx == m_joints_number)
      throw std::invalid_argument("Component Joint name '" +m_component_joint_name+"' is not a elemente of '" + robot_name+"/joint_names'");
      
    m_torques.resize(m_joints_number);
    m_torques.setZero();
  };
  virtual Eigen::VectorXd getTorque(const Eigen::Ref<Eigen::VectorXd>& q,  const Eigen::Ref<Eigen::VectorXd>& Dq, const Eigen::Ref<Eigen::VectorXd>& DDq) = 0;
  virtual Eigen::MatrixXd getRegressor(const Eigen::Ref<Eigen::VectorXd>& q,  const Eigen::Ref<Eigen::VectorXd>& Dq, const Eigen::Ref<Eigen::VectorXd>& DDq) = 0;
  virtual unsigned int getParametersNumber() = 0;
  
};


typedef  boost::shared_ptr<itia::dynamics::ComponentBase> ComponentPtr;
}
}

#endif