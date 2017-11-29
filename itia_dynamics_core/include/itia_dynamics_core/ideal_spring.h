#ifndef __ITIA_DYNAMICS_CORE__SPRING__
#define __ITIA_DYNAMICS_CORE__SPRING__

# include <itia_dynamics_core/base_component.h>


namespace itia 
{
namespace dynamics 
{

class IdealSpring: public ComponentBase
{
protected:
  double m_elasticity;
  double m_offset;
  
public:
  IdealSpring(const std::string& joint_name, const std::string& robot_name): ComponentBase(joint_name, robot_name)
  {
    
    std::vector<double> coefficients;
    if (!nh.getParam(robot_name+"/" + joint_name + "/spring/coefficients" , coefficients))
      throw std::invalid_argument(robot_name+"/" + joint_name + "/spring/coefficients NOT FOUND");
    
    if (coefficients.size() != 2)
      throw std::invalid_argument(robot_name+"/" + joint_name + "/spring/coefficients has wrong dimensions, " +  std::to_string(coefficients.size()));

    m_elasticity   = coefficients.at(0);
    m_offset = coefficients.at(1);
    
    m_regressor.resize(m_joints_number, 2);
    m_regressor.setZero();
  };
  virtual Eigen::VectorXd getTorque(const Eigen::Ref<Eigen::VectorXd>& q,  const Eigen::Ref<Eigen::VectorXd>& Dq, const Eigen::Ref<Eigen::VectorXd>& DDq)
  {
    m_torques(m_component_joint_number) = m_elasticity*q(m_joints_number)+m_offset;
    return m_torques;
  };
  virtual Eigen::MatrixXd getRegressor(const Eigen::Ref<Eigen::VectorXd>& q,  const Eigen::Ref<Eigen::VectorXd>& Dq, const Eigen::Ref<Eigen::VectorXd>& DDq)
  {
    m_regressor(m_component_joint_number, 0) = q(m_component_joint_number);
    m_regressor(m_component_joint_number, 1) = 1;
    return m_regressor;
  }
  
  virtual unsigned int getParametersNumber() {return 2;};
};

}
}

#endif