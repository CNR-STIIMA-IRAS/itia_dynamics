#ifndef __ITIA_DYNAMICS_CORE__FRICTION_POLY2__
#define __ITIA_DYNAMICS_CORE__FRICTION_POLY2__

# include <itia_dynamics_core/base_component.h>


namespace itia 
{
namespace dynamics 
{

class SecondOrderPolynomialFriction: public ComponentBase
{
protected:
  Eigen::Matrix<double, 3, 1> m_friction_coefficients;
  double m_Dq_threshold;
  double m_Dq_max;
  
  void computeRegressor(const Eigen::Ref<Eigen::VectorXd>& Dq )
  {
    double omega = std::min(std::max(Dq(m_component_joint_number), -m_Dq_max), m_Dq_max);
    double sign_Dq = 0;
    if (omega == 0)
      sign_Dq = 0;
    else if (omega > m_Dq_threshold)
      sign_Dq = 1.0;
    else if (omega < -m_Dq_threshold)
      sign_Dq = -1.0;
    else
      sign_Dq = omega/m_Dq_threshold;
    
    m_regressor(m_component_joint_number, 0) = sign_Dq;
    m_regressor(m_component_joint_number, 1) = omega;
    m_regressor(m_component_joint_number, 2) = pow(omega, 2.0) *sign_Dq;
    
  }
  
public:
  SecondOrderPolynomialFriction(const std::string& joint_name, const std::string& robot_name): ComponentBase(joint_name, robot_name)
  {
    
    std::vector<double> coefficients;
    if (!nh.getParam(robot_name+"/" + joint_name + "/friction/coefficients" , coefficients))
      throw std::invalid_argument(robot_name+"/" + joint_name + "/friction/coefficients NOT FOUND");
    
    if (coefficients.size() != 5)
      throw std::invalid_argument(robot_name+"/" + joint_name + "/friction/coefficients has wrong dimensions");
    
    m_friction_coefficients(0) = coefficients.at(0);
    m_friction_coefficients(1) = coefficients.at(1);
    m_friction_coefficients(2) = coefficients.at(2);
    
    if (coefficients.at(3) <1e-6)
    {
      ROS_WARN("[ SecondOrderPolynomialFriction ] Dq_threshold (coefficient[3] must be greater than 1.0e-6, set default value = 1.0e-6");
      m_Dq_threshold = 1.0e-6;
    }
    else
      m_Dq_threshold = coefficients.at(3);
    
    if (coefficients.at(4) <= 0)
    {
      ROS_WARN("[ SecondOrderPolynomialFriction ] Dq_max (coefficient[4] must be stricly positive, set default value = 1.0e+6");
      m_Dq_max = 1.0e+6;
    }
    else
      m_Dq_max = coefficients.at(4);
    
    m_regressor.resize(m_joints_number, 3);
    m_regressor.setZero();
    
    
    
  };
    
  virtual Eigen::VectorXd getTorque(const Eigen::Ref<Eigen::VectorXd>& q,  const Eigen::Ref<Eigen::VectorXd>& Dq, const Eigen::Ref<Eigen::VectorXd>& DDq)
  {
    computeRegressor(Dq);
    m_torques(m_component_joint_number) = m_regressor.row(m_component_joint_number)*m_friction_coefficients;
    return m_torques;
  };
  virtual Eigen::MatrixXd getRegressor(const Eigen::Ref<Eigen::VectorXd>& q,  const Eigen::Ref<Eigen::VectorXd>& Dq, const Eigen::Ref<Eigen::VectorXd>& DDq)
  {
    computeRegressor(Dq);
    return m_regressor;
  }
  
  virtual unsigned int getParametersNumber() {return 3;};
};

}
}



#endif