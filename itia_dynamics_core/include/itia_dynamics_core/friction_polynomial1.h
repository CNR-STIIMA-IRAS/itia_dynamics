
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

#ifndef __ITIA_DYNAMICS_CORE__FRICTION_POLY1__
#define __ITIA_DYNAMICS_CORE__FRICTION_POLY1__

# include <itia_dynamics_core/base_component.h>


namespace itia 
{
namespace dynamics 
{

class FirstOrderPolynomialFriction: public ComponentBase
{
protected:
  Eigen::Matrix<double, 2, 1> m_friction_coefficients;
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
    
  }
  
public:
  FirstOrderPolynomialFriction(const std::string& joint_name, const std::string& robot_name): ComponentBase(joint_name, robot_name)
  {
    
    std::vector<double> coefficients;
    if (!nh.getParam(robot_name+"/" + joint_name + "/friction/coefficients" , coefficients))
      throw std::invalid_argument(robot_name+"/" + joint_name + "/friction/coefficients NOT FOUND");
    
    if (coefficients.size() != 4)
      throw std::invalid_argument(robot_name+"/" + joint_name + "/friction/coefficients has wrong dimensions, " + std::to_string(coefficients.size()));
    
    m_friction_coefficients(0) = coefficients.at(0);
    m_friction_coefficients(1) = coefficients.at(1);
    
    if (coefficients.at(2) <1e-6)
    {
      ROS_WARN("[ FirstOrderPolynomialFriction ] Dq_threshold (coefficient[3] must be greater than 1.0e-6, set default value = 1.0e-6");
      m_Dq_threshold = 1.0e-6;
    }
    else
      m_Dq_threshold = coefficients.at(2);
    
    if (coefficients.at(3) <= 0)
    {
      ROS_WARN("[ FirstOrderPolynomialFriction ] Dq_max (coefficient[4] must be stricly positive, set default value = 1.0e+6");
      m_Dq_max = 1.0e+6;
    }
    else
      m_Dq_max = coefficients.at(3);
    
    m_regressor.resize(m_joints_number, 2);
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
  
  virtual unsigned int getParametersNumber() {return 2;};
};

}
}



#endif