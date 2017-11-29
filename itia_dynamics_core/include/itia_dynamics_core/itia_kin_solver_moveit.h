
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

#ifndef __ITIA_DYNAMICS_KIN_SOLVER_MOVEIT__
#define __ITIA_DYNAMICS_KIN_SOLVER_MOVEIT__
# include <itia_dynamics_core/itia_primitives.h>
# include <ros/ros.h>

# include <moveit/robot_model/robot_model.h>
# include <moveit/robot_state/robot_state.h>
# include <itia_dynamics_core/itia_kin_solver_base.h>


namespace itia
{
namespace dynamics 
{
  class KinSolverMoveit : public KinSolverBase
  {
  public:
    KinSolverMoveit(moveit::core::RobotModelPtr& robot_model);
    virtual itia::CMotion fkine(const itia::JMotion& jmotion, const std::string& tool_frame, const std::string& base_frame) {};
    virtual itia::JMotion ikine(const itia::CMotion& cmotion, const std::string& tool_frame, const std::string& base_frame,  itia::JMotion jmotion_0, int* result) {};
    
  protected:
    moveit::core::RobotModelPtr& m_robot_model;
    
    
    moveit::core::JointModelGroup* m_jmg;
//     kinematics::KinematicsBase m_kb;
    robot_state::RobotStatePtr m_kinematic_state;
  };
  
}
}






# endif