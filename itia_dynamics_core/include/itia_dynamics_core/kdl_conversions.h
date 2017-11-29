
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

#ifndef __ITIA_DYN__KDL_CONV__
#define __ITIA_DYN__KDL_CONV__
#include "cmotion_conversion.h"
# include <eigen_conversions/eigen_kdl.h>
# include <eigen_conversions/eigen_msg.h>
# include <kdl_conversions/kdl_msg.h>
# include <kdl/frameacc.hpp>
# include <kdl/framevel.hpp>
# include <kdl/jntarray.hpp>
# include <kdl/jntarrayvel.hpp>
# include <kdl/jntarrayacc.hpp>

namespace itia
{
namespace dynamics 
{




inline void frameToCMotion(const KDL::Frame& frame, CMotion* cmotion)
{
  if (cmotion->cols() < 1)
    cmotion->resize(7, 1);
  cmotion->setZero();
  
  geometry_msgs::Pose pose_msg;
  Eigen::Affine3d pose;
  tf::poseKDLToMsg(frame, pose_msg);
  tf::poseMsgToEigen(pose_msg, pose);
  
  cmotion->col(0) = itia::dynamics::AffineToVector(pose);
};

inline void frameVelToCMotion(const KDL::FrameVel& frame_vel, CMotion* cmotion)
{
  cmotion->resize(7, 2);
  cmotion->setZero();
  
  itia::Twist twist;
  twist.resize(6, 1);
  twist.setZero();
  
  geometry_msgs::Pose pose_msg;
  Eigen::Affine3d pose;
  tf::poseKDLToMsg(frame_vel.GetFrame(), pose_msg);
  tf::poseMsgToEigen(pose_msg, pose);
  
  
  Eigen::Matrix<double, 6, 1> vel_twist;
  tf::twistKDLToEigen(frame_vel.GetTwist(), vel_twist);
  twist.col(0) = vel_twist;
  
  *cmotion = itia::dynamics::AffineToCMotion(pose, twist);
};

inline void frameAccToCMotion(const KDL::FrameAcc& frame_acc, CMotion* cmotion)
{
  cmotion->resize(7, 3);
  cmotion->setZero();
  
  itia::Twist twist;
  twist.resize(6, 2);
  twist.setZero();
  
  geometry_msgs::Pose pose_msg;
  Eigen::Affine3d pose;
  tf::poseKDLToMsg(frame_acc.GetFrame(), pose_msg);
  tf::poseMsgToEigen(pose_msg, pose);
  
  
  Eigen::Matrix<double, 6, 1> vel_twist;
  tf::twistKDLToEigen(frame_acc.GetTwist(), vel_twist);
  
  Eigen::Matrix<double, 6, 1> acc_twist;
  tf::twistKDLToEigen(frame_acc.GetAccTwist(), acc_twist);
  
  twist.col(0) = vel_twist;
  twist.col(1) = acc_twist;  
  *cmotion = itia::dynamics::AffineToCMotion(pose, twist);
}

inline void CMotionToFrame(const itia::CMotion cmotion, KDL::Frame* frame)
{
  itia::Twist twist;
  Eigen::Affine3d pose = itia::dynamics::CMotionToAffine(cmotion, &twist);
  geometry_msgs::Pose pose_msg;
  tf::poseEigenToMsg(pose, pose_msg);
  tf::poseMsgToKDL(pose_msg, *frame);
}

inline void CMotionToFrameVel(const itia::CMotion cmotion, KDL::FrameVel* frame_vel)
{
  itia::Twist twist;
  ROS_INFO_STREAM("twist (PRE)\n" <<  twist);
  Eigen::Affine3d pose = itia::dynamics::CMotionToAffine(cmotion, &twist);
  geometry_msgs::Pose pose_msg;
  tf::poseEigenToMsg(pose, pose_msg);
  KDL::Frame frame;
  tf::poseMsgToKDL(pose_msg, frame);
  frame_vel->M.R = frame.M;
  frame_vel->p.p = frame.p;
  
  ROS_INFO_STREAM("cmotion \n" << cmotion);
  ROS_INFO_STREAM("twist \n" <<  twist);
  twist.setZero();                                          // BUG
  tf::vectorEigenToKDL(twist.block(0, 0, 3, 1), frame_vel->p.v);
  tf::vectorEigenToKDL(twist.block(3, 0, 3, 1), frame_vel->M.w);
};

inline void CMotionToFrameAcc(const itia::CMotion cmotion, KDL::FrameAcc* frame_acc)
{
  itia::Twist twist;
  Eigen::Affine3d pose = itia::dynamics::CMotionToAffine(cmotion, &twist);
  geometry_msgs::Pose pose_msg;
  tf::poseEigenToMsg(pose, pose_msg);
  KDL::Frame frame;
  tf::poseMsgToKDL(pose_msg, frame);
  frame_acc->M.R = frame.M;
  frame_acc->p.p = frame.p;
  
  tf::vectorEigenToKDL(twist.block(0, 0, 3, 1), frame_acc->p.v);
  tf::vectorEigenToKDL(twist.block(3, 0, 3, 1), frame_acc->M.w);
  tf::vectorEigenToKDL(twist.block(0, 1, 3, 1), frame_acc->p.dv);
  tf::vectorEigenToKDL(twist.block(3, 1, 3, 1), frame_acc->M.dw);
};

inline void JointToJMotion(const KDL::JntArray& joint, itia::JMotion* jmotion)
{
  jmotion->resize(joint.data.rows(), 1);
  jmotion->col(0) = joint.data;
}

inline void JointToJMotion(const KDL::JntArrayVel& joint, itia::JMotion* jmotion)
{
  jmotion->resize(joint.q.data.rows(), 2);
  jmotion->col(0) = joint.q.data;
  jmotion->col(1) = joint.qdot.data;
}

inline void JointToJMotion(const KDL::JntArrayAcc& joint, itia::JMotion* jmotion)
{
  jmotion->resize(joint.q.data.rows(), 3);
  jmotion->col(0) = joint.q.data;
  jmotion->col(1) = joint.qdot.data;
  jmotion->col(2) = joint.qdotdot.data;
}

inline void JMotionToJoint(const itia::JMotion& jmotion, KDL::JntArray* joint)
{
  joint->data = jmotion.col(0);
}

inline void JMotionToJoint(const itia::JMotion& jmotion, KDL::JntArrayVel* joint)
{
  joint->q.data = jmotion.col(0);
  joint->qdot.data = jmotion.col(1);
}

inline void JMotionToJoint(const itia::JMotion& jmotion, KDL::JntArrayAcc* joint)
{
  joint->q.data = jmotion.col(0);
  joint->qdot.data = jmotion.col(1);
  joint->qdotdot.data = jmotion.col(2);
}

inline void WrenchToKDLWrench(const itia::Wrench& wrench, KDL::Wrench* kdl_wrench)
{
  
  for (int idx = 0;idx<6;idx++)
    (*kdl_wrench)[idx] = wrench(idx);
}

}
}

# endif