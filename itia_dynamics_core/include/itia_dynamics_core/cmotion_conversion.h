
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

#ifndef __ITIA_DYN__QUATERNIONS__UTILS__
#define __ITIA_DYN__QUATERNIONS__UTILS__

# include <boost/math/special_functions/binomial.hpp>
# include <itia_mutils/norm_derivative.h>
# include <itia_dynamics_core/itia_primitives.h>

#include<Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Quaterniond)


namespace itia
{
namespace dynamics
{

inline Eigen::Affine3d VectorToAffine(const Eigen::VectorXd& vec)
{
  Eigen::Affine3d T;
  T.translation() =vec.block(0, 0, 3, 1);
  Eigen::Quaterniond q;
  q.coeffs() = vec.block(3, 0, 4, 1);
  T.linear() = q.toRotationMatrix();
  return T;
};

inline Eigen::Affine3d CMotionToAffine(const CMotion& motion, Twist* twists)
{
  unsigned int n_derivative = motion.cols()-1;
  if (n_derivative<0)
  {
    ROS_ERROR("derivative order = %d", n_derivative);
    throw std::out_of_range("negative derivative order");
  }
  twists->resize(6, n_derivative);
  twists->setZero();
  Eigen::Affine3d T;
  T.translation() =motion.block(0, 0, 3, 1);
  
  std::vector<Eigen::Quaterniond> qvec(n_derivative+1);
  
//   ROS_INFO_STREAM("motion \n" <<  motion.transpose()); 
  qvec.at(0).coeffs() = motion.block(3, 0, 4, 1);
  Eigen::Quaterniond qinv = qvec.at(0).inverse();
  T.linear() = qvec.at(0).toRotationMatrix();
  
  // Dq    = 0.5*w*q;
  // DDq   = 0.5*(Dw*q+w*Dq)
  // DDDq  = 0.5*(DDw*q+2*Dw*Dq+w*DDq)
  // DDDDq = 0.5*(DDDw*q+3*DDw*Dq+3*Dw*DDq+w*DDDq)
  
  // w      = (2*Dq )*qinv
  // Dw     = (2*DDq-w*Dq)*qinv = 2* DDq*qinv - w * (Dq*qinv) = 2* DDq*qinv - 0.5* w * w
  // DDw    = (2*DDDq-2*Dw*Dq-w*DDq)*qinv
  // DDDw   = (2*DDDDq-3*DDw*Dq-3*Dw*DDq-w*DDDq)*qinv

  std::vector<Eigen::Quaterniond> derw(n_derivative);
  
  
  for (unsigned int idx_der = 1;idx_der <= n_derivative;idx_der++)
  {
    qvec.at(idx_der).coeffs() = motion.block(3, idx_der, 4, 1);
    derw.at(idx_der-1).coeffs() =2*(qvec.at(idx_der) *qinv).coeffs();
//     ROS_INFO_STREAM("qvec \n" << qvec.at(idx_der).coeffs());
//     ROS_INFO_STREAM("qinv \n" << qinv.coeffs());
//     ROS_INFO_STREAM("derw \n" << derw.at(idx_der-1).coeffs());
    
    for (unsigned int idx = 1;idx<(idx_der);idx++)
    {
      double coeff = boost::math::binomial_coefficient<double>(idx_der-1, idx);
      int idx_1 = idx_der-idx-1;
      derw.at(idx_der-1).coeffs()-= coeff * (derw.at(idx_1) * qvec.at(idx) *qinv).coeffs();
    }
//     ROS_INFO_STREAM("derw \n" << derw.at(idx_der-1).coeffs());
    twists->block(0, idx_der-1, 3, 1) = motion.block(0, idx_der, 3, 1);
    twists->block(3, idx_der-1, 3, 1) = derw.at(idx_der-1).coeffs().block(0, 0, 3, 1);
  }
  
  return T;
};

inline Eigen::VectorXd AffineToVector(const Eigen::Affine3d& T)
{
  Eigen::VectorXd pose(7);
  pose.block(0, 0, 3, 1) =T.translation();
  pose.block(3, 0, 4, 1) =Eigen::Quaterniond(T.linear()).coeffs();
  return pose;
}

inline itia::CMotion AffineToCMotion(const Eigen::Affine3d& T, const itia::Twist& twists)
{
  unsigned int n_derivative = twists.cols();
  itia::CMotion motion(7, n_derivative+1);
  motion.col(0) = AffineToVector(T);
  
  std::vector<Eigen::Quaterniond> qvec(n_derivative+1);
  std::vector<Eigen::Quaterniond> derw(n_derivative);
  qvec.at(0).coeffs() = motion.block(3, 0, 4, 1);
  
  // Dq    = 0.5*w*q;
  // DDq   = 0.5*(Dw*q+w*Dq)
  // DDDq  = 0.5*(DDw*q+2*Dw*Dq+w*DDq)
  // DDDDq = 0.5*(DDDw*q+3*DDw*Dq+3*Dw*DDq+w*DDDq)
  
  for (unsigned int idx_der = 1;idx_der <= n_derivative;idx_der++)
  {
    derw.at(idx_der-1) = Eigen::Quaterniond( 0, twists(3, idx_der-1), twists(4, idx_der-1), twists(5, idx_der-1) );
    qvec.at(idx_der) =Eigen::Quaterniond( 0, 0, 0, 0);
    for (unsigned int idx = 0;idx<(idx_der);idx++)
    {
      double coeff = boost::math::binomial_coefficient<double>(idx_der-1, idx);
      int idx_1 = idx_der-idx-1;
      qvec.at(idx_der).coeffs()+= coeff * (derw.at(idx) * qvec.at(idx_1)).coeffs();
    }
    qvec.at(idx_der).coeffs() *=0.5;
    motion.block(0, idx_der, 3, 1) =twists.block(0, idx_der-1, 3, 1);
    motion.block(3, idx_der, 4, 1) =qvec.at(idx_der).coeffs();
  }
  
  
  return motion;
}

inline itia::CMotion PolyVectorToQuatVector(const itia::CMotion& motion, const double& toll = 1e-2)
{
  
  // q    =    p * 1/norm(p) 
  // Dq   =   Dp * 1/norm(p) +     p * D(1/norm(p))
  // DDq  =  DDp * 1/norm(p) +  2*Dp * D(1/norm(p)) +    p * DD(1/norm(p))
  // DDDq = DDDp * 1/norm(p) + 3*DDp * D(1/norm(p)) + 3*Dp * DD(1/norm(p)) + p * DDD(1/norm(p))
  
  // (p)^c   = (p'*p)^c
  // D(p)^c  = c*(p'*p)^(c-1) + (p'*p)^c * 2*p'*Dp
  // DD(p)^c = c*(p'*p)^(c-2) + D(p'*p)^c * 2*p'*Dp + (p'*p)^c * 2*Dp'*Dp + (p'*p)^c * 2*p'*DDp

  itia::CMotion quat_motion = motion;
  int nder = motion.cols()-1;
  
  Eigen::VectorXd p = motion.block(3, 0, 4, 1);
  
  double norm=p.norm();
  if (norm<toll){
    printf("[PolynomialInterpolator::normalizePolyValuve] quaternion norm too low!\n");
    norm=toll;
  }
  
  // q=p/norm;
  quat_motion.block(3, 0, 4, 1) = p/norm;
  if (nder>0)
  {
    Eigen::VectorXd dp = motion.block(3, 1, 4, 1);
    double dnorm=p.dot(dp)/norm;
    // Dq=( Dp*norm - p*Dnorm ) /norm^2
    quat_motion.block(3, 1, 4, 1) = (dp*norm-p*dnorm)/pow(norm, 2);
    if (nder>1)
    {
      Eigen::VectorXd ddp = motion.block(3, 2, 4, 1);
      double ddnorm=( (p.dot(ddp)+dp.dot(dp))*norm - (p.dot(dp))*dnorm )/pow(norm,2.0);
      // DDq=( (DDp*norm-p*DDnorm )*norm^2 - (Dp*norm - p*Dnorm)*2*norm*Dnorm)/norm^4
      quat_motion.block(3, 2, 4, 1) = ((ddp*norm-p*ddnorm) *pow(norm, 2) - (dp*norm - p*dnorm)*2*norm*dnorm)/pow(norm, 4);
      if (nder>2)
      {
        printf("[itia_mutils/quaternion.h] Not implemented yet\n");
        for (int idx = 2;idx<nder;idx++)
          quat_motion.block(3, idx+1, 4, 1).setZero();
      }
    }
    return quat_motion;
  }
  
};



}
}


#endif