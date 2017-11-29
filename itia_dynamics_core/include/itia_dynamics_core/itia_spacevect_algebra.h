
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

#ifndef __ITIA_SPACEVECT_ALGEBRA__
#define __ITIA_SPACEVECT_ALGEBRA__

# include <Eigen/Geometry>
# include <Eigen/Core>

namespace Eigen
{
  typedef Matrix<double, 6,  1> Vector6d;
  typedef Matrix<double, 6, -1> Matrix6Xd;
  typedef Matrix<double, 6, 6> Matrix66d;
}


namespace itia 
{
namespace dynamics 
{

inline Eigen::Matrix3d skew(const Eigen::Vector3d& vec)
{
  Eigen::Matrix3d mat;
  mat << 0, -vec(2), vec(1),
    vec(2), 0, -vec(0),
    -vec(1), vec(0), 0;
  return mat;
}

inline Eigen::Vector3d unskew(const Eigen::Matrix3d& mat)
{
  return Eigen::Vector3d(mat(2, 1), mat(0, 2), mat(1, 0));
}

/*
 * for twist
 */
inline void spatialCrossProduct(const Eigen::Vector6d& vet1, const Eigen::Vector6d& vet2, Eigen::Vector6d* res)
{
  res->block(3, 0, 3, 1) = ((Eigen::Vector3d)(vet1.block(3, 0, 3, 1))).cross((Eigen::Vector3d)vet2.block(3, 0, 3, 1));
  res->block(0, 0, 3, 1) = ((Eigen::Vector3d)(vet1.block(3, 0, 3, 1))).cross((Eigen::Vector3d)vet2.block(0, 0, 3, 1))+
                          ((Eigen::Vector3d)(vet1.block(0, 0, 3, 1))).cross((Eigen::Vector3d)vet2.block(3, 0, 3, 1));
}

/*
 * for twist
 */
inline Eigen::Vector6d spatialCrossProduct(const Eigen::Vector6d& vet1, const Eigen::Vector6d& vet2)
{
  Eigen::Vector6d res;
  spatialCrossProduct(vet1, vet2, &res);
  return res;
}

/*
 * for wrench
 */
inline void spatialDualCrossProduct(const Eigen::Vector6d& vet1, const Eigen::Vector6d& vet2, Eigen::Vector6d* res)
{
  res->block(3, 0, 3, 1) = ((Eigen::Vector3d)(vet1.block(3, 0, 3, 1))).cross((Eigen::Vector3d)vet2.block(3, 0, 3, 1))+
                           ((Eigen::Vector3d)(vet1.block(0, 0, 3, 1))).cross((Eigen::Vector3d)vet2.block(0, 0, 3, 1));
  res->block(0, 0, 3, 1) = ((Eigen::Vector3d)(vet1.block(3, 0, 3, 1))).cross((Eigen::Vector3d)vet2.block(0, 0, 3, 1));
}

/*
 * for wrench
 */
inline Eigen::Vector6d spatialDualCrossProduct(const Eigen::Vector6d& vet1, const Eigen::Vector6d& vet2)
{
  Eigen::Vector6d res;
  spatialDualCrossProduct(vet1, vet2, &res);
  return res;
}

/*
 * for twist
 */
inline void spatialTranslation(const Eigen::Vector6d& vet, const Eigen::Vector3d& tra, Eigen::Vector6d* res)
{
 (*res) = vet;
 res->block(0, 0, 3, 1) = vet.block(0, 0, 3, 1) + ((Eigen::Vector3d)(vet.block(3, 0, 3, 1))).cross(tra);
}

/*
 * for twist
 */
inline Eigen::Vector6d spatialTranslation(const Eigen::Vector6d& vet, const Eigen::Vector3d& tra)
{
  Eigen::Vector6d res;
  spatialTranslation(vet, tra, &res);
  return res;
}

/*
 * for wrench
 */
inline void spatialDualTranslation(const Eigen::Vector6d& vet, const Eigen::Vector3d& tra, Eigen::Vector6d* res)
{
 (*res) = vet;
 res->block(3, 0, 3, 1) = vet.block(3, 0, 3, 1) + ((Eigen::Vector3d)(vet.block(0, 0, 3, 1))).cross(tra);
}

/*
 * for wrench
 */
inline Eigen::Vector6d spatialDualTranslation(const Eigen::Vector6d& vet, const Eigen::Vector3d& tra)
{
  Eigen::Vector6d res;
  spatialDualTranslation(vet, tra, &res);
  return res;
}

/*
 * for twist and wrench
 */
inline void spatialRotation(const Eigen::Vector6d& vet, const Eigen::Matrix3d& rot, Eigen::Vector6d* res)
{
  (*res) << rot*vet.block(0, 0, 3, 1), rot*vet.block(3, 0, 3, 1);
}

/*
 * for twist and wrench
 */
inline Eigen::Vector6d spatialRotation(const Eigen::Vector6d& vet, const Eigen::Matrix3d& rot)
{
  Eigen::Vector6d res;
  spatialRotation(vet, rot, &res);
  return res;
}

/*
 * for twist
 */
inline void spatialTranformation(const Eigen::Vector6d& vet, const Eigen::Affine3d& T, Eigen::Vector6d* res)
{
  (*res) << T.linear()*vet.block(0, 0, 3, 1) + ((Eigen::Vector3d)(vet.block(3, 0, 3, 1))).cross(T.translation()), 
            T.linear()*vet.block(3, 0, 3, 1);
}

/*
 * for twist
 */
inline Eigen::Vector6d spatialTranformation(const Eigen::Vector6d& vet, const Eigen::Affine3d& T)
{
  Eigen::Vector6d res;
  spatialTranformation(vet, T, &res);
  return res;
}

/*
 * for wrench
 */
inline void spatialDualTranformation(const Eigen::Vector6d& vet, const Eigen::Affine3d& T, Eigen::Vector6d* res)
{
  (*res) << T.linear()*vet.block(0, 0, 3, 1), 
            T.linear()*vet.block(3, 0, 3, 1) + ((Eigen::Vector3d)(vet.block(0, 0, 3, 1))).cross(T.translation());
}

/*
 * for wrench
 */
inline Eigen::Vector6d spatialDualTranformation(const Eigen::Vector6d& vet, const Eigen::Affine3d& T)
{
  Eigen::Vector6d res;
  spatialDualTranformation(vet, T, &res);
  return res;
}

inline void computeSpatialInertiaMatrix(const Eigen::Ref<Eigen::MatrixXd>& inertia, const Eigen::Ref<Eigen::Vector3d> cog, const double& mass, Eigen::Ref<Eigen::Matrix<double, 6, 6>> spatial_inertia)
{
  Eigen::MatrixXd cog_skew = itia::dynamics::skew(cog);
  spatial_inertia.block(0, 0, 3, 3) = mass*Eigen::MatrixXd::Identity(3, 3);
  spatial_inertia.block(0, 3, 3, 3) = mass*cog_skew.transpose();

  spatial_inertia.block(3, 0, 3, 3) = mass*cog_skew;
  spatial_inertia.block(3, 3, 3, 3) = inertia + mass*(cog_skew*cog_skew.transpose());
  
}


}
}

# endif 
