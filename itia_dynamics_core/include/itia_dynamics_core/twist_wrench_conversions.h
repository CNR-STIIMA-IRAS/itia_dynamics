#ifndef __ITIA_DYN__TWIST_CONV__
#define __ITIA_DYN__TWIST_CONV__
# include <itia_dynamics_core/itia_primitives.h>

namespace itia
{
namespace dynamics
{


/* 
 * Twist changeRefPointTwists(const Twist twists_a_in_c, const Eigen::Affine3d& Q_b_a_in_c)
 * compute twist of body b in frame c starting from twist of body a in frame c and the rigid trasformation from a to b in frame c
 * T_b_in_c <- T_b_in_c
 */
inline Twist changeRefPointTwists(const Twist twists_a_in_c, const Eigen::Affine3d& Q_b_a_in_c)
{
  // w_b_in_c = w_a_in_c  
  // v_b_in_c = v_a_in_c + w_a_in_c X distance_ba_in_c
  
  Twist twists_b_in_c = twists_a_in_c;
  Eigen::Vector3d dist_ab_in_c = - Q_b_a_in_c.translation();
  
  for (int idx_c = 0;idx_c < twists_a_in_c.cols();idx_c++)
  {
    Eigen::Vector3d tmp = twists_a_in_c.block(3, idx_c, 3, 1);
    twists_b_in_c.block(0, idx_c, 3, 1) = twists_a_in_c.block(0, idx_c, 3, 1) + tmp.cross(dist_ab_in_c);
  }
  return twists_b_in_c;
}

/*
 * Wrench changeRefPointWrenches(const Wrench wrenches_a_in_c, const Eigen::Affine3d& Q_b_a_in_c)
 * compute wrench applied to body b in frame c starting from wrench applied to body a in frame c and the rigid trasformation from a to b in frame c
 * W_b_in_c <- W_a_in_c
 */
inline Wrench changeRefPointWrenches(const Wrench wrenches_a_in_c, const Eigen::Affine3d& Q_b_a_in_c)
{
  // t_b_in_c = t_a_in_c + distance_ab_in_c X f_a_in_c ->
  // t_b_in_c = t_a_in_c + f_a_in_c X distance_ba_in_c
  // f_b_in_c = f_a_in_c
  
  Wrench wrenches_b_in_c = wrenches_a_in_c;
  Eigen::Vector3d dist_ab_in_c = - Q_b_a_in_c.translation();
  
  for (int idx_c = 0;idx_c < wrenches_a_in_c.cols();idx_c++)
  {
    Eigen::Vector3d tmp = wrenches_a_in_c.block(0, idx_c, 3, 1);
    wrenches_b_in_c.block(3, idx_c, 3, 1) = wrenches_a_in_c.block(3, idx_c, 3, 1) + tmp.cross(dist_ab_in_c);
  }
  return wrenches_b_in_c;
}

/*
 * Twist changeFrameTwists(const Twist twists_a_in_a, const Eigen::Matrix3d& R_b_a)
 * T_a_in_b <- T_a_in_a
 */
inline Twist changeFrameTwists(const Twist twists_a_in_a, const Eigen::Matrix3d& R_b_a)
{
  // w_a_in_b = R_ba * w_a_in_a  
  // v_a_in_b = R_ba * v_a_in_a  
  
  Twist twists_a_in_b = twists_a_in_a;
  for (int idx_c = 0;idx_c < twists_a_in_a.cols();idx_c++)
  {
    twists_a_in_b.block(0, idx_c, 3, 1) = R_b_a * twists_a_in_a.block(0, idx_c, 3, 1);
    twists_a_in_b.block(3, idx_c, 3, 1) = R_b_a * twists_a_in_a.block(3, idx_c, 3, 1);
  }
  return twists_a_in_b;
}

/*
 * Wrench changeFrameWrenches(const Wrench wrenches_a_in_a, const Eigen::Matrix3d& R_b_a)
 * W_a_in_b <- W_a_in_a
 */
inline Wrench changeFrameWrenches(const Wrench wrenches_a_in_a, const Eigen::Matrix3d& R_b_a)
{
  // t_a_in_b = R_ba * t_a_in_a  
  // f_a_in_b = R_ba * f_a_in_a  
  Wrench wrenches_a_in_b = wrenches_a_in_a;
  for (int idx_c = 0;idx_c < wrenches_a_in_a.cols();idx_c++)
  {
    wrenches_a_in_b.block(0, idx_c, 3, 1) = R_b_a * wrenches_a_in_a.block(0, idx_c, 3, 1);
    wrenches_a_in_b.block(3, idx_c, 3, 1) = R_b_a * wrenches_a_in_a.block(3, idx_c, 3, 1);
  }
  return wrenches_a_in_b;
}


/*
 * Wrench changeFrameWrenches(const Wrench wrenches_a_in_a, const Eigen::Affine3d& T_b_a)
 * W_b_in_b <- W_a_in_a
 */
inline Wrench changeFrameWrenches(const Wrench wrenches_a_in_a, const Eigen::Affine3d& T_b_a)
{
  // Q_ba_in_c = T_cb * T_ac
  // Q_ba_in_b = T_bb * T_ab
  T_b_a.linear();
  Wrench wrenches_a_in_b = changeFrameWrenches(wrenches_a_in_a, T_b_a.linear());
  Wrench wrenches_b_in_b = changeRefPointWrenches(wrenches_a_in_b, T_b_a);
  return wrenches_b_in_b;
}

/*
 * Twist changeFrameTwists(const Twist wrenches_a_in_a, const Eigen::Affine3d& T_b_a)
 * W_b_in_b <- W_a_in_a
 */
inline Twist changeFrameTwists(const Twist twists_a_in_a, const Eigen::Affine3d& T_b_a)
{
  // Q_ba_in_c = T_cb * T_ac
  // Q_ba_in_b = T_bb * T_ab
  T_b_a.linear();
  Twist twists_a_in_b = changeFrameTwists(twists_a_in_a, T_b_a.linear());
  Twist twists_b_in_b = changeRefPointTwists(twists_a_in_b, T_b_a);
  return twists_b_in_b;
}

}
}
# endif