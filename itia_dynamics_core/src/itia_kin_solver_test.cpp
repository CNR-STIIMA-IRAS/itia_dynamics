#include <itia_dynamics_core/itia_kin_solver_kdl.h>
#include <sensor_msgs/JointState.h>
#include <geometry_msgs/PoseStamped.h>

#include <itia_dynamics_core/itia_urdf_parser.h>

int main(int argc, char **argv){
  ros::init(argc, argv, "itia_kin_solver_test");
  ros::NodeHandle nh;
  ros::Rate rate(1);
  
  itia::dynamics::KinSolverKDL kin_solver("robot_description");
  
  std::string base_frame = "base_link";
  std::string tool_frame = "ee_link";                 //"ee_link";
  
  sensor_msgs::JointState js;
  js.name.resize(6);
  js.position.resize(6);
  js.velocity.resize(6);
  js.name.at(0) = "shoulder_pan_joint";
  js.name.at(1) = "shoulder_lift_joint";
  js.name.at(2) = "elbow_joint";
  js.name.at(3) = "wrist_1_joint";
  js.name.at(4) = "wrist_2_joint";
  js.name.at(5) = "wrist_3_joint";
  
  urdf::Model model;
  model.initParam("robot_description");
  
  
  ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::Level::Debug);
  boost::shared_ptr<itia::dynamics::Link> root_link(new itia::dynamics::Link());  
  root_link->fromUrdf(model.root_link_);
  boost::shared_ptr<itia::dynamics::Chain> chain(new itia::dynamics::Chain(root_link, base_frame,tool_frame));
  
  chain->setInputJointsName(js.name);
  
  Eigen::VectorXd q(6);
  q.setZero();
  Eigen::VectorXd Dq(6);
  Dq.setZero();
  
  
  ros::Publisher js_pub = nh.advertise<sensor_msgs::JointState>("joint_states", 1);
  ros::Publisher pose_pub = nh.advertise<geometry_msgs::PoseStamped>("pose", 1);
  
  itia::JMotion jmotion;
  jmotion.resize(6, 2);
  jmotion.setZero();
  
  itia::CMotion cmotion;
  
 
  
  itia::Twist twist;
  
  Eigen::Affine3d pose;
  Eigen::Affine3d pose1;
  
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  double t4 = 0;
  
  ros::Time t0;
  int ntrial = 1;
  
  Eigen::Matrix6Xd jac1, jac2;
  jac1.resize(6, 6);
  jac2.resize(6, 6);
  
  boost::shared_ptr<itia::dynamics::Joint> jptr = root_link->findChildJoint(js.name.at(0));
  for (int idx = 0;idx<ntrial;idx++)
  {
    
    jmotion.col(0) = M_PI*Eigen::Matrix<double, 6, 1>::Random();
//     jmotion.col(0) = M_PI*Eigen::Matrix<double, 6, 1>::Constant(6, 1, 0);
    q = jmotion.col(0);
     
    t0 = ros::Time::now();
    pose1 = chain->getTransformation(q);
    t1 +=  (ros::Time::now()-t0).toSec() *1000;
    
    t0 = ros::Time::now();
    cmotion = kin_solver.fkine(jmotion, tool_frame, base_frame);
    pose = itia::dynamics::CMotionToAffine(cmotion, &twist);
    t2 +=  (ros::Time::now()-t0).toSec() *1000;
    
    t0 = ros::Time::now();
    jac1 = chain->getJacobian(q);
    t3 +=  (ros::Time::now()-t0).toSec() *1000;
    
    t0 = ros::Time::now();
    jac2 = kin_solver.jac(jmotion, tool_frame, base_frame);
    t4 +=  (ros::Time::now()-t0).toSec() *1000;
    
  }
  ROS_INFO("kdl computation time = %5.4f [ms]",t2/ntrial);
  ROS_INFO("    computation time = %5.4f [ms]",t1/ntrial);
  ROS_INFO("jac    computation time = %5.4f [ms]",t3/ntrial);
  ROS_INFO("jacKDL computation time = %5.4f [ms]",t4/ntrial);
  
  ROS_INFO_STREAM("      T_bt = \n" <<  pose1.matrix()) ;
  ROS_INFO_STREAM("(KDL) T_bt = \n" <<  pose.matrix());
  
  ROS_INFO_STREAM("      J = \n" <<  jac1) ;
  ROS_INFO_STREAM("(KDL) J = \n" <<  jac2);
  
  ROS_INFO("norm(J-J)=%e", (jac1-jac2).norm());
  
  return 0;
  for (int idx1 = 0;idx1<1e4;idx1++)
  {
    q.setRandom();
    chain->getTransformation(q);
    Dq.setRandom();
    std::vector< Eigen::Vector6d,Eigen::aligned_allocator<Eigen::Vector6d> > twists = chain->getTwist(q, Dq);
    double err = (twists.back()-chain->getJacobian(q) *Dq).norm();
    ROS_WARN_STREAM_COND(err>1e-12, "\n error = " << err);
  }
  return 0;
  
  Eigen::Affine3d mtx;
  Eigen::VectorXd comau_pose;                               // metri /radianti
  comau_pose.resize(6);
  comau_pose << 0.5, 0.5, 0.3, M_PI/4, M_PI/2.0, 0;
  mtx = Eigen::AngleAxisd(comau_pose(3), Eigen::Vector3d::UnitZ()) *
        Eigen::AngleAxisd(comau_pose(4), Eigen::Vector3d::UnitY()) *
        Eigen::AngleAxisd(comau_pose(5), Eigen::Vector3d::UnitZ());
  
  mtx.translation()(0) =comau_pose(0);
  mtx.translation()(1) =comau_pose(1);
  mtx.translation()(2) =comau_pose(2);
  
  twist.resize(6,1);
  twist.setZero();
  cmotion = itia::dynamics::AffineToCMotion(mtx, twist);
  
  
  int iter = 0;
  while (ros::ok())
  {
    if (iter>5)
    {
      jmotion.col(0) = M_PI*Eigen::Matrix<double, 6, 1>::Random();
      ROS_INFO_STREAM("Random " << jmotion.col(0).transpose()); 
      iter = 0;
    }
    iter++;
//     itia::CMotion cmotion = kin_solver.fkine(jmotion, tool_frame, base_frame);
//     
    itia::JMotion jmotion_rand = jmotion;
    jmotion_rand.col(0) += 0.1*Eigen::Matrix<double, 6, 1>::Random();
    itia::JMotion jmotion2;
    
    pose = itia::dynamics::CMotionToAffine(cmotion, &twist);
    geometry_msgs::PoseStamped pose_msg;
    pose_msg.header.frame_id = base_frame;
    tf::poseEigenToMsg(pose, pose_msg.pose);
    
    int res;
    jmotion2 = kin_solver.ikine(cmotion, tool_frame, base_frame, jmotion_rand, &res);
    if (res == 0)
    {
      js.header.stamp = ros::Time::now();
      pose_msg.header.stamp = ros::Time::now();
      
      for (int idx = 0;idx<6;idx++)
      {
        js.position.at(idx) =jmotion2(idx, 0);
        js.velocity.at(idx) =jmotion2(idx, 1);
        
      }
      js_pub.publish(js);
      pose_pub.publish(pose_msg);
    }
    else 
    {
      ROS_ERROR("No solutions found");
    }
    rate.sleep();
  }
  return 0;  
}