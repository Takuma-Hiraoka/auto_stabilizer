#ifndef AutoStabilizerROSBridge_H
#define AUtoStabilizerROSBridge_H

#include <rtm/Manager.h>
#include <rtm/DataFlowComponentBase.h>
#include <rtm/DataOutPort.h>
#include <rtm/DataInPort.h>

#include <auto_stabilizer_msgs/idl/AutoStabilizer.hh>
#include <auto_stabilizer_msgs/SteppableRegion.h>
#include <auto_stabilizer_msgs/LandingPosition.h>

#include <ros/ros.h>

class AutoStabilizerROSBridge : public RTC::DataFlowComponentBase{
protected:
  ros::NodeHandle nh;

  auto_stabilizer_msgs::TimedLandingPosition m_landingTarget_;
  RTC::InPort <auto_stabilizer_msgs::TimedLandingPosition> m_landingTargetIn_;
  ros::Publisher landing_target_pub_;

  ros::Subscriber steppable_region_sub_;
  ros::Subscriber landing_pose_sub_;
  auto_stabilizer_msgs::TimedSteppableRegion m_steppableRegion_;
  RTC::OutPort <auto_stabilizer_msgs::TimedSteppableRegion> m_steppableRegionOut_;
  auto_stabilizer_msgs::TimedLandingPosition m_landingPose_;
  RTC::OutPort <auto_stabilizer_msgs::TimedLandingPosition> m_landingPoseOut_;
public:
  AutoStabilizerROSBridge(RTC::Manager* manager);
  virtual RTC::ReturnCode_t onInitialize();
  virtual RTC::ReturnCode_t onExecute(RTC::UniqueId ec_id);

  void onSteppableRegionCB(const auto_stabilizer_msgs::SteppableRegion::ConstPtr& msg);
  void onLandingPoseCB(const auto_stabilizer_msgs::LandingPosition::ConstPtr& msg);
  
};

extern "C"
{
  void AutoStabilizerROSBridgeInit(RTC::Manager* manager);
};

#endif // AutoStabilizerROSBridge_H
