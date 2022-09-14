#ifndef AutoStabilizerROSBridge_H
#define AUtoStabilizerROSBridge_H

#include <rtm/Manager.h>
#include <rtm/DataFlowComponentBase.h>
#include <rtm/DataOutPort.h>
#include <rtm/DataInPort.h>

#include <auto_stabilizer_msgs/idl/AutoStabilizer.hh>

#include <ros/ros.h>

class AutoStabilizerROSBridge : public RTC::DataFlowComponentBase{
protected:
  ros::NodeHandle nh;

  auto_stabilizer_msgs::TimedLandingPosition m_landingTarget_;
  RTC::InPort <auto_stabilizer_msgs::TimedLandingPosition> m_landingTargetIn_;
  ros::Publisher pub_;

  ros::Subscriber sub_;
  auto_stabilizer_msgs::TimedSteppableRegion m_steppableRegion_;
  RTC::OutPort <auto_stabilizer_msgs::TimedSteppableRegion> m_steppableRegionOut_;
public:
  AutoStabilizerROSBridge(RTC::Manager* manager);
  virtual RTC::ReturnCode_t onInitialize();
  virtual RTC::ReturnCode_t onExecute(RTC::UniqueId ec_id);
  
};

extern "C"
{
  void AutoStabilizerROSBridgeInit(RTC::Manager* manager);
};

#endif // AutoStabilizerROSBridge_H
