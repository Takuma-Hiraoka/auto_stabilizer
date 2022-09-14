#include "AutoStabilizerROSBridge.h"

AutoStabilizerROSBridge::AutoStabilizerROSBridge(RTC::Manager* manager):
  RTC::DataFlowComponentBase(manager),
  m_steppableRegionOut_("steppableRegionOut", m_steppableRegion_),
  m_landingTargetIn_("landingTargetIn", m_landingTarget_)
{
}

RTC::ReturnCode_t AutoStabilizerROSBridge::onInitialize(){
  addOutPort("steppableRegionOut", m_steppableRegionOut_);
  addInPort("landingTargetIn", m_landingTargetIn_);

  ros::NodeHandle pnh("~");

  return RTC::RTC_OK;
}

RTC::ReturnCode_t AutoStabilizerROSBridge::onExecute(RTC::UniqueId ec_id){
  ros::spinOnce();
  if(this->m_landingTargetIn_.isNew()){
    this->m_landingTargetIn_.read();
  }
  return RTC::RTC_OK;
}

static const char* AutoStabilizerROSBridge_spec[] = {
  "implementation_id", "AutoStabilizerROSBridge",
  "type_name",         "AutoStabilizerROSBridge",
  "description",       "AutoStabilizerROSBridge component",
  "version",           "0.0",
  "vendor",            "Takuma-Hiraoka",
  "category",          "example",
  "activity_type",     "DataFlowComponent",
  "max_instance",      "10",
  "language",          "C++",
  "lang_type",         "compile",
  ""
};

extern "C"{
    void AutoStabilizerROSBridgeInit(RTC::Manager* manager) {
        RTC::Properties profile(AutoStabilizerROSBridge_spec);
        manager->registerFactory(profile, RTC::Create<AutoStabilizerROSBridge>, RTC::Delete<AutoStabilizerROSBridge>);
    }
};
