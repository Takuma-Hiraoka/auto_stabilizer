#include "Stabilizer.h"
#include "MathUtil.h"
#include <cnoid/Jacobian>
#include <cnoid/EigenUtil>
#include <cnoid/src/Body/InverseDynamics.h>

void Stabilizer::initStabilizerOutput(const GaitParam& gaitParam,
                                      cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Vector3& o_stTargetZmp, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage) const{

  o_stOffsetRootRpy.reset(cnoid::Vector3::Zero());
  o_stTargetZmp = gaitParam.refZmpTraj[0].getStart();
  for(int i=0;i<o_stServoPGainPercentage.size();i++){
    o_stServoPGainPercentage[i].reset(100.0);
    o_stServoDGainPercentage[i].reset(100.0);
  }
}

bool Stabilizer::execStabilizer(const GaitParam& gaitParam, double dt, bool useActState,
                                cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose) const{
  // - root attitude control
  // - 現在のactual重心位置から、目標ZMPを計算
  // - 目標ZMPを満たすように目標足裏反力を計算
  // - 目標反力を満たすように重力補償+仮想仕事の原理

  // root attitude control
  this->moveBasePosRotForBodyRPYControl(dt, gaitParam, useActState,// input
                                        o_stOffsetRootRpy, o_stTargetRootPose); // output

  return true;
}

bool Stabilizer::calcResolvedAccelationControl(const GaitParam& gaitParam, double dt, bool useActState, cnoid::BodyPtr& actRobotTqc, 
				     cnoid::Vector3& o_stTargetZmp, std::vector<cnoid::Vector6>& o_stEETargetWrench,
				     std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage,
					       Eigen::VectorXd& prev_q, Eigen::VectorXd& prev_dq, std::vector<cnoid::Vector6>& eePoseDiff_prev, std::vector<cnoid::Position>& eeTargetPosed, std::vector<cnoid::Position>& eeTargetPosedd, cnoid::Vector6& prev_rootd) const{
  // - 現在のactual重心位置から、目標ZMPを計算
  // - 目標位置姿勢を満たすように分解加速度制御. 重心が倒立振子で加速された場合のトルクを計算
  // - 目標のZMPを満たすように目標足裏反力を計算、仮想仕事の原理で足し込む

  // - 現在のactual重心位置から、目標ZMPを計算
  cnoid::Vector3 tgtForce; // generate frame
  cnoid::Vector3 tgtCogAcc; // generate frame
  this->calcZMP(gaitParam, dt, useActState, // input
                o_stTargetZmp, tgtForce, tgtCogAcc); // output

  cnoid::Vector6 tgtWrench; // generate frame
  this->calcTorque(dt, gaitParam, useActState, actRobotTqc, tgtCogAcc,
		   o_stServoPGainPercentage, o_stServoDGainPercentage,
		   tgtWrench,
		   prev_q, prev_dq, eePoseDiff_prev, eeTargetPosed, eeTargetPosedd, prev_rootd);

  // 目標ZMPを満たすように目標EndEffector反力を計算
  this->calcWrench(gaitParam, o_stTargetZmp, tgtWrench, useActState,// input
                   actRobotTqc, o_stEETargetWrench); // output

  // 出力トルク制限
  // ラチェッティングトルクはこれよりも小さい
  for(int i=0;i<actRobotTqc->numJoints();i++){
    mathutil::clamp(actRobotTqc->joint(i)->u(), this->torque_limit[i]);
  }

  return true;
};

bool Stabilizer::moveBasePosRotForBodyRPYControl(double dt, const GaitParam& gaitParam, bool useActState,
                                                 cpp_filters::TwoPointInterpolator<cnoid::Vector3>& o_stOffsetRootRpy, cnoid::Position& o_stTargetRootPose) const{

  // stOffsetRootRpyを計算
  if(!useActState){
      o_stOffsetRootRpy.interpolate(dt);
  }else{
    cnoid::Vector3 stOffsetRootRpy = o_stOffsetRootRpy.value(); // gaitParam.footMidCoords frame

    cnoid::Matrix3 rootRErrorGenerateFrame = gaitParam.refRobot->rootLink()->R() * gaitParam.actRobot->rootLink()->R().transpose(); // generate frame
    cnoid::Matrix3 rootRError = gaitParam.footMidCoords.value().linear().transpose() * rootRErrorGenerateFrame/*generate frame*/ * gaitParam.footMidCoords.value().linear(); // gaitParam.footMidCoords frame
    cnoid::Vector3 rootRpyError = cnoid::rpyFromRot(rootRError); // gaitParam.footMidCoords frame

    for (size_t i = 0; i < 2; i++) {
      stOffsetRootRpy[i] += (this->bodyAttitudeControlGain[i] * rootRpyError[i] - 1.0/this->bodyAttitudeControlTimeConst[i] * stOffsetRootRpy[i]) * dt;
      stOffsetRootRpy[i] = mathutil::clamp(stOffsetRootRpy[i], this->bodyAttitudeControlCompensationLimit[i]);
    }
    stOffsetRootRpy[2] = 0.0;

    o_stOffsetRootRpy.reset(stOffsetRootRpy);
  }

  // stTargetRootPoseを計算
  o_stTargetRootPose.translation() = gaitParam.refRobot->rootLink()->p();
  o_stTargetRootPose.linear() /*generate frame*/= gaitParam.footMidCoords.value().linear() * cnoid::rotFromRpy(o_stOffsetRootRpy.value()/*gaitParam.footMidCoords frame*/) * gaitParam.footMidCoords.value().linear().transpose() * gaitParam.refRobot->rootLink()->R()/*generate frame*/;
  return true;
}

bool Stabilizer::calcZMP(const GaitParam& gaitParam, double dt, bool useActState,
                         cnoid::Vector3& o_tgtZmp, cnoid::Vector3& o_tgtForce, cnoid::Vector3& o_tgtCogAcc) const{
  cnoid::Vector3 cog = useActState ? gaitParam.actCog : gaitParam.genCog;
  cnoid::Vector3 cogVel = useActState ? gaitParam.actCogVel.value() : gaitParam.genCogVel;
  cnoid::Vector3 DCM = cog + cogVel / gaitParam.omega;
  const std::vector<cnoid::Position>& EEPose = useActState ? gaitParam.actEEPose : gaitParam.abcEETargetPose;

  cnoid::Vector3 tgtZmp;
  if(gaitParam.footstepNodesList[0].isSupportPhase[RLEG] || gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    tgtZmp = footguidedcontroller::calcFootGuidedControl(gaitParam.omega,gaitParam.l,DCM,gaitParam.refZmpTraj);
    if(tgtZmp[2] >= gaitParam.actCog[2]) tgtZmp = gaitParam.actCog - cnoid::Vector3(gaitParam.l[0],gaitParam.l[1], 0.0); // 下向きの力は受けられないので
    else{
      // truncate zmp inside polygon. actual robotの関節角度を用いて計算する
      std::vector<cnoid::Vector3> vertices; // generate frame. 支持点の集合
      for(int i=0;i<NUM_LEGS;i++){
        if(!gaitParam.footstepNodesList[0].isSupportPhase[i]) continue;
        for(int j=0;j<gaitParam.legHull[i].size();j++){
          cnoid::Vector3 p = EEPose[i]*gaitParam.legHull[i][j];
          if(p[2] > gaitParam.actCog[2] - 1e-2) p[2] = gaitParam.actCog[2] - 1e-2; // 重心よりも支持点が高いと射影が破綻するので 
          vertices.push_back(p);
        }
      }
      tgtZmp = mathutil::calcInsidePointOfPolygon3D(tgtZmp,vertices,gaitParam.actCog - cnoid::Vector3(gaitParam.l[0],gaitParam.l[1], 0.0));
      // TODO. 角運動量オフセット.
    }
  }else{ // 跳躍期
    tgtZmp = cog - cnoid::Vector3(gaitParam.l[0],gaitParam.l[1], 0.0);
  }
  cnoid::Vector3 tgtCog,tgtCogVel,tgtCogAcc,tgtForce;
  footguidedcontroller::updateState(gaitParam.omega,gaitParam.l,cog,cogVel,tgtZmp,gaitParam.genRobot->mass(),dt,
                                      tgtCog, tgtCogVel, tgtCogAcc, tgtForce);

  // tgtForceにrefEEWrenchのXY成分を足す TODO

  o_tgtZmp = tgtZmp;
  o_tgtForce = tgtForce;
  o_tgtCogAcc = tgtCogAcc;
  return true;
}

bool Stabilizer::calcWrench(const GaitParam& gaitParam, const cnoid::Vector3& tgtZmp/*generate座標系*/, const cnoid::Vector6& tgtWrench/*generate座標系 ロボットが受ける力*/, bool useActState, cnoid::BodyPtr& actRobotTqc, 
                            std::vector<cnoid::Vector6>& o_tgtEEWrench) const{
  cnoid::Vector6 tgtSupWrench = cnoid::Vector6::Zero(); // ルートリンクが支持脚から受ける必要がある外力. generate frame. cog origin.
  // manipulation arm/legの目標反力
  std::vector<cnoid::Vector6> tgtEEManipWrench(gaitParam.eeName.size(), cnoid::Vector6::Zero()); /* 要素数EndEffector数. generate frame. EndEffector origin*/
  for(int i = 0;i<gaitParam.eeName.size();i++){
    if(i < NUM_LEGS && gaitParam.isManualControlMode[i].getGoal() == 0.0) tgtEEManipWrench[i].setZero(); // 支持脚 or 遊脚. 上位からの目標反力は使わない
    else tgtEEManipWrench[i] = gaitParam.refEEWrench[i]; // manipulation arm/leg. 上位からの目標反力を使う
  }
  if(useActState) {
    cnoid::Vector6 tgtSupWrench_o = tgtWrench; // ルートリンクが支持脚から受ける必要がある外力. generate frame. generate frame origin.
    for(int i = 0;i<gaitParam.eeName.size();i++){
      tgtSupWrench_o.head<3>() -= tgtEEManipWrench[i].head<3>();
      tgtSupWrench_o.tail<3>() -= tgtEEManipWrench[i].tail<3>();
      tgtSupWrench_o.tail<3>() -= gaitParam.actEEPose[i].translation().cross(tgtEEManipWrench[i].tail<3>());
    }

    tgtSupWrench.head<3>() = tgtSupWrench_o.head<3>();
    tgtSupWrench.tail<3>() = tgtSupWrench_o.tail<3>();
    tgtSupWrench.tail<3>() += (- actRobotTqc->centerOfMass()).cross(tgtSupWrench_o.head<3>());

  } else {
    return false;
  }

  // トルク制御の目標反力
  std::vector<cnoid::Vector6> tgtEEWrench = tgtEEManipWrench; /* 要素数EndEffector数. generate frame. EndEffector origin*/

  std::vector<int> supportEE;
  if(gaitParam.footstepNodesList[0].isSupportPhase[RLEG] && !gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    supportEE = {RLEG};
  }else if(!gaitParam.footstepNodesList[0].isSupportPhase[RLEG] && gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    supportEE = {LLEG};
  }else if(!gaitParam.footstepNodesList[0].isSupportPhase[RLEG] && !gaitParam.footstepNodesList[0].isSupportPhase[LLEG]){
    // 滞空期
  }else{
    supportEE = {RLEG, LLEG};
  }


  if(supportEE.size()>0){
    /*
      legは、legから受けるwrenchの和がtgtSupWrenchを満たすように.
      各EEFのwrenchを、EEF+copOffset frame/originの6軸表現で考える.

      階層QPのタスクは次の通り
      1. 接触力制約
      2. 和がtgtSupWrench (rot).
      3. 和がtgtSupWrench (trans)
      4. ノルムの2乗和の最小化 (fzは大きくて良い.)

      rotがtransより下の優先度になることにより、擬似的なhip strategyが実現されるようにも思える. しかし、重心位置はFootGuidedControlやmodifyFootStepsで制御できるが、身体の回転を制御する手段は乏しいので、回転しすぎて破綻しないようにrotを優先して満たしたほうが良い
    */

    const int dim = 6 * supportEE.size();
    {
      // 1. 接触力制約
      // 0 <  0  0  1  0  0  0 < 1e10
      // 0 <  1  0 mt  0  0  0 < 1e10
      // 0 < -1  0 mt  0  0  0 < 1e10
      // 0 <  0  1 mt  0  0  0 < 1e10
      // 0 <  0 -1 mt  0  0  0 < 1e10
      // 0 <  0  0  d r1 r2  0 < 1e10 ;; x legHull.size()
      // 0 <  0  0 mr  0  0  1 < 1e10
      // 0 <  0  0 mr  0  0 -1 < 1e10

      this->constraintTask_->A() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->constraintTask_->b() = Eigen::VectorXd::Zero(0);
      this->constraintTask_->wa() = cnoid::VectorX::Ones(0);

      int constraintDim = 0;
      for(int i=0;i<supportEE.size();i++) constraintDim += 7+gaitParam.legHull[supportEE[i]].size();
      this->constraintTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(constraintDim,dim);
      this->constraintTask_->dl() = Eigen::VectorXd::Zero(constraintDim);
      this->constraintTask_->du() = 1e10 * Eigen::VectorXd::Ones(constraintDim);
      this->constraintTask_->wc() = cnoid::VectorX::Ones(constraintDim);
      for(int i=0, idx=0;i<supportEE.size();i++){
        int leg = supportEE[i];
        this->constraintTask_->C().insert(idx,i*6+2) = 1.0; this->constraintTask_->dl()[idx] = 50.0; idx++;
        this->constraintTask_->C().insert(idx,i*6+0) = 1.0; this->constraintTask_->C().insert(idx,i*6+2) = gaitParam.muTrans[leg]; idx++;
        this->constraintTask_->C().insert(idx,i*6+0) = -1.0; this->constraintTask_->C().insert(idx,i*6+2) = gaitParam.muTrans[leg]; idx++;
        this->constraintTask_->C().insert(idx,i*6+1) = 1.0; this->constraintTask_->C().insert(idx,i*6+2) = gaitParam.muTrans[leg]; idx++;
        this->constraintTask_->C().insert(idx,i*6+1) = -1.0; this->constraintTask_->C().insert(idx,i*6+2) = gaitParam.muTrans[leg]; idx++;
        for(int j=0;j<gaitParam.legHull[leg].size();j++){
          cnoid::Vector3 v1 = gaitParam.legHull[leg][j] - gaitParam.copOffset[leg].value(); // EEF+copOffset frame/origin
          cnoid::Vector3 v2 = gaitParam.legHull[leg][(j+1<gaitParam.legHull[leg].size())?j+1:0] - gaitParam.copOffset[leg].value(); // EEF+copOffset frame/origin
          if(v1.head<2>() == v2.head<2>()) continue;
          cnoid::Vector3 r = cnoid::Vector3(v2[1]-v1[1],v1[0]-v2[0],0).normalized();
          double d = r.dot(v1);
          this->constraintTask_->C().insert(idx,i*6+2) = d; this->constraintTask_->C().insert(idx,i*6+3) = -r[1]; this->constraintTask_->C().insert(idx,i*6+4) = r[0]; idx++;
        }
        this->constraintTask_->C().insert(idx,i*6+5) = 1.0; this->constraintTask_->C().insert(idx,2) = gaitParam.muRot[leg]; idx++;
        this->constraintTask_->C().insert(idx,i*6+5) = -1.0; this->constraintTask_->C().insert(idx,2) = gaitParam.muRot[leg]; idx++;
      }

      this->constraintTask_->w() = cnoid::VectorX::Ones(dim) * 1e-6;
      this->constraintTask_->toSolve() = false;
      this->constraintTask_->settings().verbose = 0;
    }

    {
      // 2. 和がtgtSupWrench (rot)
      Eigen::SparseMatrix<double,Eigen::ColMajor> A_colMajor(3,dim); // insertする順番がcolMajorなので
      for(int i=0;i<supportEE.size();i++){
        int leg = supportEE[i];
        cnoid::Position eePose = useActState?gaitParam.actEEPose[leg]:gaitParam.abcEETargetPose[i]; eePose.translation() += eePose.linear() * gaitParam.copOffset[leg].value();
        cnoid::Matrix3 eeR = eePose.linear();
        cnoid::Matrix3 eepCross = mathutil::cross(eePose.translation() - actRobotTqc->centerOfMass()) * eeR;
        for(int j=0;j<3;j++) {
          for(int k=0;k<3;k++) A_colMajor.insert(k,i*6+j) = eepCross(k,j);
        }
        for(int j=0;j<3;j++) {
          for(int k=0;k<3;k++) A_colMajor.insert(k,i*6+3+j) = eeR(k,j);
        }
      }
      this->tgtTorqueTask_->A() = A_colMajor;
      this->tgtTorqueTask_->b() = tgtSupWrench.tail<3>();
      this->tgtTorqueTask_->wa() = cnoid::VectorX::Ones(3);

      this->tgtTorqueTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->tgtTorqueTask_->dl() = Eigen::VectorXd::Zero(0);
      this->tgtTorqueTask_->du() = Eigen::VectorXd::Ones(0);
      this->tgtTorqueTask_->wc() = cnoid::VectorX::Ones(0);

      this->tgtTorqueTask_->w() = cnoid::VectorX::Ones(dim) * 1e-6;
      this->tgtTorqueTask_->toSolve() = true;
      this->tgtTorqueTask_->settings().check_termination = 5; // default 25. 高速化
      this->tgtTorqueTask_->settings().verbose = 0;
    }

    {
      // 3. 和がtgtSupWrench (trans)
      Eigen::SparseMatrix<double,Eigen::ColMajor> A_colMajor(3,dim); // insertする順番がcolMajorなので
      for(int i=0;i<supportEE.size();i++){
        int leg = supportEE[i];
        cnoid::Matrix3 eeR = useActState?gaitParam.actEEPose[leg].linear():gaitParam.abcEETargetPose[i].linear();
        for(int j=0;j<3;j++) {
          for(int k=0;k<3;k++) A_colMajor.insert(k,i*6+j) = eeR(k,j);
        }
      }
      this->tgtForceTask_->A() = A_colMajor;
      this->tgtForceTask_->b() = tgtSupWrench.head<3>();
      this->tgtForceTask_->wa() = cnoid::VectorX::Ones(3);

      this->tgtForceTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->tgtForceTask_->dl() = Eigen::VectorXd::Zero(0);
      this->tgtForceTask_->du() = Eigen::VectorXd::Ones(0);
      this->tgtForceTask_->wc() = cnoid::VectorX::Ones(0);

      this->tgtForceTask_->w() = cnoid::VectorX::Ones(dim) * 1e-6;
      this->tgtForceTask_->toSolve() = true;
      this->tgtForceTask_->settings().check_termination = 5; // default 25. 高速化
      this->tgtForceTask_->settings().verbose = 0;
    }

    {
      // 4. ノルムの2乗和の最小化
      this->normTask_->A() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->normTask_->b() = Eigen::VectorXd::Zero(0);
      this->normTask_->wa() = cnoid::VectorX::Ones(0);

      this->normTask_->C() = Eigen::SparseMatrix<double,Eigen::RowMajor>(0,dim);
      this->normTask_->dl() = Eigen::VectorXd::Zero(0);
      this->normTask_->du() = Eigen::VectorXd::Ones(0);
      this->normTask_->wc() = cnoid::VectorX::Ones(0);

      this->normTask_->w() = cnoid::VectorX::Ones(dim);
      for(int i=0;i<supportEE.size();i++){
        int leg = supportEE[i];
        this->normTask_->w()[i*6+0] = std::pow(1e2, 2.0);
        this->normTask_->w()[i*6+1] = std::pow(1e2, 2.0);
        this->normTask_->w()[i*6+2] = std::pow(1e0, 2.0);
        this->normTask_->w()[i*6+3] = std::pow(1e2, 2.0);
        this->normTask_->w()[i*6+4] = std::pow(1e2, 2.0);
        this->normTask_->w()[i*6+5] = std::pow(1e3, 2.0);
      }

      this->normTask_->toSolve() = true;
      this->normTask_->settings().check_termination = 5; // default 25. 高速化
      this->normTask_->settings().verbose = 0;
    }

    std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks{this->constraintTask_,this->tgtTorqueTask_,this->tgtForceTask_,this->normTask_};
    cnoid::VectorX result; // EEF+copOffset frame/origin
    if(prioritized_qp_base::solve(tasks,
                                   result,
                                   0 // debuglevel
                                   )){
      for(int i=0;i<supportEE.size();i++){
        int leg = supportEE[i];
        cnoid::Vector6 w = result.segment<6>(i*6); // EEF+copOffset frame/origin
        tgtEEWrench[leg].head<3>() += gaitParam.actEEPose[leg].linear() * w.head<3>();
        tgtEEWrench[leg].tail<3>() += gaitParam.actEEPose[leg].linear() * w.tail<3>();
        tgtEEWrench[leg].tail<3>() += (gaitParam.actEEPose[leg].linear() * gaitParam.copOffset[leg].value()).cross(gaitParam.actEEPose[leg].linear() * w.head<3>());
      }
    }
  }

  // エンドエフェクタ力を関節トルクに変換
  for(int i=0;i<gaitParam.eeName.size();i++){
    cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
    cnoid::MatrixXd J = cnoid::MatrixXd::Zero(6,jointPath.numJoints()); // generate frame. endeffector origin
    cnoid::setJacobian<0x3f,0,0,true>(jointPath,actRobotTqc->link(gaitParam.eeParentLink[i]),gaitParam.eeLocalT[i].translation(), // input
                                      J); // output
    cnoid::VectorX tau = - J.transpose() * tgtEEWrench[i];
    for(int j=0;j<jointPath.numJoints();j++){
      jointPath.joint(j)->u() += tau[j];
    }
  }
    
  o_tgtEEWrench = tgtEEWrench;
  return true;
}

bool Stabilizer::calcTorque(double dt, const GaitParam& gaitParam, bool useActState, cnoid::BodyPtr& actRobotTqc, const cnoid::Vector3& targetCogAcc,
                            std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoPGainPercentage, std::vector<cpp_filters::TwoPointInterpolator<double> >& o_stServoDGainPercentage, cnoid::Vector6& o_tgtWrench,
			    Eigen::VectorXd& prev_q, Eigen::VectorXd& prev_dq, std::vector<cnoid::Vector6>& eePoseDiff_prev, std::vector<cnoid::Position>& eeTargetPosed, std::vector<cnoid::Position>& eeTargetPosedd, cnoid::Vector6& prev_rootd) const{

  if(!useActState){
    for(int i=0;i<actRobotTqc->numJoints();i++) actRobotTqc->joint(i)->u() = 0.0;
    for(int i=0;i<actRobotTqc->numJoints();i++){
      o_stServoPGainPercentage[i].interpolate(dt);
      o_stServoDGainPercentage[i].interpolate(dt);
    }
  } else {
    cnoid::VectorX tau_g = cnoid::VectorXd::Zero(actRobotTqc->numJoints()); // 重力
    cnoid::VectorX tau_ee = cnoid::VectorXd::Zero(actRobotTqc->numJoints()); // endEffector

    // 速度・加速度を考慮しない重力補償
    {
      actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
      actRobotTqc->rootLink()->v() = cnoid::Vector3::Zero();
      actRobotTqc->rootLink()->w() = cnoid::Vector3::Zero();
      actRobotTqc->rootLink()->dv() = cnoid::Vector3(0.0,0.0,gaitParam.g);
      actRobotTqc->rootLink()->dw() = cnoid::Vector3::Zero();
      for(int i=0;i<actRobotTqc->numJoints();i++){
	actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	actRobotTqc->joint(i)->dq() = 0.0;
	actRobotTqc->joint(i)->ddq() = 0.0;
      }
      actRobotTqc->calcForwardKinematics(true, true);
      o_tgtWrench = cnoid::calcInverseDynamics(actRobotTqc->rootLink()); // actRobotTqc->joint()->u()に書き込まれる
      for(int i=0;i<actRobotTqc->numJoints();i++){
	tau_g[i] = actRobotTqc->joint(i)->u();
	actRobotTqc->joint(i)->u() = 0.0;
      }
    }

    {
      this->eeComTask_->A() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 * gaitParam.eeName.size() + 3,6 + actRobotTqc->numJoints());
      std::vector<Eigen::Triplet<double> > eeComTripletList_A;
      eeComTripletList_A.reserve(500);//適当
      this->eeComTask_->b() = Eigen::VectorXd::Zero(6 * gaitParam.eeName.size() + 3);
      this->eeComTask_->C() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 + actRobotTqc->numJoints(),6 + actRobotTqc->numJoints());
      std::vector<Eigen::Triplet<double> > eeComTripletList_C;
      // エンドエフェクタ分解加速度制御
      // ee_acc = ee_acc_ref + K (ee_p_ref - ee_p_act) + D (ee_vel_ref - ee_vel_act)
      // ee_vel = J * dq
      // ee_acc = J * ddq + dJ * dq
      // ここで dJ * dq はddq = 0 としたときのee_accだから、q, dqをもとにForwardKinematicsで求められる
      // J * ddq = ee_acc - dJ * dq を満たすddqを求めて、actualのq, dqと合わせてInverseDynamicsでuを求める
      
      {
	std::vector<cnoid::Vector6> ee_acc;
	ee_acc.resize(gaitParam.eeName.size());
	// ee_act
	{
	  for(int i=0;i<gaitParam.eeName.size();i++){
	    // ee_act_ref
	    ee_acc[i].head<3>() = (gaitParam.abcEETargetPose[i].translation() - 2 * eeTargetPosed[i].translation() + eeTargetPosedd[i].translation()) / dt / dt;
	    ee_acc[i].tail<3>() = (cnoid::rpyFromRot(gaitParam.abcEETargetPose[i].linear() * eeTargetPosed[i].linear().transpose()) - cnoid::rpyFromRot(eeTargetPosed[i].linear() * eeTargetPosedd[i].linear().transpose())) / dt / dt; // TODO
	    // K (ee_p_ref - ee_p_act) + D (ee_vel_ref - ee_vel_act)
	    cnoid::Matrix3 eeR = gaitParam.actEEPose[i].linear();
	    cnoid::Vector6 eePoseDiffLocal; // endEfector frame
	    eePoseDiffLocal.head<3>() = eeR.transpose() * (gaitParam.abcEETargetPose[i].translation() - gaitParam.actEEPose[i].translation());
	    cnoid::AngleAxis angleAxis = cnoid::AngleAxis(gaitParam.actEEPose[i].linear().transpose()*gaitParam.abcEETargetPose[i].linear());
	    eePoseDiffLocal.tail<3>() = angleAxis.angle()*angleAxis.axis();
	    cnoid::Vector6 eeVelDiffLocal = (eePoseDiffLocal - eePoseDiff_prev[i]) / dt;
	    cnoid::Vector6 eePoseDiffGainLocal;
	    cnoid::Vector6 eeVelDiffGainLocal;
	    for(int j=0;j<6;j++){
	      eePoseDiffGainLocal[j] = this->ee_K[i][j] * eePoseDiffLocal[j];
	      eeVelDiffGainLocal[j] = this->ee_D[i][j] * eeVelDiffLocal[j];
	    }
	    ee_acc[i].head<3>() += eeR * eePoseDiffGainLocal.head<3>(); // generate frame
	    ee_acc[i].tail<3>() += eeR * eePoseDiffGainLocal.tail<3>(); // generate frame
	    ee_acc[i].head<3>() += eeR * eeVelDiffGainLocal.head<3>(); // generate frame
	    ee_acc[i].tail<3>() += eeR * eeVelDiffGainLocal.tail<3>(); // generate frame
	    if (ee_acc[i].head<3>().norm() > this->ee_dv_limit) ee_acc[i].head<3>() = ee_acc[i].head<3>() / ee_acc[i].head<3>().norm() * this->ee_dv_limit;
	    if (ee_acc[i].tail<3>().norm() > this->ee_dw_limit) ee_acc[i].tail<3>() = ee_acc[i].tail<3>() / ee_acc[i].tail<3>().norm() * this->ee_dw_limit;

	    eePoseDiff_prev[i] = eePoseDiffLocal;
	  }
	} // ee_act

	std::vector<cnoid::Vector6> dJdq;
	dJdq.resize(gaitParam.eeName.size());
	// dJ * dqを求める
	{
	  actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
	  actRobotTqc->rootLink()->v() = gaitParam.actRootVel.value().head<3>();
	  actRobotTqc->rootLink()->w() = gaitParam.actRootVel.value().tail<3>();
	  actRobotTqc->rootLink()->dv() = cnoid::Vector3::Zero();
	  actRobotTqc->rootLink()->dw() = cnoid::Vector3::Zero();
	  for(int i=0;i<actRobotTqc->numJoints();i++){
	    actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	    // dqとしてactualを使うと振動する可能性があるが、referenceを使うと外力による駆動を考慮できない
	    // actRobotTqc->joint(i)->dq() = (gaitParam.genRobot->joint(i)->q() - prev_q[i]) / dt;
	    actRobotTqc->joint(i)->dq() = gaitParam.actRobot->joint(i)->dq();
	    actRobotTqc->joint(i)->ddq() = 0.0;
	  }
	  actRobotTqc->calcForwardKinematics(true, true); // actRobotTqc->link(gaitParam.eeParentLink[i])->dv(), dw()が更新される
	  actRobotTqc->calcCenterOfMass();
	  for(int i=0;i<gaitParam.eeName.size();i++){
	    cnoid::Vector3 arm = actRobotTqc->link(gaitParam.eeParentLink[i])->R() * gaitParam.eeLocalT[i].translation();
	    cnoid::Vector3 ee_dv = actRobotTqc->link(gaitParam.eeParentLink[i])->dv() + actRobotTqc->link(gaitParam.eeParentLink[i])->w().cross(actRobotTqc->link(gaitParam.eeParentLink[i])->w().cross(arm)) + actRobotTqc->link(gaitParam.eeParentLink[i])->dw().cross(arm);
	    cnoid::Vector3 ee_dw = actRobotTqc->link(gaitParam.eeParentLink[i])->dw();
	    dJdq[i].head<3>() = ee_dv;
	    dJdq[i].tail<3>() = ee_dw;
	  }
	} // dJ * dq

	{ // eeComTaskを作る
	  for(int i=0;i<gaitParam.eeName.size();i++){

	    // AはJacobian
	    cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
	    cnoid::MatrixXd J = cnoid::MatrixXd::Zero(6,jointPath.numJoints()); // generate frame. endeffector origin
	    cnoid::setJacobian<0x3f,0,0,true>(jointPath,actRobotTqc->link(gaitParam.eeParentLink[i]),gaitParam.eeLocalT[i].translation(), // input
					      J); // output

	    // 該当する箇所に代入
	    for (int j=0;j<jointPath.numJoints();j++) {
	      for(int k=0;k<6;k++){
		eeComTripletList_A.push_back(Eigen::Triplet<double>(k + 6*i,6+jointPath.joint(j)->jointId(),J(k,j))); // insertすると時間がかかる
	      }
	    }
	  
	    for(int j=0;j<6;j++) eeComTripletList_A.push_back(Eigen::Triplet<double>(j + 6*i,j,1.0)); // insertすると時間がかかる
	    cnoid::Vector3 dp = (actRobotTqc->link(gaitParam.eeParentLink[i])->T() * gaitParam.eeLocalT[i]).translation() - actRobotTqc->rootLink()->p();
	    eeComTripletList_A.push_back(Eigen::Triplet<double>(0 + 6*i,4, dp[2]));
	    eeComTripletList_A.push_back(Eigen::Triplet<double>(0 + 6*i,5,-dp[1]));
	    eeComTripletList_A.push_back(Eigen::Triplet<double>(1 + 6*i,3,-dp[2]));
	    eeComTripletList_A.push_back(Eigen::Triplet<double>(1 + 6*i,5, dp[0]));
	    eeComTripletList_A.push_back(Eigen::Triplet<double>(2 + 6*i,3, dp[1]));
	    eeComTripletList_A.push_back(Eigen::Triplet<double>(2 + 6*i,4,-dp[0]));

	    // bはee_acc - dJ * dq
	    this->eeComTask_->b().block(i*6,0,6,1) = ee_acc[i] - dJdq[i];
	  
	  }

	  for(int i=0;i<6+actRobotTqc->numJoints();i++){
	    eeComTripletList_C.push_back(Eigen::Triplet<double>(i,i,1.0)); // insertすると時間がかかる
	  }
	} // eeTask
      } // ee 分解加速度制御

      // 重心分解加速度制御
      // com_acc = com_acc_ref + K (com_p_ref - com_p_act) + D (com_vel_ref - com_vel_act)
      // com_vel = J * dq
      // com_acc = J * ddq + dJ * dq
      // ここで dJ * dq はddq = 0 としたときのcom_accだから、q, dqをもとにForwardKinematicsで求められる
      // J * ddq = com_acc - dJ * dq を満たすddqを求めて、actualのq, dqと合わせてInverseDynamicsでuを求める
      {
	cnoid::Vector3 com_acc;
	// com_act
	{
	  cnoid::Vector3 actCog = actRobotTqc->centerOfMass();
	  for(int i=0;i<3;i++){
	    com_acc[i] = targetCogAcc[i] + this->com_K[i] * (gaitParam.genCog[i] - actCog[i]) + this->com_D[i] * (gaitParam.genCogVel[i] - gaitParam.actCogVel.value()[i]); // TODO genCogではなくstで出されたcogをつかうこと
	  }
	}

	cnoid::Vector3 dJdq;
	// dJ * dqを求める
	{
	  dJdq = cnoid::Vector3::Zero();
	  actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
	  actRobotTqc->rootLink()->v() = gaitParam.actRootVel.value().head<3>();
	  actRobotTqc->rootLink()->w() = gaitParam.actRootVel.value().tail<3>();
	  actRobotTqc->rootLink()->dv() = cnoid::Vector3::Zero();
	  actRobotTqc->rootLink()->dw() = cnoid::Vector3::Zero();
	  for(int i=0;i<actRobotTqc->numJoints();i++){
	    actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	    // dqとしてactualを使うと振動する可能性があるが、referenceを使うと外力による駆動を考慮できない
	    // actRobotTqc->joint(i)->dq() = (gaitParam.genRobot->joint(i)->q() - prev_q[i]) / dt;
	    actRobotTqc->joint(i)->dq() = gaitParam.actRobot->joint(i)->dq();
	    actRobotTqc->joint(i)->ddq() = 0.0;
	  }
	  actRobotTqc->calcForwardKinematics(true, true);
	  actRobotTqc->calcCenterOfMass();
	  cnoid::Vector6 virtualWrench = cnoid::calcInverseDynamics(actRobotTqc->rootLink());
	  dJdq = virtualWrench.head<3>() / gaitParam.actRobot->mass();
	  dJdq += virtualWrench.tail<3>().cross(actRobotTqc->centerOfMass()) / gaitParam.actRobot->mass();
	}

	{ // comTaskを作る
	  cnoid::MatrixXd CMJ = cnoid::MatrixXd::Zero(3,6 + actRobotTqc->numJoints());
	  cnoid::calcCMJacobian(actRobotTqc,nullptr,CMJ); //仕様でrootは後ろにつくので注意
	  for (int i=0;i<actRobotTqc->numJoints();i++) {
	    for(int j=0;j<3;j++) {
	      eeComTripletList_A.push_back(Eigen::Triplet<double>(6 * gaitParam.eeName.size() + j,6+i,CMJ(j,i)));
	    }
	  }
	  for (int i=0;i<6;i++) {
	    for(int j=0;j<3;j++) {
	      eeComTripletList_A.push_back(Eigen::Triplet<double>(6 * gaitParam.eeName.size() + j,i,CMJ(j,i + actRobotTqc->numJoints())));
	    }
	  }
	  this->eeComTask_->b().tail<3>() = com_acc - dJdq;
	}
      }
      this->eeComTask_->A().setFromTriplets(eeComTripletList_A.begin(), eeComTripletList_A.end());
      this->eeComTask_->C().setFromTriplets(eeComTripletList_C.begin(), eeComTripletList_C.end());
      this->eeComTask_->wa() = Eigen::VectorXd::Ones(6 * gaitParam.eeName.size() + 3);
      this->eeComTask_->dl() = -Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * this->defaultDdqLimit;
      this->eeComTask_->du() = Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * this->defaultDdqLimit;
      for (int i = 0;i<actRobotTqc->numJoints();i++) {
	this->eeComTask_->dl()[6+i] = -this->ddq_limit[i];
	this->eeComTask_->du()[6+i] = this->ddq_limit[i];
      }
      this->eeComTask_->wc() = cnoid::VectorX::Ones(6 + actRobotTqc->numJoints());
      this->eeComTask_->w() = cnoid::VectorX::Ones(6 + actRobotTqc->numJoints()) * 1e-6;
      this->eeComTask_->toSolve() = true;
      this->eeComTask_->settings().check_termination = 15; // default 25. 高速化
      this->eeComTask_->settings().verbose = 0;

      // jointPDTask
      // (角)加速度でのPD制御
      {
	// jointPDTaskを作る
	this->jointPDTask_->A() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 + actRobotTqc->numJoints(),6 + actRobotTqc->numJoints());
	std::vector<Eigen::Triplet<double> > tripletList_A;
	for (int i=0;i<6 + actRobotTqc->numJoints();i++){
	  tripletList_A.push_back(Eigen::Triplet<double>(i,i,1.0));
	}
	this->jointPDTask_->b() = cnoid::VectorXd::Zero(6 + actRobotTqc->numJoints());
	for (int i=0;i<3;i++){
	  this->jointPDTask_->b()[i] = (gaitParam.refRobot->rootLink()->v()[i] - prev_rootd[i]) / dt + this->refAngle_K[i] * (gaitParam.refRobot->rootLink()->p()[i] - gaitParam.actRobot->rootLink()->p()[i]) + this->refAngle_D[i] * (gaitParam.refRobot->rootLink()->v()[i] - gaitParam.actRootVel.value()[i]);
	}
	cnoid::AngleAxis rootAngleAxis = cnoid::AngleAxis(gaitParam.refRobot->rootLink()->R() * gaitParam.actRobot->rootLink()->R().transpose());
	cnoid::Vector3 rootRpyDiff = rootAngleAxis.angle()*rootAngleAxis.axis();
	for (int i=0;i<3;i++){
	  this->jointPDTask_->b()[i+3] = (gaitParam.refRobot->rootLink()->w()[i] - prev_rootd[i+3]) / dt + this->refAngle_K[i+3] * rootRpyDiff[i]  + this->refAngle_D[i+3] * (gaitParam.refRobot->rootLink()->w()[i] - gaitParam.actRootVel.value()[i+3]);
	}
	// refRobotに追従
	for (int i=0;i<actRobotTqc->numJoints();i++){
	  this->jointPDTask_->b()[i+6] = (gaitParam.refRobot->joint(i)->dq() - prev_dq[i]) / dt + this->refAngle_K[i+6] * (gaitParam.refRobot->joint(i)->q() - gaitParam.actRobot->joint(i)->q()) + this->refAngle_D[i+6] * (gaitParam.refRobot->joint(i)->dq() - gaitParam.actRobot->joint(i)->dq());
	}

	this->jointPDTask_->wa() = cnoid::VectorXd::Ones(6 + actRobotTqc->numJoints());
	this->jointPDTask_->C() = Eigen::SparseMatrix<double, Eigen::RowMajor>(6 + actRobotTqc->numJoints(),6 + actRobotTqc->numJoints());
	std::vector<Eigen::Triplet<double> > tripletList_C;
	for (int i=0;i<6 + actRobotTqc->numJoints();i++){
	  tripletList_C.push_back(Eigen::Triplet<double>(i,i,1.0));
	}
	this->jointPDTask_->A().setFromTriplets(tripletList_A.begin(), tripletList_A.end());
        this->jointPDTask_->C().setFromTriplets(tripletList_C.begin(), tripletList_C.end());
	this->jointPDTask_->dl() = - Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * this->defaultDdqLimit;
	this->jointPDTask_->du() = Eigen::VectorXd::Ones(6 + actRobotTqc->numJoints()) * this->defaultDdqLimit;
	for (int i = 0;i<actRobotTqc->numJoints();i++) {
	  this->jointPDTask_->dl()[6+i] = -this->ddq_limit[i];
	  this->jointPDTask_->du()[6+i] = this->ddq_limit[i];
	}
	this->jointPDTask_->wc() = cnoid::VectorX::Ones(6 + actRobotTqc->numJoints());
	this->jointPDTask_->w() = cnoid::VectorXd::Ones(6 + actRobotTqc->numJoints()) * 1e-6;
	this->jointPDTask_->toSolve() = true;
	this->jointPDTask_->settings().check_termination = 15; // default 25. 高速化
	this->jointPDTask_->settings().verbose = 0;
      }      

      // ddqを計算
      std::vector<std::shared_ptr<prioritized_qp_base::Task> > tasks;
      tasks.push_back(this->eeComTask_);
      tasks.push_back(this->jointPDTask_);
	
      cnoid::VectorX result;
      if(!prioritized_qp_base::solve(tasks,
				     result,
				     0 // debuglevel
				     )){
	std::cerr << "fail" << std::endl; // TODO
      }else {
	// ddqを代入
	actRobotTqc->rootLink()->T() = gaitParam.actRobot->rootLink()->T();
	actRobotTqc->rootLink()->v() = gaitParam.actRootVel.value().head<3>();
	actRobotTqc->rootLink()->w() = gaitParam.actRootVel.value().tail<3>();
	actRobotTqc->rootLink()->dv() = result.segment(0,3);//  << result[0], result[1], result[2];
	actRobotTqc->rootLink()->dw() = result.segment(3,3);// << result[3], result[4], result[5];
	for(int i=0;i<actRobotTqc->numJoints();i++){
	  actRobotTqc->joint(i)->q() = gaitParam.actRobot->joint(i)->q();
	  // dqとしてactualを使うと振動する可能性があるが、referenceを使うと外力による駆動を考慮できない
	  // actRobotTqc->joint(i)->dq() = (gaitParam.genRobot->joint(i)->q() - prev_q[i]) / dt;
	  actRobotTqc->joint(i)->dq() = gaitParam.actRobot->joint(i)->dq();
	  actRobotTqc->joint(i)->ddq() = result[6+i];
	}
	actRobotTqc->calcForwardKinematics(true, true);
	actRobotTqc->calcCenterOfMass();
	o_tgtWrench += cnoid::calcInverseDynamics(actRobotTqc->rootLink()); // actRobotTqc->joint()->u()に書き込まれる
	for(int i=0;i<actRobotTqc->numJoints();i++){
	  tau_ee[i] = actRobotTqc->joint(i)->u();
	  actRobotTqc->joint(i)->u() = 0.0;
	}
      }      
    }

    // 最終的な出力トルクを代入
    for(int i=0;i<actRobotTqc->numJoints();i++){
      actRobotTqc->joint(i)->u() = tau_g[i] +  tau_ee[i];
    }

    /*
    std::cerr << "tau_g" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << tau_g[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "tau_lip" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << tau_lip[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "tau_ee" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << tau_ee[i] << " ";
    }
    std::cerr << std::endl;

    std::cerr << "tau" << std::endl; 
    for(int i=0;i<actRobotTqc->numJoints();i++){
      std::cerr << actRobotTqc->joint(i)->u() << " ";
    }
    std::cerr << std::endl;*/
    
    // Gain
    {
      for(int i=0;i<NUM_LEGS;i++){
	cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
	if(gaitParam.isManualControlMode[i].getGoal() == 0.0) { // Manual Control off
	  if(gaitParam.footstepNodesList[0].isSupportPhase[i]){
	    double transitionTime = std::max(this->landing2SupportTransitionTime, dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	    for(int j=0;j<jointPath.numJoints();j++){
	      if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->supportPgain[i][j]) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->supportPgain[i][j], transitionTime);
	      if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->supportDgain[i][j]) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->supportDgain[i][j], transitionTime);
	    }
	  }else if(gaitParam.swingState[i] == GaitParam::DOWN_PHASE) {
	    double transitionTime = std::max(this->swing2LandingTransitionTime, dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	    for(int j=0;j<jointPath.numJoints();j++){
	      if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->landingPgain[i][j]) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->landingPgain[i][j], transitionTime);
	      if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->landingDgain[i][j]) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->landingDgain[i][j], transitionTime);
	    }
	  }else{
	    double transitionTime = std::max(this->support2SwingTransitionTime, dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	    for(int j=0;j<jointPath.numJoints();j++){
	      if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->swingPgain[i][j]) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->swingPgain[i][j], transitionTime);
	      if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != this->swingDgain[i][j]) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(this->swingDgain[i][j], transitionTime);
	    }
	  }
	}else{ // Manual Control on
	  double transitionTime = std::max(gaitParam.isManualControlMode[i].remain_time(), dt*2); // 現状, setGoal(*,dt)以下の時間でgoal指定するとwriteOutPortDataが破綻するのでテンポラリ
	  for(int j=0;j<jointPath.numJoints();j++){
	    if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != 100.0) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(100.0, transitionTime);
	    if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != 100.0) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(100.0, transitionTime);
	  }
	}
      }

      // arm
      for(int i=NUM_LEGS;i<NUM_LEGS+2;i++){
	cnoid::JointPath jointPath(actRobotTqc->rootLink(), actRobotTqc->link(gaitParam.eeParentLink[i]));
	double arm_gain = 0.0;
	for(int j=0;j<jointPath.numJoints();j++){
	  if(o_stServoPGainPercentage[jointPath.joint(j)->jointId()].getGoal() != arm_gain) o_stServoPGainPercentage[jointPath.joint(j)->jointId()].setGoal(arm_gain, 3.0);
	  if(o_stServoDGainPercentage[jointPath.joint(j)->jointId()].getGoal() != arm_gain) o_stServoDGainPercentage[jointPath.joint(j)->jointId()].setGoal(arm_gain, 3.0);
	}
      }
    }
  
    for(int i=0;i<gaitParam.genRobot->numJoints();i++){
      o_stServoPGainPercentage[i].interpolate(dt);
      o_stServoDGainPercentage[i].interpolate(dt);
    }
  } // useActState

  for(int i=0;i<gaitParam.genRobot->numJoints();i++){
    prev_q[i] = gaitParam.genRobot->joint(i)->q();
    prev_dq[i] = gaitParam.refRobot->joint(i)->dq();
  }

  for(int i=0;i<gaitParam.eeName.size();i++){
    eeTargetPosedd[i] = eeTargetPosed[i];
    eeTargetPosed[i] = gaitParam.abcEETargetPose[i];
  }

  prev_rootd.head<3>() = gaitParam.refRobot->rootLink()->v();
  prev_rootd.tail<3>() = gaitParam.refRobot->rootLink()->w();
  
  return true;
}

