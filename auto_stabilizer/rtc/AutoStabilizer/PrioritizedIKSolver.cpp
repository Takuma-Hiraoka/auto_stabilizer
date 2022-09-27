#include "PrioritizedIKSolver.h"
#include <prioritized_inverse_kinematics_solver/PrioritizedInverseKinematicsSolver.h>

bool FullbodyIKSolver::solveFullbodyIK(double dt, const GaitParam& gaitParam,
                                       cnoid::BodyPtr& genRobot) const{
  // !jointControllableの関節は指令値をそのまま入れる
  for(size_t i=0;i<genRobot->numJoints();i++){
    if(!gaitParam.jointControllable[i]) genRobot->joint(i)->q() = gaitParam.refRobot->joint(i)->q();
  }

  if(this->jlim_avoid_weight.size() != 6+genRobot->numJoints()) this->jlim_avoid_weight = cnoid::VectorX::Zero(6+genRobot->numJoints());
  cnoid::VectorX dq_weight_all = cnoid::VectorX::Zero(6+genRobot->numJoints());
  for(int i=0;i<6;i++) dq_weight_all[i] = 1.0;
  for(int i=0;i<genRobot->numJoints();i++){
    if(gaitParam.jointControllable[i]) dq_weight_all[6+i] = 1.0;
  }

  std::vector<std::shared_ptr<IK::IKConstraint> > ikConstraint0;
  std::vector<std::shared_ptr<IK::IKConstraint> > ikConstraint1;
  std::vector<std::shared_ptr<IK::IKConstraint> > ikConstraint2;
  std::vector<std::shared_ptr<IK::IKConstraint> > ikConstraint3;

  // collision

  // selfcollision
  this->selfcollisionConstraint.clear();
  for (int i=0;i<gaitParam.selfCollision.size();i++) this->selfcollisionConstraint.push_back(std::make_shared<IK::ClientCollisionConstraint>());
  for(int i=0;i<gaitParam.selfCollision.size();i++){
    this->selfcollisionConstraint[i]->A_link() = genRobot->link(gaitParam.selfCollision[i].link1);
    this->selfcollisionConstraint[i]->A_localp() = gaitParam.selfCollision[i].point1.translation();
    this->selfcollisionConstraint[i]->B_link() = genRobot->link(gaitParam.selfCollision[i].link2);
    this->selfcollisionConstraint[i]->B_localp() = gaitParam.selfCollision[i].point2.translation();
    this->selfcollisionConstraint[i]->tolerance() = 0.01; //TODO
    this->selfcollisionConstraint[i]->maxError() = 10.0*dt;
    this->selfcollisionConstraint[i]->precision() = 0.0; // 強制的にIKをmax loopまで回す
    this->selfcollisionConstraint[i]->weight() = 1.0;
    this->selfcollisionConstraint[i]->velocityDamper() = 0.1 / dt;
    this->selfcollisionConstraint[i]->direction() = gaitParam.selfCollision[i].direction21;
    ikConstraint0.push_back(this->selfcollisionConstraint[i]);
    }

  // envcollision
  // collisionがなくなったとき最後のcollisionが残り続ける
  this->envcollisionConstraint.clear();
  for (int i=0;i<gaitParam.envCollision.size();i++) this->envcollisionConstraint.push_back(std::make_shared<IK::ClientCollisionConstraint>());
  for(int i=0;i<gaitParam.envCollision.size();i++){
    if ((gaitParam.envCollision[i].link1 == "RLEG_JOINT5" )|| (gaitParam.envCollision[i].link1 == "LLEG_JOINT5" ) || (gaitParam.envCollision[i].link1 == "RLEG_JOINT4" )|| (gaitParam.envCollision[i].link1 == "LLEG_JOINT4" )) continue; // 足裏と脛の干渉は無視する
    this->envcollisionConstraint[i]->A_link() = genRobot->link(gaitParam.envCollision[i].link1);
    this->envcollisionConstraint[i]->A_localp() = gaitParam.envCollision[i].point1.translation();
    this->envcollisionConstraint[i]->B_link() = nullptr;
    this->envcollisionConstraint[i]->B_localp() = gaitParam.envCollision[i].point2.translation();
    this->envcollisionConstraint[i]->tolerance() = 0.04; //TODO
    this->envcollisionConstraint[i]->maxError() = 10.0*dt;
    this->envcollisionConstraint[i]->precision() = 0.0; // 強制的にIKをmax loopまで回す
    this->envcollisionConstraint[i]->weight() = 1.0;
    this->envcollisionConstraint[i]->velocityDamper() = 1.0 / dt;
    this->envcollisionConstraint[i]->direction() = gaitParam.envCollision[i].direction21;
    ikConstraint0.push_back(this->envcollisionConstraint[i]);
    }

  // EEF
  for(int i=0;i<gaitParam.eeName.size();i++){
    this->ikEEPositionConstraint[i]->A_link() = genRobot->link(gaitParam.eeParentLink[i]);
    this->ikEEPositionConstraint[i]->A_localpos() = gaitParam.eeLocalT[i];
    this->ikEEPositionConstraint[i]->B_link() = nullptr;
    this->ikEEPositionConstraint[i]->B_localpos() = gaitParam.abcEETargetPose[i];
    this->ikEEPositionConstraint[i]->maxError() << 10.0*dt, 10.0*dt, 10.0*dt, 10.0*dt, 10.0*dt, 10.0*dt;
    this->ikEEPositionConstraint[i]->precision() << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // 強制的にIKをmax loopまで回す
    if(i<NUM_LEGS) this->ikEEPositionConstraint[i]->weight() << 10.0, 10.0, 10.0, 10.0, 10.0, 10.0;
    else this->ikEEPositionConstraint[i]->weight() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    this->ikEEPositionConstraint[i]->eval_link() = nullptr;
    this->ikEEPositionConstraint[i]->eval_localR() = this->ikEEPositionConstraint[i]->B_localpos().linear();
    if(i<NUM_LEGS) ikConstraint1.push_back(this->ikEEPositionConstraint[i]);
    else ikConstraint3.push_back(this->ikEEPositionConstraint[i]);
  }

  // COM
  {
    this->comConstraint->A_robot() = genRobot;
    this->comConstraint->A_localp() = cnoid::Vector3::Zero();
    this->comConstraint->B_robot() = nullptr;
    this->comConstraint->B_localp() = gaitParam.genCog + gaitParam.sbpOffset;
    this->comConstraint->maxError() << 10.0*dt, 10.0*dt, 10.0*dt;
    this->comConstraint->precision() << 0.0, 0.0, 0.0; // 強制的にIKをmax loopまで回す
    this->comConstraint->weight() << 3.0, 3.0, 1.0;
    this->comConstraint->eval_R() = cnoid::Matrix3::Identity();
    ikConstraint1.push_back(this->comConstraint);
  }

  // Angular Momentum
  {
    this->angularMomentumConstraint->robot() = genRobot;
    this->angularMomentumConstraint->targetAngularMomentum() = cnoid::Vector3::Zero(); // TODO
    this->angularMomentumConstraint->maxError() << 1.0*dt, 1.0*dt, 1.0*dt;
    this->angularMomentumConstraint->precision() << 0.0, 0.0, 0.0; // 強制的にIKをmax loopまで回す
    this->angularMomentumConstraint->weight() << 1e-4, 1e-4, 0.0; // TODO
    this->angularMomentumConstraint->dt() = dt;
    this->comConstraint->eval_R() = cnoid::Matrix3::Identity();
    ikConstraint3.push_back(this->angularMomentumConstraint);
  }

  // root
  {
    this->rootPositionConstraint->A_link() = genRobot->rootLink();
    this->rootPositionConstraint->A_localpos() = cnoid::Position::Identity();
    this->rootPositionConstraint->B_link() = nullptr;
    this->rootPositionConstraint->B_localpos() = gaitParam.stTargetRootPose;
    this->rootPositionConstraint->maxError() << 10.0*dt, 10.0*dt, 10.0*dt, 10.0*dt, 10.0*dt, 10.0*dt;
    this->rootPositionConstraint->precision() << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // 強制的にIKをmax loopまで回す
    this->rootPositionConstraint->weight() << 0.0, 0.0, 0.0, 3.0, 3.0, 3.0; // 角運動量を利用するときは重みを小さく. 通常時、胴の質量・イナーシャやマスパラ誤差の大きさや、胴を大きく動かすための出力不足などによって、二足動歩行では胴の傾きの自由度を使わない方がよい
    //this->rootPositionConstraint->weight() << 0.0, 0.0, 0.0, 3e-1, 3e-1, 3e-1;
    this->rootPositionConstraint->eval_link() = nullptr;
    this->rootPositionConstraint->eval_localR() = cnoid::Matrix3::Identity();
    ikConstraint2.push_back(this->rootPositionConstraint);
  }

  // reference angle
  {
    for(size_t i=0;i<genRobot->numJoints();i++){
      if(!gaitParam.jointControllable[i]) continue;
      this->refJointAngleConstraint[i]->joint() = genRobot->joint(i);
      this->refJointAngleConstraint[i]->maxError() = 10.0 * dt; // 高優先度のmaxError以下にしないと優先度逆転するおそれ
      this->refJointAngleConstraint[i]->weight() = 1e-1; // 小さい値すぎると、qp終了判定のtoleranceによって無視されてしまう
      this->refJointAngleConstraint[i]->targetq() = gaitParam.refRobot->joint(i)->q();
      this->refJointAngleConstraint[i]->precision() = 0.0; // 強制的にIKをmax loopまで回す
      ikConstraint3.push_back(this->refJointAngleConstraint[i]);
    }
  }

  // joint limit
  {
    for(int i=0;i<genRobot->numJoints();i++){
      if(!gaitParam.jointControllable[i]) continue;
      double u = gaitParam.refRobot->joint(i)->q_upper();
      double l = gaitParam.refRobot->joint(i)->q_lower();
      for(int j=0;j<gaitParam.jointLimitTables[i].size();j++){
	u = std::min(u,gaitParam.jointLimitTables[i][j]->getUlimit());
	l = std::max(l,gaitParam.jointLimitTables[i][j]->getLlimit());
      }
      genRobot->joint(i)->setJointRange(l,u);
    }
        
    for(size_t i=0;i<genRobot->numJoints();i++){
      if(!gaitParam.jointControllable[i]) continue;
      this->jointLimitConstraint[i]->joint() = genRobot->joint(i);
      this->jointLimitConstraint[i]->maxError() = 1.0 * dt; 
      this->jointLimitConstraint[i]->weight() = 1.0; 
      this->jointLimitConstraint[i]->precision() = 0.0; // 強制的にIKをmax loopまで回す
      ikConstraint0.push_back(this->jointLimitConstraint[i]);
    }    
  }

  // joint velocity
  {    
    for(size_t i=0;i<genRobot->numJoints();i++){
      if(!gaitParam.jointControllable[i]) continue;
      this->jointVelocityConstraint[i]->joint() = genRobot->joint(i);
      this->jointVelocityConstraint[i]->dt() = dt;
      this->jointVelocityConstraint[i]->maxError() = 0.1 * dt; 
      this->jointVelocityConstraint[i]->weight() = 1.0; 
      this->jointVelocityConstraint[i]->precision() = 0.0; // 強制的にIKをmax loopまで回す
      ikConstraint0.push_back(this->jointVelocityConstraint[i]);
    }    
  }
  
  // 特異点近傍で振動するようなことは起こりにくいが、歩行動作中の一瞬だけIKがときにくい姿勢があってすぐに解ける姿勢に戻るといった場合に、その一瞬の間だけIKを解くために頑張って姿勢が大きく変化するので、危険.
  //  この現象を防ぐには、未来の情報を含んだIKを作るか、歩行動作中にIKが解きづらい姿勢を経由しないように着地位置等をリミットするか. 後者を採用
  //  歩行動作ではないゆっくりとした動作であれば、この現象が発生しても問題ない

  std::vector<cnoid::LinkPtr> variables;
  variables.push_back(genRobot->rootLink());
  for(size_t i=0;i<genRobot->numJoints();i++){
    variables.push_back(genRobot->joint(i));
  }
  std::vector<std::vector<std::shared_ptr<IK::IKConstraint> > > constraints{ikConstraint0,ikConstraint1,ikConstraint2,ikConstraint3};
  for(size_t i=0;i<constraints.size();i++){
    for(size_t j=0;j<constraints[i].size();j++){
      constraints[i][j]->debuglevel() = 0;//debug
    }
  }

  prioritized_inverse_kinematics_solver::solveIKLoop(variables,
						     constraints,
						     this->tasks,
						     1,//loop
						     1e-6, // wn
						     0, //debug
						     dt
						     );

  return true;
}
