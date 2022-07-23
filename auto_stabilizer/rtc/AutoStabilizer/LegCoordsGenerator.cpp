#include "LegCoordsGenerator.h"
#include "MathUtil.h"

#define DEBUG true

void LegCoordsGenerator::calcLegCoords(const GaitParam& gaitParam, double dt,
                                       std::vector<footguidedcontroller::LinearTrajectory<cnoid::Vector3> >& o_refZmpTraj, std::vector<cpp_filters::TwoPointInterpolatorSE3>& o_genCoords, std::vector<GaitParam::FootStepNodes>& o_footstepNodesList, std::vector<cnoid::Position>& o_srcCoords, cpp_filters::TwoPointInterpolatorSE3& o_footMidCoords) const{
  // swing期は、remainTime - supportTime - delayTimeOffset後にdstCoordsに到達するようなantececdent軌道を生成し(genCoords.getGoal()の値)、その軌道にdelayTimeOffset遅れで滑らかに追従するような軌道(genCoords.value()の値)を生成する.
  //   rectangle以外の軌道タイプや跳躍についてはひとまず考えない TODO
  //   srcCoordsとdstCoordsを結ぶ軌道を生成する. srcCoordsの高さ+[0]とdstCoordsの高さ+[1]の高い方(heightとおく)に上げるようなrectangle軌道を生成する
  //     dstに到達するまでの時間が(height-dstCoords)/touchVel以下の場合、dstまで直線で移動する
  //     そうで無い場合、今の位置がheight - eps(ひとまず0)よりも低ければ、heightの高さまで上げてからdstCoordsの上空(XYと回転はdstCoordsと同じ)まで直線で移動する軌道を,残りの自由に使える時間で実行する
  //                     今の位置がheight - eps(ひとまず0)以上であれば、dstCoordsの上空(XYと回転はdstCoordsと同じ)まで直線で移動する軌道を,残りの自由に使える時間で実行する
  //                     両者はチャタリングしそうだが、delayTimeOffset遅れで滑らかに追従するので大丈夫
  // support期は、現FootStepNodesの終了時にdstCoordsに到達するような軌道を線形補間によって生成する.

  // refZmpTrajを更新し進める
  std::vector<footguidedcontroller::LinearTrajectory<cnoid::Vector3> > refZmpTraj = gaitParam.refZmpTraj;
  {
    std::vector<cnoid::Position> stepCoords{gaitParam.genCoords[RLEG].value(),gaitParam.genCoords[LLEG].value()}; // for文中の現在の脚の位置
    cnoid::Vector3 refZmp = refZmpTraj[0].getStart(); // for文中の現在のrefzmp
    refZmpTraj.clear();
    for(int i=0;i<gaitParam.footstepNodesList.size();i++){
      double remainTime = gaitParam.footstepNodesList[i].remainTime;

      // 跳躍についてはひとまず考えない TODO
      // if(gaitParam.footstepNodesList[i].remainTime > gaitParam.footstepNodesList[i].supportTime[RLEG] &&
      //    gaitParam.footstepNodesList[i].remainTime > gaitParam.footstepNodesList[i].supportTime[LLEG])

      if(remainTime > gaitParam.footstepNodesList[i].supportTime[RLEG] && remainTime <= gaitParam.footstepNodesList[i].supportTime[LLEG]){// 右脚がswing. refzmpは左脚の位置
        double swingTime = remainTime - gaitParam.footstepNodesList[i].supportTime[RLEG];
        cnoid::Position llegStartCoords = stepCoords[LLEG];
        cnoid::Position llegGoalCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{stepCoords[LLEG],gaitParam.footstepNodesList[i].dstCoords[LLEG]},std::vector<double>{swingTime, remainTime-swingTime}); // このfootstepNode終了時にdstCoordsに行くように線形補間
        cnoid::Vector3 zmpStartPos = llegStartCoords.translation() + llegStartCoords.linear()*gaitParam.copOffset[LLEG];
        cnoid::Vector3 zmpGoalPos = llegGoalCoords.translation() + llegGoalCoords.linear()*gaitParam.copOffset[LLEG];
        refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(zmpStartPos,zmpGoalPos,swingTime));
        stepCoords[RLEG] = gaitParam.footstepNodesList[i].dstCoords[RLEG];
        stepCoords[LLEG] = llegGoalCoords;
        refZmp = zmpGoalPos;
        remainTime -= swingTime;
      }else if(gaitParam.footstepNodesList[i].remainTime <= gaitParam.footstepNodesList[i].supportTime[RLEG] && gaitParam.footstepNodesList[i].remainTime > gaitParam.footstepNodesList[i].supportTime[LLEG]){ // 左脚がswing. refzmpは右脚の位置
        double swingTime = remainTime - gaitParam.footstepNodesList[i].supportTime[LLEG];
        cnoid::Position rlegStartCoords = stepCoords[RLEG];
        cnoid::Position rlegGoalCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{stepCoords[RLEG],gaitParam.footstepNodesList[i].dstCoords[RLEG]},std::vector<double>{swingTime, remainTime-swingTime}); // このfootstepNode終了時にdstCoordsに行くように線形補間
        cnoid::Vector3 zmpStartPos = rlegStartCoords.translation() + rlegStartCoords.linear()*gaitParam.copOffset[RLEG];
        cnoid::Vector3 zmpGoalPos = rlegGoalCoords.translation() + rlegGoalCoords.linear()*gaitParam.copOffset[RLEG];
        refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(zmpStartPos,zmpGoalPos,swingTime));
        stepCoords[RLEG] = rlegGoalCoords;
        stepCoords[LLEG] = gaitParam.footstepNodesList[i].dstCoords[LLEG];
        refZmp = zmpGoalPos;
        remainTime -= swingTime;
      }
      if(remainTime <= gaitParam.footstepNodesList[i].supportTime[RLEG] && remainTime <= gaitParam.footstepNodesList[i].supportTime[LLEG]){ // 両脚がsupport. refzmpは一つ前の区間と一つ後の区間を線形につなぐ
        cnoid::Position rlegGoalCoords = gaitParam.footstepNodesList[i].dstCoords[RLEG];
        cnoid::Position llegGoalCoords = gaitParam.footstepNodesList[i].dstCoords[LLEG];
        cnoid::Vector3 rlegCOP = rlegGoalCoords.translation() + rlegGoalCoords.linear()*gaitParam.copOffset[RLEG];
        cnoid::Vector3 llegCOP = llegGoalCoords.translation() + llegGoalCoords.linear()*gaitParam.copOffset[LLEG];
        cnoid::Vector3 zmpStartPos = refZmp;
        cnoid::Vector3 zmpGoalPos;
        if(i==gaitParam.footstepNodesList.size()-1 || //末尾. 以降は末尾の状態がずっと続くとして扱うので、refzmpは両足の中心
           gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[RLEG] && gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[LLEG] // 次も両脚支持. refzmpは両足の中心
           ){
          zmpGoalPos = 0.5 * rlegCOP + 0.5 * llegCOP;
        }else if(gaitParam.footstepNodesList[i+1].remainTime > gaitParam.footstepNodesList[i+1].supportTime[RLEG] && gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[LLEG]){ // 次は右脚がswing
          zmpGoalPos = llegCOP;
        }else if(gaitParam.footstepNodesList[i+1].remainTime <= gaitParam.footstepNodesList[i+1].supportTime[RLEG] && gaitParam.footstepNodesList[i+1].remainTime > gaitParam.footstepNodesList[i+1].supportTime[LLEG]){ // 次は左脚がswing
          zmpGoalPos = rlegCOP;
        }else{// 跳躍についてはひとまず考えない TODO
        }
        refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(zmpStartPos,zmpGoalPos,gaitParam.footstepNodesList[i].remainTime));
        stepCoords[RLEG] = rlegGoalCoords;
        stepCoords[LLEG] = llegGoalCoords;
        refZmp = zmpGoalPos;
        remainTime = 0.0;
      }
    }
    // 末尾に1sぶん加える. そうしないと終端条件が厳しすぎる
    refZmpTraj.push_back(footguidedcontroller::LinearTrajectory<cnoid::Vector3>(refZmpTraj.back().getGoal(),refZmpTraj.back().getGoal(), 1.0));
    // dtだけ進める
    if(refZmpTraj[0].getTime() <= dt){
      if(refZmpTraj.size() > 1) refZmpTraj.erase(refZmpTraj.begin());
      else refZmpTraj[0] = footguidedcontroller::LinearTrajectory<cnoid::Vector3>(refZmpTraj[0].getGoal(),refZmpTraj[0].getGoal(),0.0);
    }else{
      refZmpTraj[0] = footguidedcontroller::LinearTrajectory<cnoid::Vector3>(refZmpTraj[0].getStart()+refZmpTraj[0].getSlope()*dt,refZmpTraj[0].getGoal(),refZmpTraj[0].getTime()-dt);
    }
  }

  // genCoordsを進める
  std::vector<cpp_filters::TwoPointInterpolatorSE3> genCoords = gaitParam.genCoords;
  for(int i=0;i<NUM_LEGS;i++){
    if(gaitParam.footstepNodesList[0].remainTime <= gaitParam.footstepNodesList[0].supportTime[i]) { // 支持脚
      cnoid::Position nextCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{genCoords[i].value(),gaitParam.footstepNodesList[0].dstCoords[i]},
                                                           std::vector<double>{std::max(0.0,gaitParam.footstepNodesList[0].remainTime - dt), dt}); // このfootstepNode終了時にdstCoordsに行くように線形補間
      genCoords[i].reset(nextCoords);
    }else{ // 遊脚
      double swingTime = gaitParam.footstepNodesList[0].remainTime - gaitParam.footstepNodesList[0].supportTime[i];
      double height = std::max(gaitParam.srcCoords[i].translation()[2] + gaitParam.footstepNodesList[0].stepHeight[i][0],
                               gaitParam.footstepNodesList[0].dstCoords[i].translation()[2] + gaitParam.footstepNodesList[0].stepHeight[i][1]); // 足上げ高さ. generate frame
      cnoid::Position srcCoords = gaitParam.srcCoords[i];
      cnoid::Position dstCoords = gaitParam.footstepNodesList[0].dstCoords[i];
      cnoid::Position antecedentCoords = genCoords[i].getGoal(); // 今のantecedent軌道の位置
      if(swingTime <= this->delayTimeOffset) {
        genCoords[i].setGoal(dstCoords, swingTime);
        genCoords[i].interpolate(dt);
      }else{
        double antecedentRemainTime = swingTime - this->delayTimeOffset; // antecedent軌道はこの時間後にdstCoordsにつく
        double touchDownTime = (antecedentCoords.translation()[2] - dstCoords.translation()[2]) / this->touchVel; // 地面につくのに要する時間
        cnoid::Position nextCoords;
        if(antecedentRemainTime <= touchDownTime){ // 地面に直線的につきにいく
          nextCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{antecedentCoords,dstCoords},
                                               std::vector<double>{std::max(0.0,antecedentRemainTime - dt), dt}); // このfootstepNode終了時にdstCoordsに行くように線形補間
        }else{
          double preTouchDownRemainTime = antecedentRemainTime - touchDownTime; // dstCoords の上空につくまでの時間
          cnoid::Position preTouchCoords = dstCoords; // dstCoords の上空
          preTouchCoords.translation()[2] = height;
          if(antecedentCoords.translation()[2] >= height){ // dstCoords の上空まで直線で移動する
            nextCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{antecedentCoords,preTouchCoords},
                                                 std::vector<double>{std::max(0.0,preTouchDownRemainTime - dt), dt}); // preTouchDownRemainTime後にpreTouchCoordsに行くように線形補間
          }else{ // heightの高さまで上げてからdstCoords の上空まで直線で移動する. 回転は2つの区間の間一定の速度で動き続ける
            cnoid::Vector3 viaPos = antecedentCoords.translation(); // heightの高さまで上げた位置
            viaPos[2] = height;
            double viaTime = preTouchDownRemainTime * (viaPos - antecedentCoords.translation()).norm() / ((viaPos - antecedentCoords.translation()).norm() + (preTouchCoords.translation() - viaPos).norm()); // viaPosにつくまでの時間
            nextCoords.translation() = mathutil::calcMidPos(std::vector<cnoid::Vector3>{antecedentCoords.translation(),viaPos},
                                                            std::vector<double>{std::max(0.0,viaTime - dt), dt}); // viaTime後にviaPosに行くように線形補間
            nextCoords.linear() =  mathutil::calcMidRot(std::vector<cnoid::Matrix3>{antecedentCoords.linear(),preTouchCoords.linear()},
                                                        std::vector<double>{std::max(0.0,preTouchDownRemainTime - dt), dt}); // preTouchDownRemainTime後にpreTouchCoordsの傾きになるように線形補間
          }
        }
        genCoords[i].setGoal(nextCoords, this->delayTimeOffset);
        genCoords[i].interpolate(dt);
      }
    }
  }

  // footMidCoordsを進める
  cpp_filters::TwoPointInterpolatorSE3 footMidCoords = gaitParam.footMidCoords;
  {
    cnoid::Position rleg = mathutil::orientCoordToAxis(gaitParam.footstepNodesList[0].dstCoords[RLEG], cnoid::Vector3::UnitZ());
    cnoid::Position lleg = mathutil::orientCoordToAxis(gaitParam.footstepNodesList[0].dstCoords[LLEG], cnoid::Vector3::UnitZ());
    cnoid::Position midCoords;
    if(GaitParam::isSupportPhaseEnd(gaitParam.footstepNodesList[0], RLEG) && GaitParam::isSupportPhaseEnd(gaitParam.footstepNodesList[0], LLEG)){ // 両足支持で終わる
      midCoords = mathutil::calcMidCoords(std::vector<cnoid::Position>{rleg, lleg}, std::vector<double>{1.0, 1.0});
    }else if(GaitParam::isSupportPhaseEnd(gaitParam.footstepNodesList[0], RLEG)){ // 右足支持で終わる
      midCoords = rleg; midCoords.translation() -= midCoords.linear() * gaitParam.defaultTranslatePos[RLEG];
    }else{ // 左足支持で終わる
      midCoords = lleg; midCoords.translation() -= midCoords.linear() * gaitParam.defaultTranslatePos[LLEG];
    } // footstepNodesの定義より、両方の足が遊脚の状態で終わることはない
    footMidCoords.setGoal(midCoords, gaitParam.footstepNodesList[0].remainTime);
    footMidCoords.interpolate(dt);
  }

  // footstepNodesListを進める
  std::vector<GaitParam::FootStepNodes> footstepNodesList = gaitParam.footstepNodesList;
  std::vector<cnoid::Position> srcCoords = gaitParam.srcCoords;
  footstepNodesList[0].remainTime = std::max(0.0, footstepNodesList[0].remainTime - dt);
  if(footstepNodesList[0].remainTime <= 0.0 && footstepNodesList.size() > 1){
    for(int i=0;i<NUM_LEGS;i++) srcCoords[i] = gaitParam.genCoords[i].value();
    footstepNodesList.erase(footstepNodesList.begin()); // vectorではなくlistにするべき?
  }

  o_refZmpTraj = refZmpTraj;
  o_genCoords = genCoords;
  o_footstepNodesList = footstepNodesList;
  o_srcCoords = srcCoords;
  o_footMidCoords = footMidCoords;

}

void LegCoordsGenerator::calcCOMCoords(const GaitParam& gaitParam, double dt, double g, double mass, cnoid::Vector3& o_genNextCog, cnoid::Vector3& o_genNextCogVel) const{
  double w = std::sqrt(g/gaitParam.dz); // TODO refforceZ
  cnoid::Vector3 l = cnoid::Vector3::Zero();
  l[2] = gaitParam.dz;
  cnoid::Vector3 genDCM = gaitParam.genCog + w * gaitParam.genCogVel;
  cnoid::Vector3 genZmp = footguidedcontroller::calcFootGuidedControl(w,l,genDCM,gaitParam.refZmpTraj);
  // check zmp in polygon TODO
  cnoid::Vector3 genNextCog,genNextCogVel,genNextForce;
  footguidedcontroller::updateState(w,l,gaitParam.genCog,gaitParam.genCogVel,genZmp,mass,dt,
                                    genNextCog, genNextCogVel, genNextForce);
  o_genNextCog = genNextCog;
  o_genNextCogVel = genNextCogVel;
  std::cerr << gaitParam.refZmpTraj[0].getStart().transpose() << " : " << genNextCog.transpose() << std::endl;
  return;
}
