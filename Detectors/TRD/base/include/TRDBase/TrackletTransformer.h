// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef O2_TRD_TRACKLETTRANSFORMER_H
#define O2_TRD_TRACKLETTRANSFORMER_H

#include "TRDBase/Geometry.h"
#include "DataFormatsTRD/Tracklet64.h"
#include "DataFormatsTRD/CalibratedTracklet.h"
#include "DataFormatsTRD/CalVdriftExB.h"

namespace o2
{
namespace trd
{

class TrackletTransformer
{
 public:
  TrackletTransformer() = default;
  ~TrackletTransformer() = default;

  void init();

  float getXCathode() { return mXCathode; }
  float getXAnode() { return mXAnode; }
  float getXDrift() { return mXDrift; }
  float getXtb0() { return mXtb0; }

  void setXCathode(float x) { mXCathode = x; }
  void setXAnode(float x) { mXAnode = x; }
  void setXDrift(float x) { mXDrift = x; }
  void setXtb0(float x) { mXtb0 = x; }

  void setCalVdriftExB(const CalVdriftExB* cal) { mCalVdriftExB = cal; };

  float calculateY(int hcid, int column, int position);

  float calculateZ(int padrow);

  float calculateDy(int hcid, int slope);

  float calibrateX(double x);

  std::array<float, 3> transformL2T(int hcid, std::array<double, 3> spacePoint);

  CalibratedTracklet transformTracklet(Tracklet64 tracklet);

  double getTimebin(int detector, double x);

 private:
  Geometry* mGeo{nullptr};
  const PadPlane* mPadPlane{nullptr};

  float mXCathode;
  float mXAnode;
  float mXDrift;
  float mXtb0;

  const CalVdriftExB* mCalVdriftExB{nullptr};
};

} // namespace trd
} // namespace o2

#endif
