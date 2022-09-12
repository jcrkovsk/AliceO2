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

/// \file PadCalibCCDBBuilder.h
/// \brief Krypton calibration
/// \author Jana Crkovska

#ifndef O2_TRD_KRCALIBRATION_H
#define O2_TRD_KRCALIBRATION_H

#include "DataFormatsTRD/Constants.h"
#include "Rtypes.h"
#include "TH2.h"
#include "TTree.h"
#include <cstdlib>
#include <numeric>
#include <vector>

namespace o2
{
namespace trd
{

class PadCalibCCDBBuilder
{
 public:
  PadCalibCCDBBuilder(){};
  ~PadCalibCCDBBuilder(){};

  void CheckIfIsolatedHotPadCandidate(TH2F* hDet, std::vector<int> coordinates, float upperLimit = 1.5, int areaContainedWithin = 4);
  void CheckIfSmallerCloserToCenter(TH2F* hDet, std::vector<int> coordinates, float allowedDifference);
  std::vector<int> CompareGain(TH2F* hDet, int column, int row, int shiftcolumn, int shiftrow, float allowedDifference);
  float ComputeDetectorAverage(TH2F* hDet);
  float ComputeDistance(std::vector<float> pad1, std::vector<float> pad2);
  TH2F* CreateNormalizedMap(TH2F* hDet, TString sNewName = "");
  void FillInTheGap(TH2F* hDet, int column, int row, float newGain);
  TH2F* FillTheMap(TH2F* hDet, TString sNewName = "", int nbuffer = 3);
  //////////////////////////////////
  TH2F* FillTheMapMirror(TH2F* hDet, TString sNewName = "", int nbuffer = 3);
  TH2F* FillTheMapAverage(TH2F* hDet, TString sNewName = "", int nbuffer = 3);
  //////////////////////////////////
  std::vector<std::vector<int>> FindEmpty(TH2F* hDetectorMap);
  std::vector<std::vector<int>> FindInhomogeneities(TH2F* hDet, float allowedDifference);
  float GetAverageFromNeighbors(TH2F* hDet, int column, int row, int nbuffer = 3);
  TH2F* GetDetectorMap(TTree* tree, int nDet, float mingain = 0, float maxgain = 10'000, TString sDetName = "");
  bool IsHotAreaIsolated(TH2F* hDet, int column, int row, int matrixSize = 1);
  int IsolatedHotPadsContainmentSize(TH2F* hDet, int column, int row);
  void PopulateEmptyNormalizedMap(TH2F* hDet, float valueToSet = -1);
  void RemoveEdges(TH2F* hDet, int nsize = 2);
  void RemoveExtremePads(TH2F* hDet, float upperLimit = 2., float lowerLimit = 0.5);
  void ReplacePadCloserToCenter(TH2F* hDet, int column, int row);
  void ReplaceIsolatedHotPads(TH2F* hDet, int column, int row, int nsize);
  void SetTreeBranches(TTree* tree);
  void SmoothenTheDetector(TH2F* hDet, float allowedDifference = 1000);
  TH2F* TransformMapIntoAbsoluteValues(TH2F* hDet, TString sName = "");

 private:
  float det;
  float col;
  float row;
  float adc;
  float chi;
  float sgm;
  float amp;

  ClassDefNV(PadCalibCCDBBuilder, 1);
};

} // namespace trd
} // namespace o2

#endif // O2_TRD_KRCALIBRATION_H