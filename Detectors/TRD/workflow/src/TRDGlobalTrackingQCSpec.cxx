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

/// \file TRDGlobalTrackingQCSpec.cxx
/// \brief Quality control for global tracking (residuals etc)
/// \author Ole Schmidt

#include "TRDWorkflow/TRDGlobalTrackingQCSpec.h"
#include "Framework/Task.h"
#include "Framework/ConfigParamRegistry.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Headers/DataHeader.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "TRDQC/Tracking.h"
#include <TFile.h>
#include <TTree.h>

using namespace o2::framework;
using namespace o2::globaltracking;
using GTrackID = o2::dataformats::GlobalTrackID;

namespace o2
{
namespace trd
{

class TRDGlobalTrackingQC : public Task
{
 public:
  TRDGlobalTrackingQC(std::shared_ptr<DataRequest> dr) : mDataRequest(dr) {}
  ~TRDGlobalTrackingQC() override = default;
  void init(InitContext& ic) final;
  void run(ProcessingContext& pc) final;
  void endOfStream(framework::EndOfStreamContext& ec) final;

 private:
  std::shared_ptr<DataRequest> mDataRequest;
  Tracking mQC;
};

void TRDGlobalTrackingQC::init(InitContext& ic)
{
  //-------- init geometry and field --------//
  o2::base::GeometryManager::loadGeometry();
  o2::base::Propagator::initFieldFromGRP();
  std::unique_ptr<o2::parameters::GRPObject> grp{o2::parameters::GRPObject::loadFrom()};
  mQC.init();
}

void TRDGlobalTrackingQC::run(ProcessingContext& pc)
{

  RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());
  mQC.setInput(recoData);
  mQC.run();
}

void TRDGlobalTrackingQC::endOfStream(EndOfStreamContext& ec)
{
  // for now just dump the results to a file, should become an option to be used in debugging mode
  auto fOut = std::make_unique<TFile>("trdQC.root", "recreate");
  auto tree = std::make_unique<TTree>("qc", "Track based QC for TRD");
  auto vec = mQC.getTrackQC();
  auto vecPtr = &vec;
  tree->Branch("trackQC", &vecPtr);
  tree->Fill();
  tree->Write();
  tree.reset();
  fOut->Close();
}

DataProcessorSpec getTRDGlobalTrackingQCSpec(o2::dataformats::GlobalTrackID::mask_t src)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();

  std::cout << src << std::endl;
  if (GTrackID::includesSource(GTrackID::Source::ITSTPC, src)) {
    LOGF(info, "Found ITS-TPC tracks as input, loading ITS-TPC-TRD");
    src |= GTrackID::getSourcesMask("ITS-TPC-TRD");
  }
  if (GTrackID::includesSource(GTrackID::Source::TPC, src)) {
    LOGF(info, "Found TPC tracks as input, loading TPC-TRD");
    src |= GTrackID::getSourcesMask("TPC-TRD");
  }
  GTrackID::mask_t srcClu = GTrackID::getSourcesMask("TRD"); // we don't need all clusters, only TRD tracklets
  std::cout << src << std::endl;
  dataRequest->requestTracks(src, false);
  dataRequest->requestClusters(srcClu, false);

  return DataProcessorSpec{
    "trd-tracking-qc",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<TRDGlobalTrackingQC>(dataRequest)},
    Options{}};
}

} // namespace trd
} // namespace o2
