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
//
// file RawEventData.h class  for RAW data format
// Alla.Maevskaya@cern.ch
// with Artur.Furs@cern.ch
//
#ifndef ALICEO2_FIT_RAWEVENTDATA_H_
#define ALICEO2_FIT_RAWEVENTDATA_H_

#include <CommonDataFormat/InteractionRecord.h>
#include <Framework/Logger.h>

namespace o2
{
namespace fit
{

struct EventHeader {
  static constexpr size_t PayloadSize = 10;
  static constexpr size_t PayloadPerGBTword = 10;
  static constexpr size_t MinNelements = 1;
  static constexpr size_t MaxNelements = 1;
  union {
    uint64_t word[2] = {};
    struct {
      uint64_t bc : 12;
      uint64_t orbit : 32;
      uint64_t phase : 3;
      uint64_t errorPhase : 1;
      uint64_t reservedField1 : 16;
      uint64_t reservedField2 : 8;
      uint64_t nGBTWords : 4;
      uint64_t startDescriptor : 4;
      uint64_t reservedField3 : 48;
    };
  };
  InteractionRecord getIntRec() const { return InteractionRecord{(uint16_t)bc, (uint32_t)orbit}; }
  uint16_t getBC() const { return static_cast<uint16_t>(bc); }
  uint32_t getOrbit() const { return static_cast<uint32_t>(orbit); }
  void setIntRec(const InteractionRecord& intRec)
  {
    bc = intRec.bc;
    orbit = intRec.orbit;
  }
  void print() const;
};

struct EventData {
  static constexpr size_t PayloadSize = 5;
  static constexpr size_t PayloadPerGBTword = 10;
  static constexpr size_t MinNelements = 0;
  static constexpr size_t MaxNelements = 12;
  //
  static constexpr int BitFlagPos = 25; // position of first bit flag(numberADC)

  union {
    uint64_t word = {0};
    struct {
      int64_t time : 12;
      int64_t charge : 13;
      uint64_t numberADC : 1, // 25 bit
        isDoubleEvent : 1,
        isTimeInfoNOTvalid : 1,
        isCFDinADCgate : 1,
        isTimeInfoLate : 1,
        isAmpHigh : 1,
        isEventInTVDC : 1,
        isTimeInfoLost : 1,
        reservedField : 3,
        channelID : 4;
    };
  };
  void generateFlags()
  {
    numberADC = std::rand() % 2;
    isDoubleEvent = 0;
    isTimeInfoNOTvalid = 0;
    isCFDinADCgate = 1;
    isTimeInfoLate = 0;
    isAmpHigh = 0;
    isEventInTVDC = 1;
    isTimeInfoLost = 0;
  }
  uint8_t getFlagWord() const
  {
    return uint8_t(word >> BitFlagPos);
  }
  void print() const;
};

struct TCMdata {
  static constexpr size_t PayloadSize = 10;
  static constexpr size_t PayloadPerGBTword = 10;
  static constexpr size_t MinNelements = 1;
  static constexpr size_t MaxNelements = 1;
  uint64_t orA : 1,        // 0 bit (0 byte)
    orC : 1,               // 1 bit
    sCen : 1,              // 2 bit
    cen : 1,               // 3 bit
    vertex : 1,            // 4 bit
    laser : 1,             // 5 bit
    dataIsValid : 1,       // 6 bit
    outputsAreBlocked : 1, // 7 bit
    nChanA : 7,            // 8 bit(1 byte)
    reservedField2 : 1,    // 15 bit
    nChanC : 7,            // 16 bit(2 byte)
    reservedField3 : 1;    // 23 bit
  int64_t amplA : 17,      // 24 bit (3 byte)
    reservedField4 : 1,    // 41 bit
    amplC : 17,            // 42 bit.
    reservedField5 : 1,    // 59 bit.
    // in standard case(without __atribute__((packed)) macros, or packing by using union)
    // here will be empty 4 bits, end next field("timeA") will start from 64 bit.
    timeA : 9,           // 60 bit
    reservedField6 : 1,  // 69 bit
    timeC : 9,           // 70 bit
    reservedField7 : 1,  // 79 bit
    reservedField8 : 48; // 80 bit

  void print() const;
} __attribute__((__packed__));

struct TCMdataExtended {
  static constexpr size_t PayloadSize = 4;
  static constexpr size_t PayloadPerGBTword = 10;
  static constexpr size_t MinNelements = 0;
  static constexpr size_t MaxNelements = 20;
  union {
    uint32_t word[1] = {};
    uint32_t triggerWord;
  };

  void print() const;
};

} // namespace fit
} // namespace o2
#endif
