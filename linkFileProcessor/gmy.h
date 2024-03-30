#ifndef __GMY_H__
#define __GMY_H__

constexpr uint32_t HemeLbMagicNumber = 0x686c6221;
constexpr uint32_t GmyMagicNumber = 0x676d7904;
constexpr uint32_t GmyNativeMagicNumber = 0x676d7905;
constexpr uint32_t GmyVersionNumber = 4;
constexpr uint32_t GmyNativeVersionNumber = 1;
constexpr size_t PreambleBytes = 32;

constexpr size_t HeaderRecordSize = 12;

struct NonEmptyHeaderRecord {
  uint64_t   blockNumber;  // 8 bytes
  uint64_t   fileOffset;   // 8 bytes
  uint32_t sites;          // 4 bytes
  uint32_t bytes;           // 4 bytes
  uint32_t uncompressedBytes; // 4 bytes
  uint32_t weights;			  // 4 bytes  Total 32 bytes  
  NonEmptyHeaderRecord() : blockNumber(0), fileOffset(0),
						   sites(0), bytes(0), uncompressedBytes(0),
						   weights(0) { }
};

constexpr size_t NonemptyHeaderRecordSize = 32;

struct OutputPreambleInfo {
  uint32_t HemeLBMagic;
  uint32_t GmyNativeMagic;

  uint32_t Version;
  uint32_t BlocksX;

  uint32_t BlocksY;
  uint32_t BlocksZ;

  uint32_t BlockSize;
  uint32_t MaxCompressedBytes;

  uint32_t MaxUncompressedBytes;
  uint32_t HeaderOffset;

  uint64_t NonEmptyBlocks;
  uint64_t DataOffset;
  OutputPreambleInfo() : HemeLBMagic(HemeLbMagicNumber),
                         GmyNativeMagic(GmyNativeMagicNumber),
                         Version(GmyNativeVersionNumber),
                         BlocksX(0), BlocksY(0), BlocksZ(0),
                         MaxCompressedBytes(0), MaxUncompressedBytes(0),
                         HeaderOffset(56),
                         NonEmptyBlocks(0),
                         DataOffset(0) {}

};

struct OutputLink {
    uint32_t linkType;
    uint32_t configID;
    float wallDistance;
};

struct OutputSite {
    uint64_t x, y, z;
    bool hasWallNormal;
    float normalX, normalY, normalZ;
    OutputLink links[26];
};

#endif

