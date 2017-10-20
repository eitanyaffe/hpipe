#ifndef KMER_H
#define KMER_H

using namespace std;

#include "util.h"

// 32nt per word
#ifndef KWORDS
#define KWORDS 1
#endif

inline const char base2char(uint8_t base) {
  base &= 3;
  switch (base) {
  case 0: return('A');
  case 1: return('C');
  case 2: return('G');
  case 3: return('T');
  default: exit(1);
  }
};

inline const uint8_t char2base(const char ch) {
  switch (ch) {
  case 'A': return(0);
  case 'C': return(1);
  case 'G': return(2);
  case 'T': return(3);
  default: exit(1);
  }
};

////////////////////////////////////////////////////////////
// kmer
////////////////////////////////////////////////////////////

// kmer
struct Kmer
{
  typedef uint64_t word_t;

  ////////////////////////////////////////////////////////////
  // static
  ////////////////////////////////////////////////////////////

  static uint32_t kSize;
  static const uint32_t kBitsPerWord = sizeof(word_t) * 8;
  static const uint32_t kCharsPerWord = kBitsPerWord / 2;

  static const uint32_t kBitsForKmer = KWORDS * kBitsPerWord;
  static const uint32_t kMaxSize = kBitsForKmer / 2;
  static void set_k(int k) {
    massert(k <= (int)kMaxSize, "k=%d exceeds maximal size %d", k, (int)kMaxSize);
    kSize = k;
  };

  ////////////////////////////////////////////////////////////
  // instance
  ////////////////////////////////////////////////////////////

  uint64_t m_data[KWORDS];

  Kmer();
  Kmer(const Kmer& kmer);
  Kmer(const string& seq);

  inline uint8_t get_base(int index) const {
    return (m_data[index / kCharsPerWord] >> ((kCharsPerWord - 1 - index % kCharsPerWord) << 1)) & 3;
  };

  inline void set_base(uint32_t index, uint8_t ch) {
    ch &= 3;
    unsigned offset = (kCharsPerWord - 1 - index % kCharsPerWord) << 1;
    m_data[index / kCharsPerWord] = (m_data[index / kCharsPerWord] & ~(word_t(3) << offset)) | (word_t(ch) << offset);
  }

  inline string get_sequence() const {
    string result = "";
    for (unsigned int i=0; i<kSize; i++)
      result += base2char(get_base(i));
    return result;
  };
  inline void set_sequence(const string& sequence) {
    if (sequence.length() != kSize) {
      cout << "sequence length " << sequence.length() << " must be equal to k length " << kSize << endl;
      exit(1);
    }
    for (unsigned int i=0; i<kSize; i++)
      set_base(i, char2base(sequence[i]));
  };
};

bool operator==(const Kmer& lhs, const Kmer& rhs);
bool operator<(const Kmer& lhs, const Kmer& rhs);

template <class T> class KHash
{
 public:
  T operator()(const Kmer& kmer) const {
    T result = 7;
    for (unsigned int i = 0; i < KWORDS; i++)
      result = result*31 + kmer.m_data[i];
    return result;
  }
};

#endif
