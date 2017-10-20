#include "Kmer.h"

uint32_t Kmer::kSize;

////////////////////////////////////////////////////////////
// kmer
////////////////////////////////////////////////////////////

Kmer::Kmer()
{
  for (int i=0; i < KWORDS; ++i)
    m_data[i] = 0;
}

Kmer::Kmer(const string& seq)
{
  for (int i=0; i < KWORDS; ++i)
    m_data[i] = 0;
  set_sequence(seq);
}

Kmer::Kmer(const Kmer& kmer)
{
  for (int i=0; i < KWORDS; ++i)
    m_data[i] = kmer.m_data[i];
}

bool operator==(const Kmer& lhs, const Kmer& rhs)
{
  for (int i=0; i<KWORDS; i++)
    if (lhs.m_data[i] != rhs.m_data[i])
      return false;
  return true;
}

bool operator<(const Kmer& lhs, const Kmer& rhs)
{
  for (int i=0; i<KWORDS; i++)
    if (lhs.m_data[i] != rhs.m_data[i])
      return (lhs.m_data[i] < rhs.m_data[i]);
  return false;
}
