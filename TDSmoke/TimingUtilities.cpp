#include "TimingUtilities.h"

namespace timingutils 
{

scalar seconds()
{
  struct timeval tv;
  gettimeofday(&tv,NULL);
  scalar time = ((scalar)tv.tv_sec) + ((scalar)tv.tv_usec)*1.0e-6;
  return time;
}

}
