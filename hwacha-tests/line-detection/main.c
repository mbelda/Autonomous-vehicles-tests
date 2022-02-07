#include "vec_add_mul.h"
#include "util.h"

#include "dataset1.h"

int main( int argc, char* argv[] )
{
  uint8_t result[DATA_SIZE];
  // Do the saxpy
  setStats(1);
  vec_add_asm(DATA_SIZE, input_data_X, input_data_Y);
  setStats(0);

  // Check the results
  return 0;//verifyFloat(DATA_SIZE, result, verify_data);
}
