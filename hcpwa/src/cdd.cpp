#include <cdd.hpp>

void GlobalInit() {
  dd_set_global_constants();
}

void GlobalFree() {
  dd_free_global_constants();
}
