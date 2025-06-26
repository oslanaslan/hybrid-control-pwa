#include <cddwrap/cdd.hpp>
#include <setoper.h>
#include <cddmp.h>
#include <cdd.h>

namespace cddwrap {
void global_init() {
  dd_set_global_constants();
}

void global_free() {
  dd_free_global_constants();
}
}  // namespace cddwrap