#include <stdbool.h>
#include "oerrno.h"

bool ofail() {
    return 0 != _ogeterrno();
}
