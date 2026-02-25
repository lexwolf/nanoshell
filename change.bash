git grep -l '^[[:space:]]*#include[[:space:]]*"headers/cup\.H"[[:space:]]*$' -- '*.c' '*.cc' '*.cpp' '*.cxx' '*.h' '*.hh' '*.hpp' '*.H' \
| while IFS= read -r f; do
    perl -0777 -i -pe '
      s{
        ^([ \t]*)\#include[ \t]*"headers/cup\.H"[ \t]*\r?\n
      }{$1#define CUP_BACKEND_QUASI_STATIC\n$1#include "headers/cup.H"\n}mgx
    ' "$f"
  done