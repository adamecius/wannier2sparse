# FINDINGS — issues spotted while documenting (NOT fixed)

Per the documentation brief, code is **comment-only**: anything that looks like a
bug or smell is recorded here (symptom, `file:line`, why) and left untouched for a
separate change. Nothing below has been fixed.

## Code smells

1. **`using namespace std;` in public headers** — `include/wannier_parser.hpp:26`,
   `include/tbmodel.hpp:29`, `include/hopping_list.hpp:30`.
   *Why:* this is a header-only library; the directive leaks the whole `std`
   namespace into every translation unit that includes these headers (and
   transitively into downstream `find_package` consumers), risking name clashes.
   Fix would be to qualify names (`std::`) and drop the directive — a non-comment
   change, so deferred.

2. **`hopping_storage::operator[](const string&)` aborts on a missing key** —
   `include/hopping_list.hpp:91` and `:108` (`assert(false);` after the linear
   search, then `return this->front();`).
   *Why:* in an `NDEBUG` build the `assert` is compiled out and control falls
   through to `return this->front()`, silently returning the wrong element (or UB
   if the container is empty) instead of failing. A lookup miss should throw or be
   unrepresentable. (Related to the historical map→flat-vector storage migration,
   where tests indexed `hoppings[get_tag(...)]`.)

## Gate status

3. `docs/api-comments` is based directly on `origin/master` and contains only the
   documentation commits, so `bash tools/doc_guard.sh origin/master` is **green**
   (comments in `.hpp/.cpp` + `README.md` + this file + the gate itself). The
   unrelated Plan 11 work lives on its own branch.
