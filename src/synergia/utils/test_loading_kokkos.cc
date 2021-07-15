#include <dlfcn.h>
#include <stdio.h>

#include <Kokkos_Core.hpp>

int main()
{
    void* lib_kokkoscore = dlopen("libkokkoscore.so", RTLD_NOW);
    void* pfun = nullptr;
    if (lib_kokkoscore) {
      // Try calling the Kokkos::initialize() function used in the Python code from a
      // C++ program that is not linked to the library.
      pfun = dlsym(lib_kokkoscore, "_ZN6Kokkos10initializeENS_13InitArgumentsE");
      if (!pfun) {
        printf("Failed to find Kokkos::initialize\n");
      } else {
	using fun_t = void(*)(Kokkos::InitArguments);
	auto call_ctor = reinterpret_cast<fun_t>(pfun);
	(*call_ctor)(Kokkos::InitArguments());
      }
    }
    if (lib_kokkoscore) dlclose(lib_kokkoscore);
    if (lib_kokkoscore == nullptr) return 1;
    if (pfun == nullptr) return 2;
    return 0;
}
