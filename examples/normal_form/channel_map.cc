
#include "synergia/lattice/madx_reader.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/propagator.h"
#include "synergia/utils/cereal_files.h"
#include "synergia/utils/utils.h"

#include <iostream>

void
run()
{
  Logger screen(0, LoggerV::DEBUG);
  Logger log(0, LoggerV::INFO_TURN);

  // read lattice
  MadX_reader reader;
  auto lattice = reader.get_lattice("model", "foborodobo128.madx");

  // tune the lattice
  lattice.set_all_string_attribute("extractor_type", "libff");
  Lattice_simulator::tune_circular_lattice(lattice);

  lattice.print(screen);

  // one turn map
  const int order =
    3; // higher orders use more memory and may not compile on many systems
  auto mapping = Lattice_simulator::get_one_turn_map<order>(lattice);
  // std::cout << mapping.to_json() << "\n";

  // save
  json_save(mapping, "mapping.json");

  // load
  TMapping<Trigon<double, order, 6>> m2;
  json_load(m2, "mapping.json");

  // ref
  auto const& ref = lattice.get_reference_particle();
  auto energy = ref.get_total_energy();
  auto momentum = ref.get_momentum();
  auto mass = ref.get_mass();

  // construct the normal from from loaded mapping
  auto nf = NormalForm<order>(mapping, energy, momentum, mass);

  // save the nf object
  json_save(nf, "nf.json");

  // load
  NormalForm<order> nf2;
  json_load(nf2, "nf.json");

  // print something so it doesnt get optimized away by the compiler
  std::cout << "nf.f[0] = " << nf.get_f()[1].to_json() << "\n";
  std::cout << "nf2.f[0] = " << nf2.get_f()[1].to_json() << "\n";

  for (int comp = 0; comp < 6; ++comp) {
    auto trigon = mapping[comp];
    std::cout << "\n\nComponent " << comp << "\n";

    trigon.each_term([&trigon](int i, auto const& inds, auto term) {
      std::cout << "power = " << inds.size();

      std::array<int, trigon.dim> exp{0};
      for (auto const& x : inds) ++exp[x];

      std::cout << ", exp = [";
      for (int x = 0; x < exp.size() - 1; ++x) std::cout << exp[x] << ", ";
      std::cout << exp.back() << "]";

      std::cout << ", term = " << term << "\n";
    });
  }
}

int
main(int argc, char** argv)
{
  synergia::initialize(argc, argv);

  run();

  synergia::finalize();
  return 0;
}
