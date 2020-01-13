
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"

#include "synergia/bunch/diagnostics_worker.h"
#include "synergia/bunch/diagnostics_loss.h"

#include "synergia/bunch/diagnostics_py.h"

namespace py = pybind11;
using namespace py::literals;


PYBIND11_MODULE(bunch, m)
{
    py::enum_<ParticleGroup>(m, "ParticleGroup", py::arithmetic())
        .value("regular",   ParticleGroup::regular)
        .value("spectator", ParticleGroup::spectator)
        ;

    py::enum_<LongitudinalBoundary>(m, "LongitudinalBoundary", py::arithmetic())
        .value("open",           LongitudinalBoundary::open)
        .value("periodic",       LongitudinalBoundary::periodic)
        .value("aperture",       LongitudinalBoundary::aperture)
        .value("bucket_barrier", LongitudinalBoundary::bucket_barrier)
        ;

    py::class_<HostParticles>(m, "Particles", py::buffer_protocol())
        .def_buffer([](HostParticles const& p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),       // pointer to buffer
                sizeof(double), // size of one scalar
                py::format_descriptor<double>::format(),
                2,              // num of dimensions
                { p.extent(0), p.extent(1) },    // dimensions
                { p.stride(0) * sizeof(double), 
                  p.stride(1) * sizeof(double) } // strides (in bytes)
            );
        })
        ;

    py::class_<karray1d>(m, "karray1d", py::buffer_protocol())
        .def_buffer([](karray1d const& p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                1,
                { p.extent(0) },
                { p.stride(0) * sizeof(double) }
            );
        })
        ;

    py::class_<karray2d_row>(m, "karray2d_row", py::buffer_protocol())
        .def_buffer([](karray2d_row const& p) -> py::buffer_info {
            return py::buffer_info(
                p.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                2,
                { p.extent(0), p.extent(1) },
                { p.stride(0)*sizeof(double), p.stride(1)*sizeof(double) }
            );
        })
        ;


    // Bunch
    py::class_<Bunch>(m, "Bunch")
        .def( py::init< Reference_particle const&,
                        int, double, Commxx,
                        int, int, int, int >(),
                "Construct a Bunch object.",
                "reference_particle"_a,
                "total_num"_a,
                "real_num"_a,
                "comm"_a = Commxx(),
                "total_spectator_num"_a = 0,
                "bunch_index"_a = 0,
                "bucket_index"_a = 0,
                "array_index"_a = 0 )

        .def( "read_file",
                &Bunch::read_file,
                "Read particle data from file.",
                "filename"_a, 
                "particle_group"_a = ParticleGroup::regular )

        .def( "checkout_particles",
                &Bunch::checkout_particles,
                "Copy the particle array from device memory to host memory.",
                "particle_group"_a = ParticleGroup::regular )

        .def( "checkin_particles",
                &Bunch::checkin_particles,
                "Copy the particle array from device memory to host memory.",
                "particle_group"_a = ParticleGroup::regular )

        .def( "get_host_particles",
                py::overload_cast<ParticleGroup>(&Bunch::get_host_particles),
                "Get host particles.",
                "particle_group"_a = ParticleGroup::regular )

        .def( "add_diagnostics",
                []( Bunch& self,
                    std::shared_ptr<Diagnostics> const& diag, 
                    std::string const& name,
                    std::string const& filename,
                    std::string const& local_dir) {
                        // if the diagnostics is an inherited python type,
                        // reg a ref of the python instance to keep it alive so the
                        // __dict__ object which contains the actual python methods
                        // is available for serializing and calling from C++
                        PyDiagnostics* p = dynamic_cast<PyDiagnostics*>(diag.get());
                        if (p) { p->reg_self(); }

                        self.add_diagnostics(diag, name, filename, local_dir);
                },
                "Add a diagnostics to the bunch object.",
                "diag"_a, "name"_a, "filename"_a, "local_dir"_a = "" )

        .def( "diag_type",
                &Bunch::diag_type,
                "Get the type of the named diagnostics.",
                "name"_a )

        .def( "diag_update",
                &Bunch::diag_update,
                "Performs the update on the named diagnostics.",
                "name"_a )

        .def( "dump",
                &Bunch::dump,
                "Dump." )

        .def( "load",
                &Bunch::load,
                "Load." )
        ;


    // Core_diagnostics
    py::class_<Core_diagnostics>(m, "Core_diagnostics")
        .def_static( "calculate_mean_ka",
                &Core_diagnostics::calculate_mean,
                "Calculate the mean for the bunch.",
                "bunch"_a )

        .def_static( "calculate_std_ka",
                [](Bunch const& bunch, py::buffer b) {
                    py::buffer_info info = b.request();
                    Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                      array((double*)info.ptr, info.shape[0]);
                    return Core_diagnostics::calculate_std(bunch, array);
                },
                "Calculate the standard deviation for the bunch.",
                "bunch"_a, "mean"_a )
        ;

    // Diagnostics_calculator base class
    py::class_<Diagnostics, PyDiagnostics, std::shared_ptr<Diagnostics>>(m, "Diagnostics")
        .def( py::init<std::string const&, bool>(),
                "Construct a Diagnostics object.",
                "type"_a = "py_diag", 
                "serial"_a = true )

        .def( "type",
                &Diagnostics::type,
                "Get the type of the Diagnostics." )

#if 0
        .def( "__reduce__",
                [](py::object const& self) {
                    std::cout << "diag reduce\n";
                    return py::tuple(py::cpp_function([](){
                        return PyDiagnostics("reduced pydiag"); 
                    }), py::tuple() );
                } )
#endif

        .def( py::pickle(
                [](py::object self) { 
                    return py::make_tuple(self.attr("__dict__")); 
                },
                [](py::tuple const& t) { 
                    PyDiagnostics cpp_diag("restored pydiag");
                    auto py_diag = t[0].cast<py::dict>();
                    return std::make_pair(cpp_diag, py_diag);
                } ) )

        ;

    // Diagnostics_dummy
    py::class_<Diagnostics_dummy, Diagnostics, std::shared_ptr<Diagnostics_dummy>>(m, "Diagnostics_dummy")
        .def( py::init<>() )

        ;

}


