
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "synergia/bunch/bunch.h"

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
        .def_buffer([](HostParticles & p) -> py::buffer_info {
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

        ;


}


