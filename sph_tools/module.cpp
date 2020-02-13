#include "h_proj.h"

#include <exception>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

using array_t = py::array_t<double, py::array::c_style | py::array::forcecast>;

array_t h_proj_3d(array_t pos, array_t phi, array_t h, size_t L) {
  if (pos.ndim() != 2 || pos.shape(1) != 3) {
    throw std::invalid_argument("position must be an Nx3 array");
  }
  const size_t N = pos.shape(0);

  if (phi.ndim() != 1 || phi.shape(0) != N) {
    throw std::invalid_argument("phi must be a vector of length N");
  }

  if (h.ndim() != 1 || h.shape(0) != pos.shape(0)) {
    throw std::invalid_argument("h must be a vector of length N");
  }

  array_t result({L, L, L});

  const double * pos_in = pos.data();
  const double * phi_in = phi.data();
  const double * h_in = h.data();
  double * out = result.mutable_data();
  {
    py::gil_scoped_release gil;

    double mins[3], maxs[3];
    for (size_t j = 0; j < 3; ++j) {
      mins[j] = maxs[j] = pos_in[j];
    }
    for (size_t i = 1; i < N; ++i) {
      for  (size_t j = 0; j < 3; ++j) {
        mins[j] = std::min(mins[j], pos_in[i*3 + j]);
        maxs[j] = std::max(maxs[j], pos_in[i*3 + j]);
      }
    }

    double scale = 1. / (maxs[0] - mins[0]);
    for (size_t j = 1; j < 3; ++j) {
      scale = std::min(scale, 1./(maxs[j] - mins[j]));
    }
    scale *= L;
    std::unique_ptr<double[]> r_in(new double[N*3]);
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        r_in[i*3 + j] = (pos_in[i*3 + j] - mins[j]) * scale;
      }
    }

    std::unique_ptr<double[]> h_rescaled(new double[N]);
    for (size_t i = 0; i < N; ++i) {
      h_rescaled[i] = h_in[i] * scale;
    }
    h_proj_3d_core(r_in.get(), phi_in, h_rescaled.get(), L, N, out);
  }
  return result;
}

const char * h_proj_3d_DS = R"""(h_proj_3d(pos, phi, h, L)

Deposit an SPH particle quantity onto a regular 3D grid

Parameters
----------
pos : ndarray
    Particle positions (Nx3 array)
phi : ndarray
    Scalar quantity to deposit (N array)
h : ndarray
    SPH particle smoothing lengths (N array)
L : int
    Number of cells along each axis in the resulting grid

Returns
-------
out : ndarray
    The quantity phi deposited onto an LxLxL grid

)""";

PYBIND11_MODULE(sph_tools, m)
{
  using namespace pybind11::literals;

  m.def("h_proj_3d", h_proj_3d, h_proj_3d_DS, "pos"_a, "phi"_a, "h"_a, "L"_a);
}
