CXX := g++
CXXFLAGS := -std=c++17 -O3 -march=native -mtune=native -I include
# Detect OpenMP support - handle both g++ and clang++
# On macOS, g++ is actually clang++, so check for libomp
# Try to find libomp via brew
LIBOMP_PATH := $(shell brew --prefix libomp 2>/dev/null)
ifeq ($(LIBOMP_PATH),)
  # Try alternative path
  LIBOMP_PATH := $(shell test -d /opt/homebrew/opt/libomp && echo /opt/homebrew/opt/libomp || echo "")
endif
ifneq ($(LIBOMP_PATH),)
  # libomp found - use it for clang
  # Check if include directory exists, if not try to find omp.h elsewhere
  ifeq ($(shell test -d $(LIBOMP_PATH)/include && echo "yes"),yes)
    LIBOMP_INCLUDE := -I$(LIBOMP_PATH)/include
  else
    # Try to find omp.h in common locations
    OMP_H_PATH := $(shell find $(LIBOMP_PATH) -name "omp.h" 2>/dev/null | head -1)
    ifneq ($(OMP_H_PATH),)
      LIBOMP_INCLUDE := -I$(shell dirname $(OMP_H_PATH))
    else
      # No include found, disable OpenMP
      LIBOMP_PATH :=
    endif
  endif
  ifneq ($(LIBOMP_PATH),)
    OPENMP_FLAG := -Xpreprocessor -fopenmp
    OPENMP_LIB := -L$(LIBOMP_PATH)/lib -lomp
    CXXFLAGS += $(OPENMP_FLAG) $(LIBOMP_INCLUDE)
  endif
endif
# If OpenMP not found above, try standard -fopenmp (works for GCC)
# If that also fails, OpenMP will be disabled (pragmas will be ignored)
ifeq ($(LIBOMP_PATH),)
  # Test if -fopenmp works (GCC) or fails (clang without libomp)
  OPENMP_TEST := $(shell echo 'int main(){return 0;}' | $(CXX) -x c++ -fopenmp - -o /dev/null 2>&1)
  ifeq ($(OPENMP_TEST),)
    # -fopenmp works (GCC)
    OPENMP_FLAG := -fopenmp
    OPENMP_LIB :=
    CXXFLAGS += $(OPENMP_FLAG)
  else
    # -fopenmp doesn't work, OpenMP disabled
    # Pragmas will be ignored due to #ifdef _OPENMP guards
    OPENMP_FLAG :=
    OPENMP_LIB :=
  endif
endif
LDLIBS += $(OPENMP_LIB)
# Optional: uncomment to see vectorization reports
# CXXFLAGS += -ftree-vectorize -fopt-info-vec
# Optional GUI (SFML) support; default off
GUI ?= 0
# Optional: enable slice export from GUI with S key
# Default to enabled for data export functionality
EXPORT_NPY ?= 1
ifeq ($(EXPORT_NPY),1)
  CXXFLAGS += -DEXPORT_NPY
endif
SRCS := src/equations.cpp src/dynamics.cpp src/tornado_sim.cpp src/radiation.cpp src/boundary_layer.cpp src/turbulence.cpp src/numerics.cpp src/simd_utils.cpp \
         src/advection/advection.cpp \
         src/radar.cpp \
         src/radar/base/radar_base.cpp \
         src/radar/factory.cpp \
         src/radar/schemes/reflectivity/reflectivity.cpp \
         src/radar/schemes/velocity/velocity.cpp \
         src/radar/schemes/zdr/zdr.cpp \
         src/microphysics/base/thermodynamics.cpp \
         src/microphysics/factory.cpp \
         src/microphysics/schemes/kessler/kessler.cpp \
         src/microphysics/schemes/lin/lin.cpp \
         src/microphysics/schemes/thompson/thompson.cpp \
         src/microphysics/schemes/milbrandt/milbrandt.cpp \
         src/dynamics/factory.cpp \
         src/dynamics/schemes/supercell/supercell.cpp \
         src/dynamics/schemes/tornado/tornado.cpp \
         src/radiation/base/radiative_transfer.cpp \
         src/radiation/factory.cpp \
         src/radiation/schemes/simple_grey/simple_grey.cpp \
         src/boundary_layer/base/surface_fluxes.cpp \
         src/boundary_layer/factory.cpp \
         src/boundary_layer/schemes/slab/slab.cpp \
         src/boundary_layer/schemes/ysu/ysu.cpp \
         src/boundary_layer/schemes/mynn/mynn.cpp \
         src/turbulence/base/eddy_viscosity.cpp \
         src/turbulence/factory.cpp \
         src/turbulence/schemes/smagorinsky/smagorinsky.cpp \
         src/turbulence/schemes/tke/tke.cpp \
         src/numerics/advection/factory.cpp \
         src/numerics/advection/schemes/tvd/tvd.cpp \
         src/numerics/advection/schemes/weno5/weno5.cpp \
         src/numerics/diffusion/factory.cpp \
         src/numerics/diffusion/schemes/explicit/explicit.cpp \
         src/numerics/diffusion/schemes/implicit/implicit.cpp \
         src/numerics/time_stepping/factory.cpp \
         src/numerics/time_stepping/schemes/rk3/rk3.cpp \
         src/numerics/time_stepping/schemes/rk4/rk4.cpp \
         src/chaos/chaos.cpp \
         src/chaos/base/random_generator.cpp \
         src/chaos/base/perturbation_field.cpp \
         src/chaos/base/correlation_filter.cpp \
         src/chaos/factory.cpp \
         src/chaos/schemes/none/none.cpp \
         src/chaos/schemes/initial_conditions/initial_conditions.cpp \
         src/chaos/schemes/boundary_layer/boundary_layer.cpp \
         src/chaos/schemes/full_stochastic/full_stochastic.cpp \
         src/terrain.cpp \
         src/terrain/base/topography.cpp \
         src/terrain/factory.cpp \
         src/terrain/schemes/bell/bell.cpp \
         src/terrain/schemes/schar/schar.cpp \
         src/terrain/schemes/none.cpp
CPPFLAGS :=
LDLIBS :=
ifeq ($(GUI),1)
  # Detect SFML via pkg-config if available; fallback to Homebrew prefix
  PKG_CONFIG := $(shell command -v pkg-config 2>/dev/null)
  ifeq ($(PKG_CONFIG),)
    SFML_PREFIX ?= $(shell brew --prefix sfml 2>/dev/null)
    SFML_CFLAGS := -I$(SFML_PREFIX)/include
    SFML_LIBS := -L$(SFML_PREFIX)/lib -lsfml-graphics -lsfml-window -lsfml-system
  else
    SFML_CFLAGS := $(shell pkg-config --cflags sfml-graphics)
    SFML_LIBS := $(shell pkg-config --libs sfml-graphics)
  endif
  CPPFLAGS += $(SFML_CFLAGS) -DENABLE_GUI=1
  LDLIBS += $(SFML_LIBS)
  SRCS += src/gui.cpp
endif
BIN := bin/tornado_sim

$(BIN): $(SRCS) | bin
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(SRCS) $(LDLIBS) -o $(BIN)

bin:
	mkdir -p bin

.PHONY: run clean
run: $(BIN)
	./$(BIN)

clean:
	rm -f $(BIN)


