CXX := g++
CXXFLAGS := -std=c++17 -O2 -I include
# Optional GUI (SFML) support; default off
GUI ?= 0
# Optional: enable slice export from GUI with S key
# Default to enabled for data export functionality
EXPORT_NPY ?= 1
ifeq ($(EXPORT_NPY),1)
  CXXFLAGS += -DEXPORT_NPY
endif
SRCS := src/equations.cpp src/dynamics.cpp src/tornado_sim.cpp src/radiation.cpp src/boundary_layer.cpp src/turbulence.cpp src/numerics.cpp \
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
         src/chaos/schemes/full_stochastic/full_stochastic.cpp
# Terrain files excluded due to compilation issues - TODO for future contributors
#         src/terrain/base/topography.cpp \
#         src/terrain/factory.cpp \
#         src/terrain/schemes/bell/bell.cpp \
#         src/terrain/schemes/schar/schar.cpp \
#         src/terrain/schemes/none.cpp
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


