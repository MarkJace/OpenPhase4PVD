PYTHON ?= /home/majie/trainFolder/openphaseFolder/op_env/bin/python
OPENPHASE ?= /home/majie/softwares/OpenPhase
CXX ?= g++

# op_env supplies the runtime ABI.  This existing environment is used only to
# locate the header-only pybind11 dependency, which is absent from op_env.
PYBIND11_PYTHON ?= /home/majie/modelFolder/mpfm_code/mpfm_env/bin/python

EXT_SUFFIX := $(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))")
TARGET := openphase_vapor$(EXT_SUFFIX)
PYBIND11_INCLUDES := $(shell $(PYBIND11_PYTHON) -m pybind11 --includes)

.PHONY: all clean

all: $(TARGET)

$(TARGET): openphase_vapor.cpp
	$(CXX) -O3 -DNDEBUG -std=c++17 -fPIC -shared -fopenmp \
		$(PYBIND11_INCLUDES) -I$(OPENPHASE)/include \
		$< -o $@ -L$(OPENPHASE)/lib -lOpenPhase \
		-Wl,-rpath,$(OPENPHASE)/lib -Wl,-z,origin

clean:
	rm -f openphase_vapor*.so
