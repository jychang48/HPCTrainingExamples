HIPC=hipcc

# All supported radius values
ALL_RADII := 1 2 3 4 6 8
# All supported vec values (VEC_EXP: 0=vec1, 1=vec2, 2=vec4)
ALL_VECS := 0 1 2

# Track if filters were explicitly specified
RAD_SPECIFIED := $(if $(rad),yes,no)
VEC_SPECIFIED := $(if $(vec),yes,no)

# If rad is specified, use only that; otherwise use all
ifdef rad
  RADII := $(rad)
else
  RADII := $(ALL_RADII)
endif

# If vec is specified, use only that; otherwise use all
ifdef vec
  VECS := $(vec)
else
  VECS := $(ALL_VECS)
endif

# Directories
SRC_DIR=src
INC_DIR=include
BUILD_DIR=build
BIN_DIR=bin

# GPU Architecture: can be specified via command line (e.g., make ARCH=gfx1100)
# If not specified, auto-detect from available GPU using rocminfo
# Falls back to gfx90a (MI250X) if detection fails
ifndef ARCH
  DETECTED_ARCH := $(shell rocminfo 2>/dev/null | grep -o -m1 'gfx[0-9a-z]*')
  ifeq ($(DETECTED_ARCH),)
    ARCH := gfx90a
    $(info Auto-detection failed, using default: $(ARCH))
  else
    ARCH := $(DETECTED_ARCH)
    $(info Auto-detected GPU architecture: $(ARCH))
  endif
else
  $(info Using specified GPU architecture: $(ARCH))
endif

# Optional: save temp files (assembly, IR) for debugging
# Usage: make SAVE_TEMPS=1
ifdef SAVE_TEMPS
  TEMPS_FLAG=--save-temps
else
  TEMPS_FLAG=
endif

# Generate all target configurations (binaries go in bin/)
TARGETS := $(foreach r,$(RADII),$(foreach v,$(VECS),$(BIN_DIR)/baseline_R$(r)_vec$(v)_$(ARCH).x))

# Function to get config for a specific radius/vec combination (for binary naming)
define get_config
R$1_vec$2_$(ARCH)
endef

# Function to get build subdir for a specific radius/vec combination
# Structure: build/rad{radius}/vec{vec}/{arch}
define get_build_subdir
$(BUILD_DIR)/rad$1/vec$2/$(ARCH)
endef

# Function to get object files for a specific radius/vec combination
define get_objs
$(call get_build_subdir,$1,$2)/main.o $(call get_build_subdir,$1,$2)/baseline_kernel.o $(call get_build_subdir,$1,$2)/optimized_v1_kernel.o
endef

# Function to get HIPFLAGS for a specific radius/vec combination
define get_hipflags
-O3 -g --offload-arch=$(ARCH) -std=c++17 -I$(CURDIR)/$(INC_DIR) -DRADIUS=$1 -DVEC_EXP=$2 $(TEMPS_FLAG)
endef

# Default target: build all configurations
all: $(TARGETS)

# Rule to build each target
define build_target_rule
$(BIN_DIR)/baseline_R$1_vec$2_$(ARCH).x: $(call get_objs,$1,$2) | $(BIN_DIR)
	$$(HIPC) $$(call get_hipflags,$1,$2) $$^ -o $$@
endef

# Generate build rules for each radius/vec combination
$(foreach r,$(RADII),$(foreach v,$(VECS),$(eval $(call build_target_rule,$r,$v))))

# Rule to compile main.cpp for a specific radius/vec combination
define compile_main_rule
$(call get_build_subdir,$1,$2)/main.o: $(SRC_DIR)/main.cpp | $(call get_build_subdir,$1,$2)
	cd $$(call get_build_subdir,$1,$2) && $$(HIPC) $$(call get_hipflags,$1,$2) -c $$(CURDIR)/$$< -o main.o
endef

# Generate compile rules for each radius/vec combination
$(foreach r,$(RADII),$(foreach v,$(VECS),$(eval $(call compile_main_rule,$r,$v))))

# Rule to compile baseline_kernel.hip for a specific radius/vec combination
define compile_kernel_rule
$(call get_build_subdir,$1,$2)/baseline_kernel.o: $(SRC_DIR)/baseline_kernel.hip | $(call get_build_subdir,$1,$2)
	cd $$(call get_build_subdir,$1,$2) && $$(HIPC) $$(call get_hipflags,$1,$2) -c $$(CURDIR)/$$< -o baseline_kernel.o
endef

# Generate compile rules for each radius/vec combination
$(foreach r,$(RADII),$(foreach v,$(VECS),$(eval $(call compile_kernel_rule,$r,$v))))

# Rule to compile optimized_v1_kernel.hip for a specific radius/vec combination
define compile_optimized_v1_kernel_rule
$(call get_build_subdir,$1,$2)/optimized_v1_kernel.o: $(SRC_DIR)/optimized_v1_kernel.hip | $(call get_build_subdir,$1,$2)
	cd $$(call get_build_subdir,$1,$2) && $$(HIPC) $$(call get_hipflags,$1,$2) -c $$(CURDIR)/$$< -o optimized_v1_kernel.o
endef

# Generate compile rules for each radius/vec combination
$(foreach r,$(RADII),$(foreach v,$(VECS),$(eval $(call compile_optimized_v1_kernel_rule,$r,$v))))

# Rule to create build subdirectory
define mkdir_rule
$(call get_build_subdir,$1,$2):
	mkdir -p $$@
endef

# Generate mkdir rules for each radius/vec combination
$(foreach r,$(RADII),$(foreach v,$(VECS),$(eval $(call mkdir_rule,$r,$v))))

# Create bin directory
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Clean: if rad/vec specified, clean only those; otherwise clean everything
clean:
ifeq ($(RAD_SPECIFIED)$(VEC_SPECIFIED),nono)
	@echo "Cleaning everything"
	rm -rf $(BUILD_DIR) $(BIN_DIR) *.x
else
	@echo "Cleaning configurations: rad=$(RADII), vec=$(VECS)"
	$(foreach r,$(RADII),$(foreach v,$(VECS),rm -rf $(call get_build_subdir,$r,$v) $(BIN_DIR)/baseline_R$(r)_vec$(v)_$(ARCH).x;))
endif

.PHONY: all clean
