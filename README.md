# chroma_sdb-h5-convert
Converting binary chroma output into hdf5 for perambulators and meson elementals

Reference files for conversion: 
source/chroma/lib/meas/inline/hadron/inline_meson_matelem_colorvec_superb_w.cc

Make sure to add the conversion binaries to the CMakeLists.txt file in $CHROMA-DIR/chroma-distillation/source/chroma/mainprogs/main/CMakeLists.txt

Make sure the following modules are present in the config file being sourced 
```
module load OpenMPI/4.1.5
module load MPI-settings/CUDA
module load UCX-settings/RC-CUDA
```
Note: Always set the following for the meson elemental conversion ``SB:Local'' needs to be set 
```
//! Meson operator, colorstd::vector source and sink with momentum projection
    struct ValMesonElementalOperator_t : public SB::Tensor<2, SB::ComplexD> {
      int type_of_data; /*!< Flag indicating type of data (maybe trivial) */
      ValMesonElementalOperator_t(int n = num_vecs - 1, int type_of_data = COLORVEC_MATELEM_TYPE_DERIV)
	: SB::Tensor<2, SB::ComplexD>("ij", {n, n}, SB::OnHost, SB::Local),
	  type_of_data(type_of_data)
      {
      }
    };
```

add this version of src/chroma/lib/meas/inline/hadron/inline_create_colorvecs_superb.cc to enable eigenvalue printing. Also add the env variable to shell script when launching eigs 
```
SB_LOG=1
```
