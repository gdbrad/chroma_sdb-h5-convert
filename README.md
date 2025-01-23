# chroma_sdb-h5-convert
Converting binary chroma output into hdf5 for perambulators and meson elementals

Reference files for conversion: 
source/chroma/lib/meas/inline/hadron/inline_meson_matelem_colorvec_superb_w.cc

Note: Always set the following 
'''

//! Meson operator, colorstd::vector source and sink with momentum projection
    struct ValMesonElementalOperator_t : public SB::Tensor<2, SB::ComplexD> {
      int type_of_data; /*!< Flag indicating type of data (maybe trivial) */
      ValMesonElementalOperator_t(int n = num_vecs - 1, int type_of_data = COLORVEC_MATELEM_TYPE_DERIV)
	: SB::Tensor<2, SB::ComplexD>("ij", {n, n}, SB::OnHost, SB::Local),
	  type_of_data(type_of_data)
      {
      }
    };
'''
