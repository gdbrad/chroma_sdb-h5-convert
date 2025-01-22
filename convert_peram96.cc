/*! \file
 *  \brief Main program to run all measurement codes.
 */
#include "chroma.h"

#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_matelem.h"
#include <array>
#include <cstddef>
#include <qdp_scalarsite_defs.h>
#include <qdp_stdio.h>
#include <qdp_word.h>
#include <string>
#include "qdp_hdf5.h"

const size_t num_vecs=96;

using namespace Chroma;
extern "C" { 
 void _mcleanup();
}


/*
 * Input 
 */
struct Params_t
{
  multi1d<int>    nrow;
  std::string     inline_measurement_xml;
};

struct Inline_input_t
{
  Params_t        param;
  GroupXML_t      cfg;
  QDP::Seed       rng_seed;
};


void read(XMLReader& xml, const std::string& path, Params_t& p) 
{
  XMLReader paramtop(xml, path);
  read(paramtop, "nrow", p.nrow);

  XMLReader measurements_xml(paramtop, "InlineMeasurements");
  std::ostringstream inline_os;
  measurements_xml.print(inline_os);
  p.inline_measurement_xml = inline_os.str();
  QDPIO::cout << "InlineMeasurements are: " << std::endl;
  QDPIO::cout << p.inline_measurement_xml << std::endl;
}


void read(XMLReader& xml, const std::string& path, Inline_input_t& p) 
{
  try {
    XMLReader paramtop(xml, path);
      
    read(paramtop, "Param", p.param);
    p.cfg = readXMLGroup(paramtop, "Cfg", "cfg_type");

    if (paramtop.count("RNG") > 0)
      read(paramtop, "RNG", p.rng_seed);
    else
      p.rng_seed = 11;     // default seed
  }
  catch( const std::string& e ) 
  {
    QDPIO::cerr << "Error reading XML : " << e << std::endl;
    QDP_abort(1);
  }
}

 

bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();
  foo &= GaugeInitEnv::registerAll();

  return foo;
}

//! Main program to run all measurement codes
/*! \defgroup chromamain Main program to run all measurement codes.
 *  \ingroup main
 *
 * Main program to run all measurement codes.
 */

int main(int argc, char *argv[]) 
{
  // Chroma Init stuff
  Chroma::initialize(&argc, &argv);
  
  START_CODE();
  std::string sdb_file = argv[1];

  QDPIO::cout << "Number of vectors: " << num_vecs << std::endl;

  // Extract the prefix and create the HDF5 file
  std::string h5_file_prefix = sdb_file.substr(0, sdb_file.find_last_of('.'));
  std::string h5_file = h5_file_prefix + ".h5";

  QDPIO::cout << "SDB file: " << sdb_file << std::endl;
  QDPIO::cout << "HDF5 file: " << h5_file << std::endl;

  QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  // We always use the input list
  if (Chroma::getInputFileList().empty())
    {
      Chroma::getInputFileList().push_back(Chroma::getXMLInputFileName());
    }

  // // Sanity check
  // if (!Chroma::getOutputFileList().empty())
  //   {
  //     if (Chroma::getInputFileList().size() != Chroma::getOutputFileList().size())
  // 	{
  // 	  QDPIO::cerr << "Input file list size not equal output file list size" << std::endl;
  // 	  QDP_abort(1);
  // 	}
  //   }


  for(int list_index = 0 ; list_index < Chroma::getInputFileList().size() ; list_index++) 
    {
      Chroma::setXMLInputFileName( Chroma::getInputFileList().at( list_index ) );
      
      XMLReader xml_in;

      // Input parameter structure
      Inline_input_t  input;
      try
	{
	  xml_in.open(Chroma::getXMLInputFileName());
	  read(xml_in, "/chroma", input);
	}
      catch(const std::string& e) 
	{
	  QDPIO::cerr << "CHROMA: Caught Exception reading XML: " << e << std::endl;
	  QDP_abort(1);
	}
      catch(std::exception& e) 
	{
	  QDPIO::cerr << "CHROMA: Caught standard library exception: " << e.what() << std::endl;
	  QDP_abort(1);
	}
      catch(...)
	{
	  QDPIO::cerr << "CHROMA: caught generic exception reading XML" << std::endl;
	  QDP_abort(1);
	}

      XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();

      if (list_index == 0)
	{
	  push(xml_out, "chroma");
	}

      // Write out the input
      write(xml_out, "Input", xml_in);

      if (list_index == 0)
	{
	  Layout::setLattSize(input.param.nrow);
	  Layout::create();
	}
  
      proginfo(xml_out);    // Print out basic program info

      // Initialise the RNG
      QDP::RNG::setrn(input.rng_seed);
      write(xml_out,"RNG", input.rng_seed);

      // Start up the config
      StopWatch swatch;
      swatch.reset();
      multi1d<LatticeColorMatrix> u(Nd);
      XMLReader gauge_file_xml, gauge_xml;

      // Start up the gauge field
      QDPIO::cout << "CHROMA: Attempt to read gauge field" << std::endl;
      swatch.start();
      try 
	{
	  std::istringstream  xml_c(input.cfg.xml);
	  XMLReader  cfgtop(xml_c);
	  QDPIO::cout << "CHROMA: Gauge initialization: cfg_type = " << input.cfg.id << std::endl;

	  Handle< GaugeInit >
	    gaugeInit(TheGaugeInitFactory::Instance().createObject(input.cfg.id,
								   cfgtop,
								   input.cfg.path));
	  (*gaugeInit)(gauge_file_xml, gauge_xml, u);
	}
      catch(std::bad_cast) 
	{
	  QDPIO::cerr << "CHROMA: caught cast error" << std::endl;
	  QDP_abort(1);
	}
      catch(std::bad_alloc) 
	{ 
	  // This might happen on any node, so report it
	  std::cerr << "CHROMA: caught bad memory allocation" << std::endl;
	  QDP_abort(1);
	}
      catch(const std::string& e) 
	{
	  QDPIO::cerr << "CHROMA: Caught Exception: " << e << std::endl;
	  QDP_abort(1);
	}
      catch(std::exception& e) 
	{
	  QDPIO::cerr << "CHROMA: Caught standard library exception: " << e.what() << std::endl;
	  QDP_abort(1);
	}
      catch(...)
	{
	  // This might happen on any node, so report it
	  std::cerr << "CHROMA: caught generic exception during gaugeInit" << std::endl;
	  //QDP_abort(1);
	  throw;
	}
      swatch.stop();

      QDPIO::cout << "CHROMA: Gauge field successfully read: time= " 
		  << swatch.getTimeInSeconds() 
		  << " secs" << std::endl;

      XMLBufferWriter config_xml;
      config_xml << gauge_xml;

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      swatch.reset();
      swatch.start();
  
      // Calculate some gauge invariant observables
      MesPlq(xml_out, "Observables", u);
      swatch.stop();
      QDPIO::cout << "CHROMA: initial plaquette measurement time=" << swatch.getTimeInSeconds() << " secs" << std::endl;


 
      // Get the measurements
      try 
	{
	  swatch.reset();
	  swatch.start();
	  std::istringstream Measurements_is(input.param.inline_measurement_xml);
	  XMLReader MeasXML(Measurements_is);
	  multi1d < Handle< AbsInlineMeasurement > > the_measurements;
	  read(MeasXML, "/InlineMeasurements", the_measurements);

	  QDPIO::cout << "CHROMA: There are " << the_measurements.size() << " measurements " << std::endl;

	  // Reset and set the default gauge field
	  InlineDefaultGaugeField::reset();
	  InlineDefaultGaugeField::set(u, config_xml);

	  // Measure inline observables 
	  push(xml_out, "InlineObservables");
	  xml_out.flush();

	  swatch.stop();
	  QDPIO::cout << "CHROMA: parsing inline measurements time=" << swatch.getTimeInSeconds() << " secs" << std::endl;
	  QDPIO::cout << "CHROMA: Doing " << the_measurements.size() 
		      <<" measurements" << std::endl;
	  swatch.reset();
	  swatch.start();
	  unsigned long cur_update = 0;
	//   for(int m=0; m < the_measurements.size(); m++) 
	//     {
	//       AbsInlineMeasurement& the_meas = *(the_measurements[m]);
	//       if( cur_update % the_meas.getFrequency() == 0 ) 
	// 	{
	// 	  // Caller writes elem rule
	// 	  push(xml_out, "elem");
	// 	  the_meas(cur_update, xml_out);
	// 	  pop(xml_out); 

	// 	  xml_out.flush();
	// 	}
	//     }
	  swatch.stop();

	  QDPIO::cout << "CHROMA: measurements: time= " 
		      << swatch.getTimeInSeconds() 
		      << " secs" << std::endl;


	  pop(xml_out); // pop("InlineObservables");

	  // Reset the default gauge field
	  InlineDefaultGaugeField::reset();
	}
      catch(std::bad_cast) 
	{
	  QDPIO::cerr << "CHROMA: caught cast error" << std::endl;
	  QDP_abort(1);
	}
      catch(std::bad_alloc) 
	{ 
	  // This might happen on any node, so report it
	  std::cerr << "CHROMA: caught bad memory allocation" << std::endl;
	  QDP_abort(1);
	}
      catch(const std::string& e) 
	{
	  QDPIO::cerr << "CHROMA: Caught Exception: " << e << std::endl;
	  QDP_abort(1);
	}
      catch(const char* e) 
	{ 
	  QDPIO::cout << "CHROMA: Caught const char * exception: " << e << std::endl;
	  QDP_abort(1);
	}
      catch(std::exception& e) 
	{
	  QDPIO::cerr << "CHROMA: Caught standard library exception: " << e.what() << std::endl;
	  QDP_abort(1);
	}
      catch(...)
	{
	  // This might happen on any node, so report it
	  std::cerr << "CHROMA: caught generic exception during measurement" << std::endl;
	  std::cerr << "Rethrowing" << std::endl;
	  throw;
	}

      // Last?
      if (list_index == Chroma::getInputFileList().size()-1 )
	{
	  pop(xml_out);
	}

      if (Chroma::getInputFileList().size() > 1)
	{
	  TheNamedObjMap::Instance().erase_all();
	}

    } // list_index

  
  if ( Chroma::getInputFileList().size() > 1 )
    {
      QDPIO::cout << "CHROMA: total number of input files processed = " << Chroma::getInputFileList().size() << std::endl;
    }
  
// HERE BE DRAGONS

  std::cout << "creating obj"<< std::endl;
  BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;
  

  std::cout << "opening file"<< std::endl;
  qdp_db.open(sdb_file, O_RDONLY, 0664);

  std::vector<SerialDBKey<KeyPropElementalOperator_t>> keys;

  std::cout << "reading keys"<< std::endl;
  qdp_db.keys(keys);

  int key_idx = 0;

  SerialDBData<ValPropElementalOperator_t> val;


  QDP::HDF5Writer h5file(h5_file);
  multi2d<double> mat(num_vecs, num_vecs);
  for(auto&& k : keys){
    
    std::cout << "Key number: " << key_idx << std::endl;
    std::cout << "\tt_slice     "<< k.key().t_slice << std::endl;
    std::cout << "\tt_source    "<< k.key().t_source << std::endl;
    std::cout << "\tspin_src   "<< k.key().spin_src << std::endl;
    std::cout << "\tspin_snk   "<< k.key().spin_snk << std::endl;
    std::cout << "\tmass_label "<< k.key().mass_label << std::endl;
    qdp_db.get(k, val);

	std::string path = "/t_source_" + std::to_string(k.key().t_source) + "/t_slice_" + std::to_string(k.key().t_slice) + "/spin_src_" + std::to_string(k.key().spin_src) + "/spin_sink_" + std::to_string(k.key().spin_snk);
	h5file.push(path);

	for(size_t i = 0; i < num_vecs; i++){
		for(size_t j = 0; j < num_vecs; j++){
			mat(i, j) = val.data().mat(i, j).elem().elem().elem().real().elem();
		}
	}
	h5file.write("real", mat);
	for(size_t i = 0; i < num_vecs; i++){
		for(size_t j = 0; j < num_vecs; j++){
			mat(i, j) = val.data().mat(i, j).elem().elem().elem().imag().elem();
		}
	}
	h5file.write("imag", mat);
	

	// for(size_t i = 0; i < 10; i++){
	// 	for(size_t j = 0; j < 10; j++){
	// 		prop_matrix[k.key().t_slice][k.key().spin_src][k.key().spin_snk][i][j] = val.data().mat(i, j);
	// 	}
	// }
    key_idx++;
  }

  h5file.close();
  snoop.stop();
  QDPIO::cout << "CHROMA: total time = "
	      << snoop.getTimeInSeconds() 
	      << " secs" << std::endl;

  QDPIO::cout << "CHROMA: ran successfully" << std::endl;

  END_CODE();

  Chroma::finalize();
  exit(0);
}


