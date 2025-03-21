### Examination of `chroma` building block data structures: Towards distilling distillation..
There is a sequential dependency chain within `chroma`:
sequential dependency:

  1. Gauge generation (hmc)
  2. Eigenvector generation (chroma single prec, laplace_eigs)
  3.  Distillation objects (baryon, meson, props, genprops): chroma double prec, harom
  4.  Corr. func. generation: redstar, colorvec, harom
  5.  Analysis


There is no cpp source file that specifies how distillation is carried out in `chroma`. Rather, the means by which distilled objects are index and stored is not obvious, this is hidden within the `harom` external library. The nature of the problem is thus: we must "reverse engineer" the process to characterize the data structures for distilled propagators and the associated perambulators on each time slice. 

This is a recognized issue, according to https://github.com/eromero-vlc/chroma-challenges:

> Minimal support to introspect eigenvectors and distillation objects and check properties; which makes debugging and finding mistakes on the input files cumbersome



We start with the input file `prop_and_pearm.ini.xml`:
```
<NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <colorvec_files><elem>./colorvec.mod</elem></colorvec_files>
        <prop_op_file>./peram.sdb</prop_op_file>
      </NamedObject>
```
However, this fails when passing to the main `chroma` executable. 

We also have the test input file `src/chroma/tests/chroma/hadron/colorvec_matelem/colorvec_matelem_test.ini.xml`. Let's pass this to the main chroma executable. As we are only interested in mesons for the time being, we can safely remove the baryon measurements in the xml file. 



We need to generate the following .sdb files to then read with `read_sdb` executable, in order of attainment:
- `meson.sdb`
- `prob.sdb`
- `peram.sdb`
I think at some  point we should have a `colorvec.sdb`? 

We cant? use `meson.ini.xml` because it is wrapped in the `harom` library which is closed source. 


1. First order of business is to generate a `colorvec.mod` file using the input file above 
2. Next we use in file `src/chroma/tests/chroma/hadron/distillation/meson.ini.xml` to generate a `meson.sdb` file? We must pass `colorvec.mod` to this ?





In turn, we will use the `src/chroma/lib/util/ferm/key_prop_matelem.cc` to transform the above objects into `KeyPropElementalOperator` and `ValPropElementalOperator`

The information we need to extract is 
$$ V^{\dagger}M^{-1}V \rightarrow \tau $$ 
where $\tau$ is the perambulator matrix on a single time slice? 

Propagators transform with tensor product structure 
$$\text{Lattice} \otimes \text{Matrix(Nc)} \otimes \text{Matrix(Ns)} \otimes \text{Complex}$$

We can work out these dimensions for ourselves; A distilled propagator stored on disk has dimensions 
$ 2 * 8 *2 *4 *4 *10 * 10 * 16$ 
with the dictionary 
$$\text{} \text{complex} * \text{snk} * \text{src} * N \times N_{\sigma} * \text{tslice}$$

The information for how an elemental operator is stored in chroma resides in 
`src/chroma/lib/util/ferm/key_prop_matelem.h`. We see that the structure is like so 
```//! Prop operator
  struct KeyPropElementalOperator_t
  {
    int                t_slice;       /*!< Propagator time slice */
    int                t_source;      /*!< Source time slice */
    int                spin_src;      /*!< Source spin index */
    int                spin_snk;      /*!< Sink spin index */
    std::string        mass_label;    /*!< A mass label */
  };
  ```


At this point, we need to actually perform contractions to obtain the correlator 
$$C_M^{(2)}(t',t) = Tr[\Phi^B(t')\tau(t',t)\Phi^A(t)\tau(t,t')]$$ 
where 
$$\Phi^A_{\alpha\beta}(t) = V^{\dagger}(t) [\Gamma^A(t)]_{\alpha\beta} V(t) \equiv V^{\dagger}(t)\mathcal{D}^A(t)V(t)S^A_{\alpha\beta}$$ 
and 
$$\tau_{\alpha\beta}(t',t) = V^{\dagger}(t')M_{\alpha\beta}^{-1}(t',t)V(t)$$ 
is the perambulator, defined by the lattice representation of the Dirac operator, $M$. 
See [https://arxiv.org/abs/0905.2160v1]. 

## Demystifying what the `harom` library is 

The `harom` library used? to build with chroma but the code is closed source; In the XML files used to generate `MESON_MATELEM_COLORVEC` in `meson_test.xml` , the entire input file is braced in  a call to `harom`. 

We need to investigate the data structure of `meson.sdb`

The corresponding key in `std::map` 






## problem statement 

Take components of `chroma` output :
- quark lines 
- perambulators 
- random vectors 
- eigenvector space 
And implement the tensor contractions. Then we project out the correlators to the respective irreps we are interested in. 

## Misc
We can cache matrix elements 


Distillution = Distillation + Dilution 


