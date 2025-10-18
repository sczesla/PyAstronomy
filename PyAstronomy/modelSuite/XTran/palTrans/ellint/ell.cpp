#include "boost/math/special_functions/ellint_1.hpp"
#include "boost/math/special_functions/ellint_2.hpp"
#include "boost/math/special_functions/ellint_3.hpp"
#include "boost/python.hpp"

using namespace std;
using namespace boost::math;
using namespace boost::python;

  struct Ell
  {
    double ell1(double k) {return ellint_1(k);};
    double ell2(double k) {return ellint_2(k);};
    double ell3(double n, double k) {return ellint_3(k,n);};
  };

  struct ell_pickle_suite : pickle_suite
  {
    static tuple getinitargs(const Ell& e)
    {
        return boost::python::make_tuple();
    }
  };


  BOOST_PYTHON_MODULE(ell)
  {
      class_<Ell>("ell")
        .def_pickle(ell_pickle_suite())
        .def("ell1", &Ell::ell1)
        .def("ell2", &Ell::ell2)
        .def("ell3", &Ell::ell3)
      ;
  }
