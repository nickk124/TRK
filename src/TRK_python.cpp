/*
 Trotter Reichart Konz (TRK) Regression Official Codebase
 Creator/Author: Nick C. Konz
 See license at https://github.com/nickk124/TRK

This file houses all of the TRK functionality that needs to be exposed to Python.

build cmd: (only tested on mac atm):
    c++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` TRK_python.cpp TRK.cpp -o trk`python3-config --extension-suffix`
*/
#include "TRK.h"
#include <pybind11/pybind11.h> // pybind header files are within ./pybind11/include/pybind11/
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace TRKLib;

// python binding functions
PYBIND11_MODULE(trk, m) { // trk is module name, m is docstring instance
    m.doc() = "TRK (Trotter Reichart Konz Worst-Case Regression) Package API Details.";

    // note: all nested classes within TRK class are public.

    // ENUMS
    py::enum_<priorType>(m, "priorType", py::arithmetic(), "Types of prior probability density functions that can be applied to model parameters.")
        .value("CUSTOM", CUSTOM, "Custom, function-defined prior probability density functions(s).")
        .value("CUSTOM_JOINT", CUSTOM_JOINT, "Custom, function-defined joint (non-independent) prior probability density function.")
        .value("GAUSSIAN", GAUSSIAN, "Gaussian (normal) prior probability density function(s).")
        .value("CONSTRAINED", CONSTRAINED, "Bounded/hard-constrained prior probability density function(s).")
        .value("MIXED", MIXED, "A mixture of gaussian (normal), hard-constrained, and uninformative (uniform/flat) prior probability density functions.")
        .export_values();

    // CLASSES
    // parameter prior probability distributions
    py::class_<Priors>(m, "Priors", R"mydelimiter(
            *class*. Class that encapsulates probabalistic priors to beß applied to model parameters when using model-fitting/functional form RCR (see :ref:`priors` for an example).

            Constructor arguments:

            Parameters
            ----------
            priorType : :class:`rcr.priorsTypes`
                The type of priors that you're applying to your model (see :class:`rcr.priorsTypes` and :ref:`priorstypes`).

            priorsPDFs : list of functions, optional 2nd argument
                An ordered list of custom model parameter priors functions; each function takes in a model parameter as argument and returns the prior probability density function for that parameter (see :ref:`priors` for an example).
            
            jointPriorsPDF: function, optional 2nd argument
                A custom joint prior probability function, i.e. it takes an argument of a vector of the model params (including slop, last), and returns the joint prior probability density (float).

            gaussianParams : 2D list/array_like, optional 2nd argument
                A list that contains lists of mu and sigma for the Gaussian prior of each param. If no prior, then just use NaNs (see :ref:`priors` for an example).
            
            paramBounds : 2D list/array_like, optional 2nd argument (or 3rd, for the case of ``rcr.MIXED_PRIORS``)
                A list that contains lists of the lower and upper hard bounds of each param. If not bounded, use NaNs, and if there's only one bound, use NaN for the other bound (see :ref:`priors` for an example).
        )mydelimiter")
        
        // constructors
        .def(py::init< priorType, std::vector < std::vector <double> > >()) // only Gaussian or only bounded/hard constraints
        .def(py::init< priorType, std::vector < std::vector <double> >, std::vector < std::vector <double> > >()) // mixed priors
        .def(py::init< priorType, std::vector < std::function <double(double)> > > ()) // custom priors
        .def(py::init< priorType, std::function <double(std::vector <double>)> >()) // custom joint priors
        .def(py::init<>())

        // members
        .def_readwrite("priorType", &Priors::priorType, R"mydelimiter(
            ``rcr.priorsTypes`` *object*. The type of priors that you're applying to your model (see :class:`rcr.priorsTypes` and :ref:`priorstypes`).
        )mydelimiter")
        .def_readwrite("gaussianParams", &Priors::gaussianParams, R"mydelimiter(
            *2D list/array_like of floats*. A list that contains lists of mu and sigma for the Gaussian prior of each param. If no prior, then just use NaNs (see :ref:`priors` for an example).
        )mydelimiter")
        .def_readwrite("paramBounds", &Priors::paramBounds, R"mydelimiter(
            *2D list/array_like of floats*. A list that contains lists of the lower and upper hard bounds of each param. If not bounded, use NaNs, and if there's only one bound, use NaN for the other bound (see :ref:`priors` for an example).
        )mydelimiter")
        .def_readwrite("priorsPDFs", &Priors::priorsPDFs, R"mydelimiter(
            *list/array_like of functions*. An ordered list of custom model parameter priors functions; each function takes in a model parameter as argument and returns the prior probability density function for that parameter (see :ref:`priors` for an example).
        )mydelimiter")
        .def_readwrite("joint_priors_PDF", &Priors::jointPriorsPDF, R"mydelimiter(
            *function*. A custom joint prior probability function, i.e. it takes an argument of a vector of the model params (including slop, last), and returns the joint prior probability density (float).
        )mydelimiter");

        // fit results
        py::class_<Results>(m, "Results", "Results from running TRK fit algorithms.")
            .def_readwrite("slop_x", &Results::slop_x, R"mydelimiter(

            )mydelimiter")

            .def_readwrite("slop_y", &Results::slop_y, R"mydelimiter(

            )mydelimiter")

            .def_readwrite("optimum_scale", &Results::optimumScale, R"mydelimiter(

            )mydelimiter")

            .def_readwrite("minimum_scale", &Results::minimumScale, R"mydelimiter(

            )mydelimiter");
}