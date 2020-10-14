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
            ``rcr.priorsTypes`` *object*. The type of priors that you're applying to your model (see :class:`trk.priorsTypes` and :ref:`priorstypes`).
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
            *float*. Best fit extrinsic scatter/sample variance :math:`\sigma_x`, or *slop*, of the model along the :math:`x` direction. (Uncertainty in the dataset that cannot soleley be accounted for by the :math:`x` error bars).
        )mydelimiter")

        .def_readwrite("slop_y", &Results::slop_y, R"mydelimiter(
            *float*. Best fit extrinsic scatter/sample variance :math:`\sigma_y`, or *slop*, of the model along the :math:`y` direction. (Uncertainty in the dataset that cannot soleley be accounted for by the :math:`y` error bars).
        )mydelimiter")

        .def_readwrite("optimum_scale", &Results::optimumScale, R"mydelimiter(
            *float*. The determined global optimal fitting scale :math:`s_0` for the model/dataset, where :math:`s_0` is a factor to multiply the `y`-axis of the data by.
        )mydelimiter")

        .def_readwrite("minimum_scale", &Results::minimumScale, R"mydelimiter(
            *float*. The determined minimum fitting scale :math:`a` for the model/dataset, where :math:`a` is the fitting scale where the :math:`x`-slop :math:`\sigma_x` goes to zero.
        )mydelimiter")

        .def_readwrite("maximum_scale", &Results::maximumScale, R"mydelimiter(
            *float*. The determined maximum fitting scale :math:`b` for the model/dataset, where :math:`b` is the fitting scale where the :math:`y`-slop :math:`\sigma_y` goes to zero.
        )mydelimiter")

        .def_readwrite("fitness", &Results::fitness, R"mydelimiter(
            *float*. The relative fitness :math:`-2\ln\mathcal{L}^\text{TRK}` of the best fit, analogous to a :math:`\chi^2` statistic, where :math:`\mathcal{L}^\text{TRK}` is the TRK likelihood function. Useful for comparing different fits obtained for the same model, or possibly model comparison.
        )mydelimiter")

        .def_readwrite("pivots", &Results::pivots, R"mydelimiter(
            *1D list/array_like of floats*. The pivot points :math:`x_p^i` of the linearized segments of the model, that are optimized to minimize the correlations between the slope and intercept of each segment.
        )mydelimiter")

        .def_readwrite("model_parameters", &Results::bestFitParams, R"mydelimiter(
            *1D list/array_like of floats*. The best fit parameters :math:`\vartheta_m` for the model (that minimize :math:`-2\ln\mathcal{L}^\text{TRK}`, where :math:`\mathcal{L}^\text{TRK}` is the TRK likelihood function).
        )mydelimiter")

        .def_readwrite("model_parameter_uncertainties", &Results::bestFitParamUncertainties, R"mydelimiter(
            *2D list/array_like of floats*. The (possibly asymmetric) :math:`(+, -)1\sigma` widths/uncertainties of the best fit model parameters :math:`\vartheta_m`. Array organized as :math:`\left[(\sigma_{\vartheta_1,-}, \sigma_{\vartheta_1,+}), \ldots\right]`. 
        )mydelimiter");

        // main class
        py::class_<TRK>(m, "TRK", R"mydelimiter(
                *class*. Master class used to initialize and run TRK fitting procedures.

                Constructor arguments:

                Parameters
                ----------
                f : function
                    Model function :math:`y(x|\vec{\theta})` to fit to data, where :math:`x` is the independent variable (float) and :math:`\vec{\theta}` is an
                    :math:`M`-dimensional list/array_like of model parameters. Arguments for ``f`` must follow this prototype:

                    Parameters
                    ----------
                    x : float
                        Independent variable of model
                    params : list/array_like, 1D
                        Parameters of model

                    Returns
                    -------
                    y : float
                        Model evaluated at the corresponding values of ``x`` and ``params``.

                df : function
                    First derivative of the model function, :math:`\frac{dy(x|\vec{\theta})}{dx}`, with respect to the independent variable :math:`x`,  where :math:`x` is the independent variable (float) and :math:`\vec{\theta}` is an
                    :math:`M`-dimensional list/array_like of model parameters. Arguments for ``df`` must follow this prototype:

                    Parameters
                    ----------
                    x : float
                        Independent variable of model
                    params : list/array_like, 1D
                        Parameters of model

                    Returns
                    -------
                    dy : float
                        First :math:`x`-derivative of the model evaluated at the corresponding values of ``x`` and ``params``.

                ddf : function
                    Second derivative of the model function, :math:`\frac{d^2y(x|\vec{\theta})}{dx^2}`, with respect to the independent variable :math:`x`,  where :math:`x` is the independent variable (float) and :math:`\vec{\theta}` is an
                    :math:`M`-dimensional list/array_like of model parameters. Arguments for ``ddf`` must follow this prototype:

                    Parameters
                    ----------
                    x : float
                        Independent variable of model
                    params : list/array_like, 1D
                        Parameters of model

                    Returns
                    -------
                    ddy : float
                        Second :math:`x`-derivative of the model evaluated at the corresponding values of ``x`` and ``params``.

                xdata : list/array_like, 1D
                    Independent variable data to fit model to.

                ydata : list/array_like, 1D
                    Dependent variable (model function evaluation) data to fit model to.

                error_x : list/array_like, 1D
                    Error bars/:math:`x`-uncertainties to be applied to dataset (see :ref:`errorbars`).

                error_y : list/array_like, 1D
                    Error bars/:math:`y`-uncertainties to be applied to dataset (see :ref:`errorbars`).

                guess : list/array_like, 1D
                    Guess for best fit values of model parameters :math:`\vec{\theta}` (for the fitting algorithm).

                weights : list/array_like, optional, 1D
                    Optional weights to be applied to dataset (see :ref:`weighting`).


                tol : float, optional
                    Default: ``1e-6``. Convergence tolerance of modified Gauss-Newton fitting algorithm.

                has_priors : bool, optional
                    Default: ``False``. Set to ``True`` if you're going to apply statistical priors to your model parameters 
                    (see :ref:`priors`; you'll also need to create an instance of :class:`rcr.Priors` and set the ``priors`` attribute of this instance of ``FunctionalForm`` equal to it).

                pivot_function : function, optional
                    Default: ``None``. Function that returns the pivot point of some linearized model (see :ref:`pivots`). Must be of the form/prototype of:

                    Parameters
                    ----------
                    xdata : list/array_like, 1D or 2D
                        :math:`n`-dimensional independent variable data to fit model to; same as above``xdata``.
                    weights : list/array_like, optional, 1D
                        Optional weights to be applied to dataset (see :ref:`weighting`).
                    f : function
                        Model function; same as above ``f``.
                    params : list/array_like, 1D
                        Parameters of model

                    Returns
                    -------
                    pivot : float or 1D list/array_like
                        Pivot point(s) of the model; (``float`` if you're using a one-dimensional model/independent variable, ``list/array_like`` if :math:`n`-dimensional.)

                    However, note that all arguments need to be actually used for the pivot point computation. For example,
                    a simple linear model :math:`y(x|b,m) = b + m(x-x_p)` has a pivot point found by :math:`x_p=\sum_iw_ix_i/\sum_iw_i`, where
                    :math:`w_i` are the weights of the datapoints.
                
                pivot_guess : float or 1D list/array_like, optional
                    Initial guess for the pivot point(s) of the model (``float`` if you're using a one-dimensional model/independent variable, ``list/array_like`` if :math:`n`-dimensional; see :ref:`pivots`).
            )mydelimiter")
            // constructors
            .def(py::init(&getTRKObject))
            .def(py::init<>())

}
