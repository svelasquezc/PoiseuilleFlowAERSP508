#ifndef ANALYTICAL_SOLUTION_HPP
#define ANALYTICAL_SOLUTION_HPP

#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <cmath>


constexpr double pi = 3.14159265358979323846;


template<int N>
class AnalyticalSolution {

private:

    using Vector = Eigen::VectorXd;
    std::map<unsigned int, double> fourierCoefficients;
    Vector heightDistribution;
    Vector steadyStateVelocity;
    double height;
    double pressureGradient;
    double density;
    double kinematicViscosity;


    Vector steadyStateSolution(){
        Vector ones;
        ones.resize(N);
        ones.array() = 1.0;
        return (-pressureGradient*
        (height/(2*density*kinematicViscosity)
        )*heightDistribution).cwiseProduct(ones - (1/height)*heightDistribution);
    }

    Vector transientVelocity(const double& time){
        Vector transient;
        transient.resize(heightDistribution.size());
        unsigned int i = 0;
        for (const auto& heightValue : heightDistribution ){
            double summation = 0;
            for (const auto& [n, Bn] : fourierCoefficients){
                summation += Bn*std::exp(
                    -std::pow(n,2)*std::pow(pi,2)*kinematicViscosity*time/std::pow(height,2)
                    )*std::sin(n*pi*heightValue/height);
            }
            transient[i] = summation;
            ++i;
        }
        return transient;
    }

    void calculateCoefficients(){
        unsigned int n = 1;
        while(true){
            fourierCoefficients[n] = 4*std::pow(height,2)*pressureGradient/
            (std::pow(n,3)*std::pow(pi,3)*density*kinematicViscosity);
            if (n>1) if (std::abs(fourierCoefficients[n] - fourierCoefficients[n-2]) < 1e-8) break;
            n+=2;
        }
    }

public:

    AnalyticalSolution(std::vector<double>& heights, double _height, double _pressureGradient, double _density, double _kinematicViscosity){

        height = _height;
        pressureGradient = _pressureGradient;
        density = _density;
        kinematicViscosity = _kinematicViscosity;

        heightDistribution = Eigen::Map<Vector>(heights.data(), heights.size());

        steadyStateVelocity = steadyStateSolution();
        calculateCoefficients();
    }

    Vector calculate(const double& time){
        return steadyStateVelocity + transientVelocity(time);
    }
};

#endif /* ANALYTICAL_SOLUTION_HPP */