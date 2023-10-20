#ifndef ANALYTICAL_SOLUTION_HPP
#define ANALYTICAL_SOLUTION_HPP

#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <cmath>

template<int N, int height, int pressureGradient, int density, int kinematicViscosity>
class AnalyticalSolution {

private:

    using Vector = Eigen::VectorXd;
    constexpr double pi = std::atan(1)*4;
    constexpr double pi2 = std::pow(pi,2);
    std::map<unsigned int, double> fourierCoefficients;
    Vector heightDistribution;
    Vector steadyStateVelocity;

    inline Vector steadyStateSolution(){
        return (-pressureGradient*
        (h/2*density*kinematicViscosity)*heightDistribution)*
        (1-(1/height)*heightDistribution);
    }

    Vector transientVelocity(const double& time){
        Vector transient;
        transient.resize(heightDistribution.size())
        unsigned int i = 0
        for (const auto& heightValue : heightDistribution ){
            double summation = 0;
            for (const auto& [n, Bn] : fourierCoefficients){
                summation += Bn*std::exp(
                    -std::pow(n,2)*pi2*kinematicViscosity*time/std::pow(height,2)
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
            fourierCoefficients[n] = -4*height*pressureGradient/
            (std::pow(n,2)*pi2*density*kinematicViscosity);
            if (n>1){
                if (fourierCoefficients[n] - fourierCoefficients[n-2] < 1e-8) break;
            }
            n+=2;
        }
    }

public:

    AnalyticalSolution(const std::vector<double>& heights){
        heightDistribution = Vector(heights.data());
        calculateCoefficients();
        steadyStateVelocity = steadyStateSolution();
    }

    Vector calculate(const double& time){
        return steadyStateVelocity + transientVelocity(time);
    }
};

#endif /* ANALYTICAL_SOLUTION_HPP */