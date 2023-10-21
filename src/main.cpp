#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <AnalyticalSolution.hpp>
#include <NumericalSolution.hpp>

#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;

auto eigenToStandard(Eigen::VectorXd& vec){
    std::vector<double> uvec(vec.data(), vec.data() + vec.size());
    return uvec;
}


int main(int, char**){

    double pressureGradient = -8;
    double kinematicViscosity = 1;
    double density = 1;
    double height = 1;
    constexpr int N = 101;
    double deltaY = height/(N-1);

    

    NumericalSolution<N> numericalPoiseuille(
        pressureGradient,
        kinematicViscosity,
        density,
        height
        );

    double time1 = 0.025, time2 = 0.05, time3=0.1, time4=0.4;
    
    numericalPoiseuille.simulate(time1);
    auto numericalVelocitiesTime1 = numericalPoiseuille.velocities();
    numericalPoiseuille.simulate(time2);
    auto numericalVelocitiesTime2 = numericalPoiseuille.velocities();
    numericalPoiseuille.simulate(time3);
    auto numericalVelocitiesTime3 = numericalPoiseuille.velocities();
    numericalPoiseuille.simulate(time4);
    auto numericalVelocitiesTime4 = numericalPoiseuille.velocities();

    auto heights = std::vector<double>(N);
    std::generate(heights.begin(), heights.end(), [n = 0, &deltaY]() mutable { return n++ * deltaY; });

    AnalyticalSolution<N> analyticalPoiseuille(heights, height, pressureGradient, density, kinematicViscosity);

    auto analyticalVelocitiesTime1 = analyticalPoiseuille.calculate(time1);
    auto analyticalVelocitiesTime2 = analyticalPoiseuille.calculate(time2);
    auto analyticalVelocitiesTime3 = analyticalPoiseuille.calculate(time3);
    auto analyticalVelocitiesTime4 = analyticalPoiseuille.calculate(time4);

    Eigen::VectorXd differenceTime1 = analyticalVelocitiesTime1 - numericalVelocitiesTime1; 
    Eigen::VectorXd differenceTime2 = analyticalVelocitiesTime2 - numericalVelocitiesTime2;
    Eigen::VectorXd differenceTime3 = analyticalVelocitiesTime3 - numericalVelocitiesTime3;
    Eigen::VectorXd differenceTime4 = analyticalVelocitiesTime4 - numericalVelocitiesTime4;

    plt::named_plot("Analytical, time = 0.025", eigenToStandard(analyticalVelocitiesTime1), heights,"b-");
    plt::named_plot("Numerical, time = 0.025", eigenToStandard(numericalVelocitiesTime1), heights, "bx");
    plt::named_plot("Analytical, time = 0.05", eigenToStandard(analyticalVelocitiesTime2), heights, "r-");
    plt::named_plot("Numerical, time = 0.05", eigenToStandard(numericalVelocitiesTime2), heights, "rx");
    plt::named_plot("Analytical, time = 0.1", eigenToStandard(analyticalVelocitiesTime3), heights, "g-");
    plt::named_plot("Numerical, time = 0.1", eigenToStandard(numericalVelocitiesTime3), heights, "gx");
    plt::named_plot("Analytical, time = 0.4", eigenToStandard(analyticalVelocitiesTime4), heights, "y-");
    plt::named_plot("Numerical, time = 0.4", eigenToStandard(numericalVelocitiesTime4), heights, "yx");
    
    plt::title("Poiseuille flow Analytical vs Numerical Solution");
    plt::xlabel("Velocity");
    plt::ylabel("Height");
    plt::legend();
    plt::save("./Comparison.pdf");
    plt::show();

    plt::named_plot("time = 0.025", eigenToStandard(differenceTime1), heights);
    plt::named_plot("time = 0.05", eigenToStandard(differenceTime2), heights);
    plt::named_plot("time = 0.1", eigenToStandard(differenceTime3), heights);
    plt::named_plot("time = 0.4", eigenToStandard(differenceTime4), heights);
    plt::xlabel("Error (Diference) in Velocity");
    plt::ylabel("Height");
    plt::title("Poiseuille flow Difference Error (Analytical - Numerical)");
    plt::legend();
    plt::save("./Difference.pdf");
    plt::show();

    return 0;
}
