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
/*
    auto numericalTime1 = eigenToStandard(numericalVelocitiesTime1);
    auto numericalTime2 = eigenToStandard(numericalVelocitiesTime2);
    auto numericalTime3 = eigenToStandard(numericalVelocitiesTime3);
    auto numericalTime4 = eigenToStandard(numericalVelocitiesTime4);

    auto analyticalTime1 = eigenToStandard(analyticalVelocitiesTime1);
    auto analyticalTime2 = eigenToStandard(analyticalVelocitiesTime2);
    auto analyticalTime3 = eigenToStandard(analyticalVelocitiesTime3);
    auto analyticalTime4 = eigenToStandard(analyticalVelocitiesTime4);
*/
    plt::named_plot("time = 0.025 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime1), heights, "bo");
    plt::named_plot("time = 0.025 (Analytical Solution)", eigenToStandard(analyticalVelocitiesTime1), heights,"b-");
    plt::named_plot("time = 0.05 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime2), heights, "ro");
    plt::named_plot("time = 0.05 (Analytical Solution)", eigenToStandard(analyticalVelocitiesTime2), heights, "r-");
    plt::named_plot("time = 0.1 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime3), heights, "go");
    plt::named_plot("time = 0.1 (Analytical Solution)", eigenToStandard(analyticalVelocitiesTime3), heights, "g-");
    plt::named_plot("time = 0.4 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime4), heights, "yo");
    plt::named_plot("time = 0.4 (Analytical Solution)", eigenToStandard(analyticalVelocitiesTime4), heights, "y-");
    
    plt::title("Poiseuille flow Analytical vs Numerical Solution");
    plt::legend();
    plt::save("./Comparison.pdf");
    plt::show();

    plt::named_plot("time = 0.025 (Difference)", eigenToStandard(differenceTime1), heights);
    plt::named_plot("time = 0.05 (Difference)", eigenToStandard(differenceTime2), heights);
    plt::named_plot("time = 0.1 (Difference)", eigenToStandard(differenceTime3), heights);
    plt::named_plot("time = 0.4 (Difference)", eigenToStandard(differenceTime4), heights);
    plt::title("Poiseuille flow Difference");
    plt::legend();
    plt::save("./Difference.pdf");
    plt::show();

    return 0;
}
