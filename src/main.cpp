#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <NumericalSolution.hpp>
#include <AnalyticalSolution.hpp>

#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;

auto eigenToStandard(auto vec){
    std::vector<double> uvec(vec.data(), vec.data() + vec.size());
    return uvec;
}


int main(int, char**){

    constexpr int pressureGradient = -8;
    constexpr int kinematicViscosity = 1;
    constexpr int density = 1;
    constexpr int height = 1;
    constexpr int N = 101;
    constexpr double deltaY = height/(N-1);

    auto heights = std::vector<double>(N);
    std::generate(heights.begin(), heights.end(), [n = 0, &deltaY]() mutable { return n++ * deltaY; });

    AnalyticalSolution<N,height,pressureGradient,density,kinematicViscosity> analyticalPoiseuille(heights);
    NumericalSolution<N> numericalPoiseuille(pressureGradient, kinematicViscosity, density, height);

    double time1 = 0.025, time2 = 0.05, time3=0.1, time4=0.4;
    
    auto analyticaVelocitiesTime1 = analyticalPoiseuille.calculate(time1);
    numericalPoiseuille.simulate(time1);
    auto numericalVelocitiesTime1 = numericalPoiseuille.velocities();

    auto analyticaVelocitiesTime2 = analyticalPoiseuille.calculate(time2);
    numericalPoiseuille.simulate(time2);
    auto numericalVelocitiesTime2 = numericalPoiseuille.velocities();

    auto analyticaVelocitiesTime3 = analyticalPoiseuille.calculate(time3);
    numericalPoiseuille.simulate(time3);
    auto numericalVelocitiesTime3 = numericalPoiseuille.velocities();

    auto analyticaVelocitiesTime4 = analyticalPoiseuille.calculate(time4);
    numericalPoiseuille.simulate(time4);
    auto numericalVelocitiesTime4 = numericalPoiseuille.velocities();

    auto differenceTime1 = analyticaVelocitiesTime1 - numericalVelocitiesTime1; 
    auto differenceTime2 = analyticaVelocitiesTime2 - numericalVelocitiesTime2;
    auto differenceTime3 = analyticaVelocitiesTime3 - numericalVelocitiesTime3;
    auto differenceTime4 = analyticaVelocitiesTime4 - numericalVelocitiesTime4;

    plt::named_plot("time = 0.025 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime1), heights);
    plt::named_plot("time = 0.05 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime2), heights);
    plt::named_plot("time = 0.1 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime3), heights);
    plt::named_plot("time = 0.4 (Numerical Solution)", eigenToStandard(numericalVelocitiesTime4), heights);
    
    plt::title("Poiseuille flow Analytical vs Numerical Solution");
    plt::legend();
    plt::show();

    return 0;
}
