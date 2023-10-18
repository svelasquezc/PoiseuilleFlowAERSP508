#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <NumericalSolution.hpp>

#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;



int main(int, char**){

    constexpr double pressureGradient = -8;
    constexpr double kinematicViscosity = 1;
    constexpr double density = 1;
    constexpr double height = 1;
    constexpr int N = 101;
    constexpr double deltaY = height/(N-1);
    double time = 0;

    auto heights = std::vector<double>(N);
    std::generate(heights.begin(), heights.end(), [n = 0, &deltaY]() mutable { return n++ * deltaY; });

    NumericalSolution<N> numericalPoiseuille(pressureGradient, kinematicViscosity, density, height);

    double time1 = 0.025, time2 = 0.05, time3=0.1, time4=0.4;
    
    numericalPoiseuille.simulate(time1);
    auto numericalVelocitiesTime1 = numericalPoiseuille.velocities();
    numericalPoiseuille.simulate(time2);
    auto numericalVelocitiesTime2 = numericalPoiseuille.velocities();
    numericalPoiseuille.simulate(time3);
    auto numericalVelocitiesTime3 = numericalPoiseuille.velocities();
    numericalPoiseuille.simulate(time4);
    auto numericalVelocitiesTime4 = numericalPoiseuille.velocities();

    plt::named_plot("time = 0.025 (Numerical Solution)", numericalVelocitiesTime1, heights);
    plt::named_plot("time = 0.05 (Numerical Solution)", numericalVelocitiesTime2, heights);
    plt::named_plot("time = 0.1 (Numerical Solution)", numericalVelocitiesTime3, heights);
    plt::named_plot("time = 0.4 (Numerical Solution)", numericalVelocitiesTime4, heights);
    
    plt::title("Poiseuille flow Analytical vs Numerical Solution");
    plt::legend();
    plt::show();

    return 0;
}
