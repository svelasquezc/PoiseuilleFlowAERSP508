#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <Eigen/Sparse>
#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;



int main(int, char**){

    using SparseMatrix = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    using Vector = Eigen::VectorXd;
    using Solver =  Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double, int> >;

    constexpr int pressureGradient = -8;
    constexpr unsigned int kinematicViscosity = 1;
    constexpr unsigned int density = 1;

    constexpr double timeDelta = 0.001;

    constexpr double deltaY = 0.01;
    constexpr auto N = 101;

    constexpr auto lambda = kinematicViscosity*timeDelta/std::pow(deltaY, 2);


    /* Initializing the whole vector as the first rhs term */
    Vector gradientVector(N);
    gradientVector.array() = -timeDelta*pressureGradient/density;
    gradientVector[0] = 0;
    gradientVector[N-1] = 0;

    Vector velocityVector(N);
    velocityVector.setZero();

    SparseMatrix lhsMatrix;
    {
        constexpr auto maxNonZeros = 3*N-4;
        auto nonZeros = std::vector<Triplet>();
        nonZeros.reserve(maxNonZeros);

        nonZeros.push_back({0,0,1});
        

        for (int i=1; i<N-1; ++i){
            nonZeros.push_back({i,i-1,-lambda/2});
            nonZeros.push_back({i,i,lambda+1});
            nonZeros.push_back({i,i+1,-lambda/2});
        }
        nonZeros.push_back({N-1,N-1,1});
        
        lhsMatrix.resize(N,N);
        lhsMatrix.reserve(maxNonZeros);
        lhsMatrix.setFromTriplets(nonZeros.begin(), nonZeros.end());
    }

    auto rhsExplicit = [N, lambda](const Vector& velocityVector){
        Vector rhs(N);
        rhs.setZero();

        for(int i=1; i<N-1;++i){
            rhs[i]= lambda*velocityVector[i-1]/2 + (1-lambda)*velocityVector[i] + lambda*velocityVector[i+1]/2;
        }

        return rhs;
    };

    Solver solver;
    solver.preconditioner().setDroptol(0.00000001);
    solver.compute(lhsMatrix);

    double time = 0;
    double givenTime = 0.4;

    while (time <= givenTime){

        velocityVector = solver.solve(gradientVector + rhsExplicit(velocityVector));

        time += timeDelta;
    }

    auto heights = std::vector<double>(N);
    std::vector<double> uvec(velocityVector.data(), velocityVector.data() + velocityVector.size());

    std::generate(heights.begin(), heights.end(), [n = 0, &deltaY]() mutable { return n++ * deltaY; });

    plt::plot(heights, uvec);
    plt::show();

    return 0;
}
