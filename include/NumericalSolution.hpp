#ifndef NUMERICAL_SOLUTION_HPP
#define NUMERICAL_SOLUTION_HPP

#include <cmath>
#include <Eigen/Sparse>

template<int N>
class NumericalSolution
{
private:

    /*
    Aliases from Eigen types
    */
    using SparseMatrix = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    using Vector = Eigen::VectorXd;
    using Solver =  Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double, int> >; // Preconditioned solver


    // Parameters 
    const double timeDelta = 0.001;
    mutable double deltaY = 0.01; // height/(N-1)

    Vector gradientVector, velocityVector;
    SparseMatrix lhsMatrix;
    const int maxNonZeros = 3*N-4;

    mutable double lambda = 0;
    double time;

    Solver solver;

    Vector rhsExplicit (const Vector& velocityVector){
        Vector rhs(N);
        rhs.setZero();

        for(int i=1; i<N-1;++i){
            rhs[i]= lambda*velocityVector[i-1]/2 + (1-lambda)*velocityVector[i] + lambda*velocityVector[i+1]/2;
        }

        return rhs;
    };

public:
    NumericalSolution(
        double pressureGradient, 
        double kinematicViscosity, 
        double density, 
        double height
    ){

        gradientVector.resize(N);
        velocityVector.resize(N);
        lhsMatrix.resize(N,N);
        deltaY = height/(N-1);

        lambda = kinematicViscosity*timeDelta/std::pow(deltaY, 2);

        gradientVector.array() = -timeDelta*pressureGradient/density;
        gradientVector[0] = 0; //Non-slip condition
        gradientVector[N-1] = 0;//Non-slip condition

        time = 0;
        velocityVector.setZero(); //Initial Condition

        /*
        Matrix initialization: As the values of pressure gradient, viscosity, density, and the
                               delta of time are constant over the whole simulation, the lambda
                               parameter is constant, therefore the matrix remains unchanged, only 
        */
        {
            auto nonZeros = std::vector<Triplet>();
            nonZeros.reserve(maxNonZeros);

            nonZeros.push_back({0,0,1});//Non-slip condition
            

            for (int i=1; i<N-1; ++i){
                nonZeros.push_back({i,i-1,-lambda/2});
                nonZeros.push_back({i,i,lambda+1});
                nonZeros.push_back({i,i+1,-lambda/2});
            }
            nonZeros.push_back({N-1,N-1,1});//Non-slip condition
            
            lhsMatrix.reserve(maxNonZeros);
            lhsMatrix.setFromTriplets(nonZeros.begin(), nonZeros.end());
        }

        //Initialize solver with matrix already defined
        solver.preconditioner().setDroptol(0.00000001);
        solver.compute(lhsMatrix);
    }

    void startSimulation(){
        time = 0;
        velocityVector.setZero();
    }

    void simulate(const double& goalTime){

        while (time <= goalTime){
            velocityVector = solver.solve(gradientVector + rhsExplicit(velocityVector));
            time += timeDelta;
        }
    }

    Vector velocities() const{
        return velocityVector;
    }

    ~NumericalSolution() = default;
};

#endif /* NUMERICAL_SOLUTION_HPP */