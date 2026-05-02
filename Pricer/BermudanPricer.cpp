#include "pch.h"
#include "BermudanPricer.h"
#include <cmath>
#include <Eigen/Dense>
#include <functional>

BermudanPricer::BermudanPricer(RandomProcess* _process, PayOff* _payoff, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps, const std::vector<double>& _exerciseDates) :
    MCPricer(_process, _payoff, _rate, _maturity, _nbSim, _nbSteps),
    exerciseDates(_exerciseDates)
{
}

BermudanPricer::~BermudanPricer()
{
}

std::vector<std::vector<int>> BermudanPricer::generateBasisDegrees(int D, int p) {


    std::vector<std::vector<int>> basisDegrees;
    std::vector<int> current_degrees(D, 0);
    
    std::function<void(int, int)> generate = [&](int dim, int remainingDegree) {
        if (dim == D - 1) {
            for (int i = 0; i <= remainingDegree; ++i) {
                current_degrees[dim] = i;
                basisDegrees.push_back(current_degrees);
            }
        } else {
            for (int i = 0; i <= remainingDegree; ++i) {
                current_degrees[dim] = i;
                generate(dim + 1, remainingDegree - i);
            }
        }
    };
    
    generate(0, p);
    return basisDegrees;
}

Eigen::VectorXd BermudanPricer::regressionFit(const std::vector<std::vector<double>>& X, const std::vector<double>& Y, const std::vector<std::vector<int>>& basisDegrees) {
    int nPaths = X.size();
    int nBasis = basisDegrees.size();
    int D = X[0].size();
    
    if (nPaths < nBasis) {
        return Eigen::VectorXd::Zero(nBasis);
    }
    
    Eigen::MatrixXd Phi(nPaths, nBasis);
    Eigen::VectorXd Yvec(nPaths);
    
    for (int row = 0; row < nPaths; ++row) {
        Yvec(row) = Y[row];
        for (int col = 0; col < nBasis; ++col) {
            double term = 1.0;
            for (int d = 0; d < D; ++d) {
                if (basisDegrees[col][d] > 0) {
                    term *= std::pow(X[row][d], basisDegrees[col][d]);
                }
            }
            Phi(row, col) = term;
        }
    }
    
    return Phi.colPivHouseholderQr().solve(Yvec);
}

double BermudanPricer::Price() {
    int nDates = exerciseDates.size();
    if (nDates == 0) return 0.0;
    
    int D = Process->GetDimension();

    std::vector<int> exerciseIndices;
    for (double d : exerciseDates) exerciseIndices.push_back(static_cast<int>((d / maturity) * nbSteps));

    std::vector<std::vector<double>> PayoffsAtDates(nbSim, std::vector<double>(nDates));
    std::vector<std::vector<std::vector<double>>> StatesAtDates(nbSim, std::vector<std::vector<double>>(nDates, std::vector<double>(D)));

    for (size_t sim = 0; sim < nbSim; ++sim) {
        Process->Simulate(0.0, maturity, nbSteps);
        
        for (int k = 0; k < nDates; ++k) {
            int idx = exerciseIndices[k];
            
            std::vector<SinglePath*> truncPaths(D);
            for (int d = 0; d < D; ++d) {
                SinglePath* original = Process->GetPath(d);
                truncPaths[d] = new SinglePath(0.0, exerciseDates[k], idx);
                const auto& vals = original->GetAllValues();
                
                int limit = std::min(idx, (int)vals.size() - 1);
                for(int step = 0; step <= limit; ++step) {
                    truncPaths[d]->InsertValue(vals[step]);
                }
                StatesAtDates[sim][k][d] = vals[limit];
            }
            
            PayoffsAtDates[sim][k] = (*Payoff)(truncPaths);
            
            for (int d = 0; d < D; ++d) {
                delete truncPaths[d];
            }
        }
    }

    std::vector<double> V(nbSim);
    for (size_t sim = 0; sim < nbSim; ++sim) {
        V[sim] = PayoffsAtDates[sim][nDates - 1];
    }

    auto basisDegrees = generateBasisDegrees(D, PolynomialDegree);
    int nBasis = basisDegrees.size();

    for (int k = nDates - 2; k >= 0; --k) {
        double dt = exerciseDates[k+1] - exerciseDates[k];
        double discount = std::exp(-rate * dt);
        
        std::vector<std::vector<double>> X;
        std::vector<double> Y;
        std::vector<int> itmPaths;
        
        for (size_t sim = 0; sim < nbSim; ++sim) {
            double payoff = PayoffsAtDates[sim][k];
            if (payoff > 0.0) {
                X.push_back(StatesAtDates[sim][k]);
                Y.push_back(V[sim] * discount);
                itmPaths.push_back(sim);
            } else {
                V[sim] = V[sim] * discount;
            }
        }
        
        if (!X.empty()) {
            Eigen::VectorXd model = regressionFit(X, Y, basisDegrees);
            
            for (size_t i = 0; i < itmPaths.size(); ++i) {
                int sim = itmPaths[i];
                double continuation = 0.0;
                for (int col = 0; col < nBasis; ++col) {
                    double term = 1.0;
                    for (int d = 0; d < D; ++d) {
                        if (basisDegrees[col][d] > 0) {
                            term *= std::pow(X[i][d], basisDegrees[col][d]);
                        }
                    }
                    continuation += model(col) * term;
                }
                
                double exercise = PayoffsAtDates[sim][k];
                if (exercise >= continuation) {
                    V[sim] = exercise;
                } else {
                    V[sim] = Y[i];
                }
            }
        }
    }

    double discount0 = std::exp(-rate * exerciseDates[0]);
    double sum = 0.0;
	double sumSq = 0.0;
    for (size_t sim = 0; sim < nbSim; ++sim) {
        double discountedValue = V[sim] * discount0;
        sum += discountedValue;
        sumSq += discountedValue * discountedValue;
    }

    double averagePrice = sum / static_cast<double>(nbSim);
    double variance = (sumSq / static_cast<double>(nbSim)) - (averagePrice * averagePrice);
    return averagePrice;
}
