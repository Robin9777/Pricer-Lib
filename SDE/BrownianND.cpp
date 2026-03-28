#include "pch.h"
#include "BrownianND.h"
#include "../SDE/SinglePath.h"
#include <cmath>
#include <stdexcept>

BrownianND::BrownianND(RandomGenerator* Gen, int dim, std::vector<std::vector<double>>* Corr) : RandomProcess(Gen, dim)
{
    this->CorrelationMatrix = Corr;
}

std::vector<std::vector<double>> BrownianND::cholesky(const std::vector<std::vector<double>>& matrix)
{
    int n = matrix.size();
    std::vector<std::vector<double>> lower(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            if (j == i) {
                for (int k = 0; k < j; k++)
                    sum += pow(lower[j][k], 2);
                lower[j][j] = sqrt(matrix[j][j] - sum);
            } else {
                for (int k = 0; k < j; k++)
                    sum += (lower[i][k] * lower[j][k]);
                if (lower[j][j] == 0) {
                    throw std::runtime_error("Matrix is not positive definite.");
                }
                lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
            }
        }
    }
    return lower;
}

void BrownianND::Simulate(double startTime, double endTime, size_t nbSteps)
{
    // Clean up existing paths
    for (auto p : this->Paths) {
        delete p;
    }
    this->Paths.clear();

    // Create new paths
    for (int d = 0; d < this->Dimension; d++) {
        this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));
    }

    double dt = (endTime - startTime) / nbSteps;
    double sqrt_dt = std::sqrt(dt);

    std::vector<std::vector<double>> L;
    bool isCorrelated = (this->CorrelationMatrix != nullptr && !this->CorrelationMatrix->empty());
    if (isCorrelated) {
        L = cholesky(*(this->CorrelationMatrix));
    }

    std::vector<double> current_values(this->Dimension, 0.0);
    for (int d = 0; d < this->Dimension; d++) {
        this->Paths[d]->InsertValue(current_values[d]);
    }

    for (size_t i = 1; i <= nbSteps; i++) {
        // Generate independent standard normal increments
        std::vector<double> Z(this->Dimension);
        for (int d = 0; d < this->Dimension; d++) {
            Z[d] = this->Generator->Generate() * sqrt_dt;
        }

        // Apply correlation if present
        std::vector<double> dW(this->Dimension, 0.0);
        if (isCorrelated) {
            for (int d = 0; d < this->Dimension; d++) {
                for (int k = 0; k <= d; k++) {
                    dW[d] += L[d][k] * Z[k];
                }
            }
        } else {
            dW = Z;
        }

        // Update paths
        for (int d = 0; d < this->Dimension; d++) {
            current_values[d] += dW[d];
            this->Paths[d]->InsertValue(current_values[d]);
        }
    }
}
