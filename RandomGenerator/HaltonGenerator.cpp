#include "pch.h"
#include "HaltonGenerator.h"

HaltonGenerator::HaltonGenerator(size_t base)
    : bases({base}), counters({0}), currentDim(0)
{
}

HaltonGenerator::HaltonGenerator(const std::vector<size_t>& _bases)
    : bases(_bases), counters(_bases.size(), 0), currentDim(0)
{
}

double HaltonGenerator::radicalInverse(size_t n, size_t base) const
{
    double result = 0.0;
    double f = 1.0;
    while (n > 0) {
        f /= static_cast<double>(base);
        result += f * static_cast<double>(n % base);
        n /= base;
    }
    return result;
}

double HaltonGenerator::Generate()
{
    size_t dim = currentDim;
    ++counters[dim];
    double val = radicalInverse(counters[dim], bases[dim]);
    currentDim = (currentDim + 1) % bases.size();
    return val;
}

void HaltonGenerator::Reset()
{
    std::fill(counters.begin(), counters.end(), 0);
    currentDim = 0;
}

std::vector<size_t> HaltonGenerator::firstNPrimes(size_t n)
{
    std::vector<size_t> primes;
    primes.reserve(n);
    for (size_t candidate = 2; primes.size() < n; ++candidate) {
        bool isPrime = true;
        for (size_t p : primes) {
            if (p * p > candidate) break;
            if (candidate % p == 0) { isPrime = false; break; }
        }
        if (isPrime) primes.push_back(candidate);
    }
    return primes;
}
