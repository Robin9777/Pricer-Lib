#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cctype>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "RandomGenerator/HaltonGenerator.h"
#include "RandomGenerator/AntitheticNormal.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
#include "Pricer/BermudanPricer.h"
#include "Pricer/VarRedMCPricer.h"
#include "Pricer/AntitheticMCPricer.h"
#include "VarRed/BasketGeomControlVariate.h"

struct OptionConfig {
    std::string type         = "european";
    std::vector<double> spots;
    std::vector<double> vols;
    std::vector<double> weights;
    double rho               = 0.0;
    double strike            = 100.0;
    double rate              = 0.05;
    double maturity          = 1.0;
    size_t nbSim             = 10000;
    size_t nbSteps           = 50;
    std::string varReduction = "none";
    std::vector<double> exerciseDates;
    size_t seed              = 42;
};

static double varFromCI(double ci, size_t n) {
    double se = ci / 1.96;
    return se * se * (double)n;
}

static void printRow(const std::string& label, const PriceResult& r, size_t n) {
    double var = varFromCI(r.confidenceInterval, n);
    std::cout << std::fixed
              << "  " << std::left << std::setw(16) << label
              << "  price=" << std::setprecision(4) << r.price
              << "  var="   << std::setprecision(1) << std::setw(8) << var
              << "  std="   << std::setprecision(4) << std::sqrt(var)
              << "  CI=+/-" << r.confidenceInterval << "\n";
}

static void printComparison(const PriceResult& base, const PriceResult& vr,
                             const std::string& label, size_t n) {
    double red = (1.0 - varFromCI(vr.confidenceInterval, n) /
                        varFromCI(base.confidenceInterval, n)) * 100.0;

    std::cout << "\n  Results:\n";
    printRow("Plain MC (base)", base, n);
    printRow(label,             vr,   n);
    std::cout << std::fixed << std::setprecision(1)
              << "\n  Variance reduction : " << red << "%\n"
              << std::setprecision(4)
              << "  Price difference   : " << (vr.price - base.price)
              << "  (should be ~0 within noise)\n";
}

static PriceResult runBaseline(const OptionConfig& cfg,
                               const std::vector<std::vector<double>>& corr) {
    int N = (int)cfg.spots.size();
    static const size_t A = 16807, C = 0, M = 2147483647;
    LinearCongruential ugen(cfg.seed, A, C, M);
    Normal normGen(0.0, 1.0, ugen);

    if (N == 1) {
        BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
        EuropeanCallPayoff payoff(cfg.strike);
        if (cfg.type == "bermudan") {
            BermudanPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity,
                                  cfg.nbSim, cfg.nbSteps, cfg.exerciseDates);
            return pricer.Price();
        }
        EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity,
                                cfg.nbSim, cfg.nbSteps);
        return pricer.Price();
    }

    BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
    EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
    if (cfg.type == "bermudan") {
        BermudanPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity,
                              cfg.nbSim, cfg.nbSteps, cfg.exerciseDates);
        return pricer.Price();
    }
    EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity,
                            cfg.nbSim, cfg.nbSteps);
    return pricer.Price();
}

static void runConfig(const OptionConfig& cfg) {
    int N = (int)cfg.spots.size();
    static const size_t A = 16807, C = 0, M = 2147483647;

    std::vector<std::vector<double>> corr(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            corr[i][j] = (i == j) ? 1.0 : cfg.rho;

    if (cfg.varReduction == "none") {
        std::cout << "\n  Results:\n";
        printRow("Plain MC", runBaseline(cfg, corr), cfg.nbSim);
        return;
    }

    if (cfg.varReduction == "qmc") {
        size_t nsteps  = (cfg.type == "bermudan") ? cfg.nbSteps : 1;
        size_t nBases  = nsteps * (size_t)N * 2;
        auto bases     = HaltonGenerator::firstNPrimes(nBases);

        OptionConfig baseCfg = cfg;
        baseCfg.nbSteps = nsteps;
        auto rBase = runBaseline(baseCfg, corr);

        HaltonGenerator ugen(bases);
        Normal normGen(0.0, 1.0, ugen);
        PriceResult rVR;

        if (N == 1) {
            BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            if (cfg.type == "bermudan") {
                BermudanPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity,
                                     cfg.nbSim, nsteps, cfg.exerciseDates);
                rVR = pricer.Price();
            } else {
                EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity, cfg.nbSim, 1);
                rVR = pricer.Price();
            }
        } else {
            BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
            EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
            EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity, cfg.nbSim, 1);
            rVR = pricer.Price();
        }

        printComparison(rBase, rVR, "QMC (Halton)", cfg.nbSim);
        return;
    }

    auto rBase = runBaseline(cfg, corr);

    if (cfg.varReduction == "antithetic") {
        LinearCongruential ugen(cfg.seed, A, C, M);
        Normal normRaw(0.0, 1.0, ugen);
        AntitheticNormal anti(normRaw);
        PriceResult rVR;

        if (N == 1) {
            BSMilstein1D proc(&anti, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            AntitheticMCPricer pricer(&proc, &payoff, anti, cfg.rate, cfg.maturity,
                                     cfg.nbSim, cfg.nbSteps);
            rVR = pricer.Price();
        } else {
            BSMilsteinND proc(&anti, cfg.spots, cfg.rate, cfg.vols, corr);
            EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
            AntitheticMCPricer pricer(&proc, &payoff, anti, cfg.rate, cfg.maturity,
                                     cfg.nbSim, cfg.nbSteps);
            rVR = pricer.Price();
        }

        printComparison(rBase, rVR, "Antithetic", cfg.nbSim);
        return;
    }

    if (cfg.varReduction == "staticcv") {
        LinearCongruential ugen(cfg.seed, A, C, M);
        Normal normGen(0.0, 1.0, ugen);
        PriceResult rVR;

        if (N == 1) {
            BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            BasketGeomControlVariate cv({1.0}, cfg.spots, {{1.0}},
                                       cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
            VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity,
                                  cfg.nbSim, cfg.nbSteps);
            rVR = pricer.Price();
        } else {
            BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
            EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
            BasketGeomControlVariate cv(cfg.weights, cfg.spots, corr,
                                       cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
            VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity,
                                  cfg.nbSim, cfg.nbSteps);
            rVR = pricer.Price();
        }

        printComparison(rBase, rVR, "Static CV", cfg.nbSim);
        return;
    }

    if (cfg.varReduction == "qmc_cv") {
        size_t nsteps = 1;
        size_t nBases = nsteps * (size_t)N * 2;
        auto bases    = HaltonGenerator::firstNPrimes(nBases);

        OptionConfig baseCfg = cfg;
        baseCfg.nbSteps = nsteps;
        auto rBase = runBaseline(baseCfg, corr);

        HaltonGenerator ugen(bases);
        Normal normGen(0.0, 1.0, ugen);
        PriceResult rVR;

        if (N == 1) {
            BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            BasketGeomControlVariate cv({1.0}, cfg.spots, {{1.0}},
                                       cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
            VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity, cfg.nbSim, 1);
            rVR = pricer.Price();
        } else {
            BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
            EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
            BasketGeomControlVariate cv(cfg.weights, cfg.spots, corr,
                                       cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
            VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity, cfg.nbSim, 1);
            rVR = pricer.Price();
        }

        printComparison(rBase, rVR, "QMC + Static CV", cfg.nbSim);
        return;
    }
}

static double readDouble(const std::string& prompt, double lo = -1e18, double hi = 1e18) {
    while (true) {
        std::cout << prompt;
        std::string line;
        if (!std::getline(std::cin, line)) { std::cin.clear(); continue; }
        try {
            size_t pos;
            double val = std::stod(line, &pos);
            while (pos < line.size() && std::isspace((unsigned char)line[pos])) ++pos;
            if (pos == line.size() && val >= lo && val <= hi) return val;
        } catch (...) {}
        std::cout << "  Invalid.";
        if (lo > -1e17) std::cout << " Min: " << lo;
        if (hi <  1e17) std::cout << " Max: " << hi;
        std::cout << "\n";
    }
}

static size_t readSizeT(const std::string& prompt, size_t lo = 1, size_t hi = 10000000) {
    while (true) {
        std::cout << prompt;
        std::string line;
        if (!std::getline(std::cin, line)) { std::cin.clear(); continue; }
        try {
            size_t pos;
            long long val = std::stoll(line, &pos);
            while (pos < line.size() && std::isspace((unsigned char)line[pos])) ++pos;
            if (pos == line.size() && val >= (long long)lo && val <= (long long)hi)
                return (size_t)val;
        } catch (...) {}
        std::cout << "  Invalid. Expected integer in [" << lo << ", " << hi << "]\n";
    }
}

static int readChoice(const std::string& prompt, int lo, int hi) {
    while (true) {
        std::cout << prompt;
        std::string line;
        if (!std::getline(std::cin, line)) { std::cin.clear(); continue; }
        try {
            size_t pos;
            int val = std::stoi(line, &pos);
            while (pos < line.size() && std::isspace((unsigned char)line[pos])) ++pos;
            if (pos == line.size() && val >= lo && val <= hi) return val;
        } catch (...) {}
        std::cout << "  Enter a number between " << lo << " and " << hi << "\n";
    }
}

static OptionConfig collectInteractive() {
    OptionConfig cfg;

    std::cout << "\n--- Option type ---\n"
              << "  1. European call\n"
              << "  2. Bermudan call\n";
    cfg.type = (readChoice("> ", 1, 2) == 1) ? "european" : "bermudan";

    std::cout << "\n--- Number of assets ---\n";
    int N = (int)readSizeT("Number of assets (1 for vanilla, 2+ for basket): ", 1, 10);

    cfg.spots.resize(N);
    cfg.vols.resize(N);
    for (int i = 0; i < N; ++i) {
        std::cout << "\n--- Asset " << (i + 1) << " ---\n";
        cfg.spots[i] = readDouble("  Spot price S0: ", 0.001, 1e9);
        cfg.vols[i]  = readDouble("  Volatility sigma: ", 0.001, 10.0);
    }

    if (N > 1) {
        std::cout << "\n--- Basket weights (must sum to 1) ---\n";
        cfg.weights.resize(N);
        double wsum = 0.0;
        for (int i = 0; i < N - 1; ++i) {
            cfg.weights[i] = readDouble("  Weight for asset " + std::to_string(i + 1) + ": ",
                                        0.001, 1.0 - wsum - 0.001 * (N - 1 - i));
            wsum += cfg.weights[i];
        }
        cfg.weights[N - 1] = 1.0 - wsum;
        std::cout << "  Weight for asset " << N << " (auto): " << cfg.weights[N - 1] << "\n";
        cfg.rho = readDouble("\nCorrelation rho between assets [-1, 1]: ", -0.999, 0.999);
    } else {
        cfg.weights = {1.0};
    }

    std::cout << "\n--- Contract ---\n";
    cfg.strike   = readDouble("Strike K: ", 0.001, 1e9);
    cfg.rate     = readDouble("Rate r (e.g. 0.05): ", 0.0, 1.0);
    cfg.maturity = readDouble("Maturity T (years): ", 0.001, 50.0);

    if (cfg.type == "bermudan") {
        std::cout << "\n--- Exercise dates ---\n";
        int nd = (int)readSizeT("Number of exercise dates: ", 1, 100);
        cfg.exerciseDates.resize(nd);
        for (int i = 0; i < nd; ++i)
            cfg.exerciseDates[i] = readDouble("  Date " + std::to_string(i + 1) + ": ",
                                              0.0, cfg.maturity);
    }

    std::cout << "\n--- Simulation ---\n";
    cfg.nbSim = readSizeT("Number of paths: ", 1, 10000000);

    std::cout << "\n--- Variance reduction ---\n";
    if (cfg.type == "european") {
        std::cout << "  0. None\n"
                  << "  1. QMC (forces nsteps=1)\n"
                  << "  2. Antithetic\n"
                  << "  3. Static CV (geometric basket)\n"
                  << "  4. QMC + Static CV\n";
        int vr = readChoice("> ", 0, 4);
        if      (vr == 0) cfg.varReduction = "none";
        else if (vr == 1) { cfg.varReduction = "qmc";    cfg.nbSteps = 1; }
        else if (vr == 2)   cfg.varReduction = "antithetic";
        else if (vr == 3)   cfg.varReduction = "staticcv";
        else              { cfg.varReduction = "qmc_cv"; cfg.nbSteps = 1; }
    } else {
        std::cout << "  0. None\n"
                  << "  1. QMC\n"
                  << "  (antithetic / static CV not available for Bermudan)\n";
        cfg.varReduction = (readChoice("> ", 0, 1) == 0) ? "none" : "qmc";
    }

    if (cfg.varReduction != "qmc")
        cfg.nbSteps = readSizeT("Number of time steps: ", 1, 10000);
    else if (cfg.type == "bermudan")
        cfg.nbSteps = readSizeT("Number of time steps: ", 1, 100);

    return cfg;
}

static std::string jsonGetString(const std::string& block, const std::string& key,
                                 const std::string& def = "") {
    std::string search = "\"" + key + "\"";
    auto pos = block.find(search);
    if (pos == std::string::npos) return def;
    pos = block.find(':', pos + search.size());
    if (pos == std::string::npos) return def;
    ++pos;
    while (pos < block.size() && std::isspace((unsigned char)block[pos])) ++pos;
    if (pos >= block.size() || block[pos] != '"') return def;
    ++pos;
    auto end = block.find('"', pos);
    if (end == std::string::npos) return def;
    return block.substr(pos, end - pos);
}

static double jsonGetDouble(const std::string& block, const std::string& key, double def = 0.0) {
    std::string search = "\"" + key + "\"";
    auto pos = block.find(search);
    if (pos == std::string::npos) return def;
    pos = block.find(':', pos + search.size());
    if (pos == std::string::npos) return def;
    ++pos;
    while (pos < block.size() && std::isspace((unsigned char)block[pos])) ++pos;
    try { return std::stod(block.substr(pos)); } catch (...) { return def; }
}

static std::vector<double> jsonGetArray(const std::string& block, const std::string& key) {
    std::vector<double> result;
    std::string search = "\"" + key + "\"";
    auto pos = block.find(search);
    if (pos == std::string::npos) return result;
    pos = block.find('[', pos + search.size());
    if (pos == std::string::npos) return result;
    auto end = block.find(']', pos);
    if (end == std::string::npos) return result;
    std::istringstream ss(block.substr(pos + 1, end - pos - 1));
    std::string token;
    while (std::getline(ss, token, ',')) {
        size_t s = token.find_first_not_of(" \t\n\r");
        if (s == std::string::npos) continue;
        try { result.push_back(std::stod(token.substr(s))); } catch (...) {}
    }
    return result;
}

static std::vector<OptionConfig> parseJson(const std::string& filename) {
    std::ifstream f(filename);
    if (!f) throw std::runtime_error("Cannot open: " + filename);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    std::vector<OptionConfig> configs;
    size_t pos = 0;
    while (pos < content.size()) {
        auto start = content.find('{', pos);
        if (start == std::string::npos) break;
        int depth = 0;
        size_t end = start;
        for (size_t i = start; i < content.size(); ++i) {
            if      (content[i] == '{') ++depth;
            else if (content[i] == '}') { --depth; if (!depth) { end = i; break; } }
        }

        std::string block = content.substr(start, end - start + 1);
        OptionConfig cfg;
        cfg.type          = jsonGetString(block, "type",         "european");
        cfg.spots         = jsonGetArray (block, "spots");
        cfg.vols          = jsonGetArray (block, "vols");
        cfg.weights       = jsonGetArray (block, "weights");
        cfg.rho           = jsonGetDouble(block, "rho",           0.0);
        cfg.strike        = jsonGetDouble(block, "strike",         100.0);
        cfg.rate          = jsonGetDouble(block, "rate",           0.05);
        cfg.maturity      = jsonGetDouble(block, "maturity",       1.0);
        cfg.nbSim         = (size_t)jsonGetDouble(block, "nbSim",  10000.0);
        cfg.nbSteps       = (size_t)jsonGetDouble(block, "nbSteps", 50.0);
        cfg.varReduction  = jsonGetString(block, "varReduction",  "none");
        cfg.exerciseDates = jsonGetArray (block, "exerciseDates");
        cfg.seed          = (size_t)jsonGetDouble(block, "seed",   42.0);

        if (cfg.spots.empty() || cfg.vols.empty())
            throw std::runtime_error("Each option needs 'spots' and 'vols'.");
        if (cfg.spots.size() != cfg.vols.size())
            throw std::runtime_error("'spots' and 'vols' must be the same length.");

        int N = (int)cfg.spots.size();
        if (cfg.weights.empty())
            cfg.weights.assign(N, 1.0 / N);
        else if ((int)cfg.weights.size() != N)
            throw std::runtime_error("'weights' must be the same length as 'spots'.");

        if (cfg.type == "bermudan" && cfg.exerciseDates.empty())
            throw std::runtime_error("Bermudan needs 'exerciseDates'.");

        configs.push_back(cfg);
        pos = end + 1;
    }

    if (configs.empty())
        throw std::runtime_error("No options found in file.");
    return configs;
}

static void printSummary(const OptionConfig& cfg, int idx = -1) {
    if (idx >= 0) std::cout << "\n=== Option " << (idx + 1) << " ===\n";
    else          std::cout << "\n=== Pricing ===\n";

    int N = (int)cfg.spots.size();
    std::cout << "  Type    : " << cfg.type << "\n"
              << "  Assets  : " << N << "\n";
    for (int i = 0; i < N; ++i)
        std::cout << "  Asset " << (i+1) << " : S0=" << cfg.spots[i]
                  << "  sigma=" << cfg.vols[i] << "  w=" << cfg.weights[i] << "\n";
    if (N > 1)
        std::cout << "  Rho     : " << cfg.rho << "\n";
    std::cout << "  K=" << cfg.strike << "  r=" << cfg.rate << "  T=" << cfg.maturity << "\n"
              << "  Paths=" << cfg.nbSim << "  Steps=" << cfg.nbSteps
              << "  VarRed=" << cfg.varReduction << "\n";
    if (!cfg.exerciseDates.empty()) {
        std::cout << "  Dates  :";
        for (double d : cfg.exerciseDates) std::cout << " " << d;
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::cout << "=== MC Option Pricer ===\n";

    if (argc >= 2) {
        try {
            auto configs = parseJson(argv[1]);
            std::cout << "Loaded " << configs.size() << " option(s) from " << argv[1] << "\n";
            for (int i = 0; i < (int)configs.size(); ++i) {
                printSummary(configs[i], i);
                runConfig(configs[i]);
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            return 1;
        }
    } else {
        while (true) {
            OptionConfig cfg;
            try { cfg = collectInteractive(); }
            catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << "\n";
                return 1;
            }
            printSummary(cfg);
            runConfig(cfg);

            std::cout << "\nPrice another? (y/n): ";
            std::string ans;
            std::getline(std::cin, ans);
            if (ans.empty() || (ans[0] != 'y' && ans[0] != 'Y')) break;
        }
    }

    std::cout << "\nDone.\n";
    return 0;
}
