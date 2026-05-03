#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cctype>
#include <sys/stat.h>

#include "RandomGenerator/LinearCongruential.h"
#include "RandomGenerator/Normal.h"
#include "RandomGenerator/HaltonGenerator.h"
#include "RandomGenerator/AntitheticNormal.h"
#include "SDE/BSMilstein1D.h"
#include "SDE/BSMilsteinND.h"
#include "Payoffs/EuropeanCallPayoff.h"
#include "Payoffs/EuroCallBasketPayOff.h"
#include "Pricer/EuropeanMCPricer.h"
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
    size_t nbSteps           = 1;
    size_t seed              = 42;
};

struct StudyConfig {
    OptionConfig option;
    std::vector<size_t> nbSimValues;
    size_t nSamples     = 10000;
    std::string outDir  = "results";
};

// ---------------------------------------------------------------------------
// JSON helpers
// ---------------------------------------------------------------------------

static std::string jsonGetString(const std::string& b, const std::string& key,
                                 const std::string& def = "") {
    std::string s = "\"" + key + "\"";
    auto pos = b.find(s);
    if (pos == std::string::npos) return def;
    pos = b.find(':', pos + s.size()); if (pos == std::string::npos) return def;
    ++pos;
    while (pos < b.size() && std::isspace((unsigned char)b[pos])) ++pos;
    if (pos >= b.size() || b[pos] != '"') return def;
    ++pos;
    auto end = b.find('"', pos);
    if (end == std::string::npos) return def;
    return b.substr(pos, end - pos);
}

static double jsonGetDouble(const std::string& b, const std::string& key, double def = 0.0) {
    std::string s = "\"" + key + "\"";
    auto pos = b.find(s);
    if (pos == std::string::npos) return def;
    pos = b.find(':', pos + s.size()); if (pos == std::string::npos) return def;
    ++pos;
    while (pos < b.size() && std::isspace((unsigned char)b[pos])) ++pos;
    try { return std::stod(b.substr(pos)); } catch (...) { return def; }
}

static std::vector<double> jsonGetArray(const std::string& b, const std::string& key) {
    std::vector<double> result;
    std::string s = "\"" + key + "\"";
    auto pos = b.find(s);
    if (pos == std::string::npos) return result;
    pos = b.find('[', pos + s.size()); if (pos == std::string::npos) return result;
    auto end = b.find(']', pos);       if (end == std::string::npos) return result;
    std::istringstream ss(b.substr(pos + 1, end - pos - 1));
    std::string token;
    while (std::getline(ss, token, ',')) {
        size_t i = token.find_first_not_of(" \t\n\r");
        if (i == std::string::npos) continue;
        try { result.push_back(std::stod(token.substr(i))); } catch (...) {}
    }
    return result;
}

static std::string extractBlock(const std::string& content, const std::string& key) {
    std::string s = "\"" + key + "\"";
    auto pos = content.find(s);
    if (pos == std::string::npos) return "";
    pos = content.find('{', pos + s.size());
    if (pos == std::string::npos) return "";
    int depth = 0; size_t end = pos;
    for (size_t i = pos; i < content.size(); ++i) {
        if      (content[i] == '{') ++depth;
        else if (content[i] == '}') { --depth; if (!depth) { end = i; break; } }
    }
    return content.substr(pos, end - pos + 1);
}

// ---------------------------------------------------------------------------
// Parse study config
// ---------------------------------------------------------------------------

static StudyConfig parseStudy(const std::string& filename) {
    std::ifstream f(filename);
    if (!f) throw std::runtime_error("Cannot open: " + filename);
    std::string content((std::istreambuf_iterator<char>(f)),
                         std::istreambuf_iterator<char>());

    StudyConfig study;

    std::string optBlock = extractBlock(content, "option");
    if (optBlock.empty()) throw std::runtime_error("Missing 'option' block.");

    auto& opt = study.option;
    opt.type     = jsonGetString(optBlock, "type",     "european");
    opt.spots    = jsonGetArray (optBlock, "spots");
    opt.vols     = jsonGetArray (optBlock, "vols");
    opt.weights  = jsonGetArray (optBlock, "weights");
    opt.rho      = jsonGetDouble(optBlock, "rho",      0.0);
    opt.strike   = jsonGetDouble(optBlock, "strike",   100.0);
    opt.rate     = jsonGetDouble(optBlock, "rate",     0.05);
    opt.maturity = jsonGetDouble(optBlock, "maturity", 1.0);
    opt.nbSteps  = (size_t)jsonGetDouble(optBlock, "nbSteps", 1.0);
    opt.seed     = (size_t)jsonGetDouble(optBlock, "seed",    42.0);

    if (opt.spots.empty() || opt.vols.empty())
        throw std::runtime_error("'spots' and 'vols' required.");
    int N = (int)opt.spots.size();
    if ((int)opt.vols.size() != N)
        throw std::runtime_error("'spots' and 'vols' must have the same length.");
    if (opt.weights.empty())
        opt.weights.assign(N, 1.0 / N);

    study.nbSimValues = [&]() -> std::vector<size_t> {
        auto v = jsonGetArray(content, "nbSim_values");
        std::vector<size_t> r;
        for (double x : v) r.push_back((size_t)x);
        return r;
    }();
    if (study.nbSimValues.empty())
        throw std::runtime_error("'nbSim_values' required.");

    study.nSamples = (size_t)jsonGetDouble(content, "nSamples", 10000.0);
    study.outDir   = jsonGetString(content, "output_dir", "results");

    return study;
}

// ---------------------------------------------------------------------------
// Price one method at a given nbSim
// ---------------------------------------------------------------------------

static PriceResult priceWith(const OptionConfig& cfg,
                              const std::string& method, size_t nbSim,
                              const std::vector<std::vector<double>>& corr) {
    int N = (int)cfg.spots.size();
    static const size_t A = 16807, C = 0, M = 2147483647;

    if (method == "qmc") {
        size_t nBases = cfg.nbSteps * (size_t)N * 2;
        auto bases = HaltonGenerator::firstNPrimes(nBases);
        HaltonGenerator ugen(bases);
        Normal normGen(0.0, 1.0, ugen);
        if (N == 1) {
            BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity, nbSim, cfg.nbSteps);
            return pricer.Price();
        }
        BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
        EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
        EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity, nbSim, cfg.nbSteps);
        return pricer.Price();
    }

    if (method == "antithetic") {
        LinearCongruential ugen(cfg.seed, A, C, M);
        Normal normRaw(0.0, 1.0, ugen);
        AntitheticNormal anti(normRaw);
        if (N == 1) {
            BSMilstein1D proc(&anti, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            AntitheticMCPricer pricer(&proc, &payoff, anti, cfg.rate, cfg.maturity,
                                     nbSim, cfg.nbSteps);
            return pricer.Price();
        }
        BSMilsteinND proc(&anti, cfg.spots, cfg.rate, cfg.vols, corr);
        EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
        AntitheticMCPricer pricer(&proc, &payoff, anti, cfg.rate, cfg.maturity,
                                  nbSim, cfg.nbSteps);
        return pricer.Price();
    }

    if (method == "staticcv") {
        LinearCongruential ugen(cfg.seed, A, C, M);
        Normal normGen(0.0, 1.0, ugen);
        if (N == 1) {
            BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            BasketGeomControlVariate cv({1.0}, cfg.spots, {{1.0}},
                                       cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
            VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity,
                                  nbSim, cfg.nbSteps);
            return pricer.Price();
        }
        BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
        EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
        BasketGeomControlVariate cv(cfg.weights, cfg.spots, corr,
                                   cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
        VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity,
                              nbSim, cfg.nbSteps);
        return pricer.Price();
    }

    if (method == "qmc_cv") {
        size_t nBases = cfg.nbSteps * (size_t)N * 2;
        auto bases = HaltonGenerator::firstNPrimes(nBases);
        HaltonGenerator ugen(bases);
        Normal normGen(0.0, 1.0, ugen);
        if (N == 1) {
            BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
            EuropeanCallPayoff payoff(cfg.strike);
            BasketGeomControlVariate cv({1.0}, cfg.spots, {{1.0}},
                                       cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
            VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity,
                                  nbSim, cfg.nbSteps);
            return pricer.Price();
        }
        BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
        EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
        BasketGeomControlVariate cv(cfg.weights, cfg.spots, corr,
                                   cfg.vols, cfg.maturity, cfg.strike, cfg.rate);
        VarRedMCPricer pricer(&proc, &payoff, &cv, cfg.rate, cfg.maturity,
                              nbSim, cfg.nbSteps);
        return pricer.Price();
    }

    // plain
    LinearCongruential ugen(cfg.seed, A, C, M);
    Normal normGen(0.0, 1.0, ugen);
    if (N == 1) {
        BSMilstein1D proc(&normGen, cfg.spots[0], cfg.rate, cfg.vols[0]);
        EuropeanCallPayoff payoff(cfg.strike);
        EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity, nbSim, cfg.nbSteps);
        return pricer.Price();
    }
    BSMilsteinND proc(&normGen, cfg.spots, cfg.rate, cfg.vols, corr);
    EuroCallBasketPayOff payoff(cfg.strike, cfg.weights);
    EuropeanMCPricer pricer(&proc, &payoff, cfg.rate, cfg.maturity, nbSim, cfg.nbSteps);
    return pricer.Price();
}

// ---------------------------------------------------------------------------
// Convergence study
// ---------------------------------------------------------------------------

static void runConvergence(const StudyConfig& study,
                           const std::vector<std::vector<double>>& corr) {
    std::string path = study.outDir + "/convergence.csv";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "method,nbSim,price,variance,ci\n";

    const std::vector<std::string> methods = {"plain", "antithetic", "staticcv", "qmc", "qmc_cv"};

    for (size_t nbSim : study.nbSimValues) {
        std::cout << "  nbSim=" << nbSim << "  ";
        for (const auto& method : methods) {
            auto r = priceWith(study.option, method, nbSim, corr);
            double se  = r.confidenceInterval / 1.96;
            double var = se * se * (double)nbSim;
            f << method << "," << nbSim << ","
              << r.price << "," << var << "," << r.confidenceInterval << "\n";
            std::cout << method << " ";
        }
        std::cout << "\n";
    }
}

// ---------------------------------------------------------------------------
// Terminal value collection (plain MC)
// ---------------------------------------------------------------------------

static void collectTerminalPlain(const StudyConfig& study,
                                  const std::vector<std::vector<double>>& corr) {
    std::string path = study.outDir + "/terminal_plain.csv";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "sT,payoff\n";

    int N = (int)study.option.spots.size();
    static const size_t A = 16807, C = 0, M = 2147483647;
    LinearCongruential ugen(study.option.seed, A, C, M);
    Normal normGen(0.0, 1.0, ugen);
    double discount = std::exp(-study.option.rate * study.option.maturity);
    double K = study.option.strike;

    if (N == 1) {
        BSMilstein1D proc(&normGen, study.option.spots[0],
                          study.option.rate, study.option.vols[0]);
        for (size_t i = 0; i < study.nSamples; ++i) {
            proc.Simulate(0.0, study.option.maturity, study.option.nbSteps);
            double sT = proc.GetPath(0)->GetAllValues().back();
            f << sT << "," << discount * std::max(sT - K, 0.0) << "\n";
        }
    } else {
        BSMilsteinND proc(&normGen, study.option.spots,
                          study.option.rate, study.option.vols, corr);
        for (size_t i = 0; i < study.nSamples; ++i) {
            proc.Simulate(0.0, study.option.maturity, study.option.nbSteps);
            double basket = 0.0;
            for (int d = 0; d < N; ++d)
                basket += study.option.weights[d] * proc.GetPath(d)->GetAllValues().back();
            f << basket << "," << discount * std::max(basket - K, 0.0) << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// Terminal value collection (antithetic pairs)
// ---------------------------------------------------------------------------

static void collectTerminalAntithetic(const StudyConfig& study,
                                       const std::vector<std::vector<double>>& corr) {
    std::string path = study.outDir + "/antithetic_pairs.csv";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write: " + path);
    f << "sT_z,payoff_z,sT_negz,payoff_negz\n";

    int N = (int)study.option.spots.size();
    static const size_t A = 16807, C = 0, M = 2147483647;
    LinearCongruential ugen(study.option.seed, A, C, M);
    Normal normRaw(0.0, 1.0, ugen);
    AntitheticNormal anti(normRaw);
    double discount = std::exp(-study.option.rate * study.option.maturity);
    double K = study.option.strike;
    size_t nPairs = study.nSamples / 2;

    if (N == 1) {
        BSMilstein1D proc(&anti, study.option.spots[0],
                          study.option.rate, study.option.vols[0]);
        for (size_t i = 0; i < nPairs; ++i) {
            anti.ResetBuffer(); anti.SetAntithetic(false);
            proc.Simulate(0.0, study.option.maturity, study.option.nbSteps);
            double sT_z = proc.GetPath(0)->GetAllValues().back();

            anti.ResetIndex(); anti.SetAntithetic(true);
            proc.Simulate(0.0, study.option.maturity, study.option.nbSteps);
            double sT_negz = proc.GetPath(0)->GetAllValues().back();

            f << sT_z     << "," << discount * std::max(sT_z     - K, 0.0) << ","
              << sT_negz  << "," << discount * std::max(sT_negz  - K, 0.0) << "\n";
        }
    } else {
        BSMilsteinND proc(&anti, study.option.spots,
                          study.option.rate, study.option.vols, corr);
        for (size_t i = 0; i < nPairs; ++i) {
            anti.ResetBuffer(); anti.SetAntithetic(false);
            proc.Simulate(0.0, study.option.maturity, study.option.nbSteps);
            double basket_z = 0.0;
            for (int d = 0; d < N; ++d)
                basket_z += study.option.weights[d] * proc.GetPath(d)->GetAllValues().back();

            anti.ResetIndex(); anti.SetAntithetic(true);
            proc.Simulate(0.0, study.option.maturity, study.option.nbSteps);
            double basket_negz = 0.0;
            for (int d = 0; d < N; ++d)
                basket_negz += study.option.weights[d] * proc.GetPath(d)->GetAllValues().back();

            f << basket_z    << "," << discount * std::max(basket_z    - K, 0.0) << ","
              << basket_negz << "," << discount * std::max(basket_negz - K, 0.0) << "\n";
        }
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./demo/run_study demo/study_config.json\n";
        return 1;
    }

    StudyConfig study;
    try {
        study = parseStudy(argv[1]);
    } catch (const std::exception& e) {
        std::cerr << "Config error: " << e.what() << "\n";
        return 1;
    }

    mkdir(study.outDir.c_str(), 0755);

    int N = (int)study.option.spots.size();
    std::vector<std::vector<double>> corr(N, std::vector<double>(N));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            corr[i][j] = (i == j) ? 1.0 : study.option.rho;

    std::cout << "=== Study ===\n"
              << "Option : " << study.option.type
              << "  N=" << N
              << "  K=" << study.option.strike
              << "  r=" << study.option.rate
              << "  T=" << study.option.maturity << "\n"
              << "Output : " << study.outDir << "/\n\n";

    std::cout << "Running convergence study...\n";
    try { runConvergence(study, corr); }
    catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }

    std::cout << "\nCollecting terminal values (plain MC, " << study.nSamples << " paths)...\n";
    try { collectTerminalPlain(study, corr); }
    catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }

    std::cout << "Collecting antithetic pairs (" << study.nSamples / 2 << " pairs)...\n";
    try { collectTerminalAntithetic(study, corr); }
    catch (const std::exception& e) { std::cerr << e.what() << "\n"; return 1; }

    std::cout << "\nDone. Results saved to " << study.outDir << "/\n";
    std::cout << "  convergence.csv\n"
              << "  terminal_plain.csv\n"
              << "  antithetic_pairs.csv\n";
    return 0;
}
