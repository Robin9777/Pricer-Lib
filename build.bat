@echo off
setlocal enabledelayedexpansion

:: Get the directory of the script
set ROOT=%~dp0
cd /d "%ROOT%"

:: Output files
set OUT_TESTS=tests\run_tests.exe
set OUT_ND=tests\run_tests_nd.exe
set OUT_QMC=tests\run_tests_qmc.exe
set OUT_VR=tests\run_tests_varred_compare.exe
set OUT_AT=tests\run_tests_antithetic.exe
set OUT_FC=tests\run_tests_full_comparison.exe
set OUT_DEMO=demo\run_demo.exe
set OUT_STUDY=demo\run_study.exe

:: Clean option
if "%1"=="clean" (
    del /f /q "%OUT_TESTS%" "%OUT_ND%" "%OUT_QMC%" "%OUT_VR%" "%OUT_AT%" "%OUT_FC%" "%OUT_DEMO%" "%OUT_STUDY%" 2>nul
    echo Cleaned.
    exit /b 0
)

:: List of library source files
set LIB_SOURCES= ^
    RandomGenerator\RandomGenerator.cpp ^
    RandomGenerator\UniformGenerator.cpp ^
    RandomGenerator\PseudoGenerator.cpp ^
    RandomGenerator\LinearCongruential.cpp ^
    RandomGenerator\EcuyerCombined.cpp ^
    RandomGenerator\ContinuousGenerator.cpp ^
    RandomGenerator\Normal.cpp ^
    RandomGenerator\Exponential.cpp ^
    RandomGenerator\DiscreteGenerator.cpp ^
    RandomGenerator\Bernoulli.cpp ^
    RandomGenerator\Binomial.cpp ^
    RandomGenerator\Poisson.cpp ^
    RandomGenerator\FiniteSet.cpp ^
    RandomGenerator\HeadTail.cpp ^
    RandomGenerator\HaltonGenerator.cpp ^
    RandomGenerator\AntitheticNormal.cpp ^
    SDE\SinglePath.cpp ^
    SDE\RandomProcess.cpp ^
    SDE\BlackScholes1D.cpp ^
    SDE\BlackScholes2D.cpp ^
    SDE\BSEuler1D.cpp ^
    SDE\BSMilstein1D.cpp ^
    SDE\BSMilstein2D.cpp ^
    SDE\BrownianD1.cpp ^
    SDE\BrownianND.cpp ^
    SDE\Heston.cpp ^
    SDE\BlackScholesND.cpp ^
    SDE\BSMilsteinND.cpp ^
    Payoffs\PayOff.cpp ^
    Payoffs\EuropeanCallPayoff.cpp ^
    Payoffs\EuroCallBasketPayOff.cpp ^
    Pricer\MCPricer.cpp ^
    Pricer\EuropeanMCPricer.cpp ^
    Pricer\BermudanPricer.cpp ^
    Pricer\VarRedMCPricer.cpp ^
    Pricer\AntitheticMCPricer.cpp ^
    VarRed\BSClosedForm.cpp ^
    VarRed\BasketGeomControlVariate.cpp ^
    VarRed\ControlVariate.cpp

:: Compilation flags
set FLAGS=-std=c++17 -O2 -I.

echo Compiling tests\run_tests...
g++ %FLAGS% %LIB_SOURCES% tests\test_main.cpp -o "%OUT_TESTS%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo Compiling tests\run_tests_nd...
g++ %FLAGS% %LIB_SOURCES% tests\test_milstein_nd.cpp -o "%OUT_ND%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo.
echo --- Running tests\run_tests ---
"%OUT_TESTS%"

echo.
echo --- Running tests\run_tests_nd ---
"%OUT_ND%"

echo.
echo Compiling tests\run_tests_qmc...
g++ %FLAGS% %LIB_SOURCES% tests\test_qmc.cpp -o "%OUT_QMC%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo --- Running tests\run_tests_qmc ---
"%OUT_QMC%"

echo.
echo Compiling tests\run_tests_varred_compare...
g++ %FLAGS% %LIB_SOURCES% tests\test_varred_compare.cpp -o "%OUT_VR%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo --- Running tests\run_tests_varred_compare ---
"%OUT_VR%"

echo.
echo Compiling tests\run_tests_antithetic...
g++ %FLAGS% %LIB_SOURCES% tests\test_antithetic.cpp -o "%OUT_AT%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo --- Running tests\run_tests_antithetic ---
"%OUT_AT%"

echo.
echo Compiling tests\run_tests_full_comparison...
g++ %FLAGS% %LIB_SOURCES% tests\test_full_comparison.cpp -o "%OUT_FC%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo --- Running tests\run_tests_full_comparison ---
"%OUT_FC%"

echo.
echo Compiling demo\run_demo...
g++ %FLAGS% %LIB_SOURCES% demo\demo.cpp -o "%OUT_DEMO%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo Demo ready. Run interactively: demo\run_demo.exe
echo            or with config:     demo\run_demo.exe demo\example_config.json

echo.
echo Compiling demo\run_study...
g++ %FLAGS% %LIB_SOURCES% demo\study.cpp -o "%OUT_STUDY%"
if %ERRORLEVEL% neq 0 exit /b %ERRORLEVEL%

echo Study ready. Run with: demo\run_study.exe demo\study_config.json
echo.
echo Build complete.
