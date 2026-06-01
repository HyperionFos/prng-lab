@echo off
g++ -O2 -std=c++17 -Wall main.cpp -o prng-lab.exe
if %errorlevel% neq 0 (
    echo.
    echo [BUILD FAILED]
    exit /b %errorlevel%
)
echo.
echo [BUILD OK] -^> prng-lab.exe