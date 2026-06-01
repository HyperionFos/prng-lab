/**
 * @file main.cpp
 * @brief ЛР3 - Генераторы псевдослучайных чисел
 *
 * Реализованы три собственных метода генерации ПСЧ, проверка выборок
 * на равномерность (критерий хи-квадрат) и случайность (тесты NIST),
 * а также сравнение скорости генерации со стандартным генератором.
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <random>
#include <string>

//  Параметры эксперимента
static const uint32_t RANGE   = 10000; ///< диапазон значений [0, RANGE), >= 5000
static const int      N       = 1000;  ///< элементов в выборке, >= 1000
static const int      SAMPLES = 20;    ///< число выборок на метод, >= 20
static const int      BINS    = 100;   ///< число интервалов для хи-квадрат

int main() {
    std::cout << "=== LR3: PRNG ===\n";
    std::cout << "RANGE   = " << RANGE   << "\n";
    std::cout << "N       = " << N       << "\n";
    std::cout << "SAMPLES = " << SAMPLES << "\n";
    std::cout << "BINS    = " << BINS    << "\n";
    std::cout << "\n(skeleton: generators not added yet)\n";
    return 0;
}