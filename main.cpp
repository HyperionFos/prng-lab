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
static const int N = 1000;  ///< элементов в выборке, >= 1000
static const int SAMPLES = 20;    ///< число выборок на метод, >= 20
static const int BINS = 100;   ///< число интервалов для хи-квадрат

/**
 * @brief Генератор ПСЧ на основе времени (метод A)
 *
 * Момент времени в наносекундах берется один раз при создании объекта
 * и используется как затравка (seed). Дальше состояние перемешивается
 * собственным микшером на основе splitmix64: время напрямую не годится,
 * т.к. при вызовах подряд таймер растет почти линейно и соседние числа
 * были бы сильно похожи. Микшер с умножениями и сдвигами убирает эту
 * линейность и дает хороший разброс бит.
 */
class ChronoSplit {
    uint64_t state; ///< внутреннее состояние генератора
public:
    /**
     * @brief Конструктор
     *
     * Если seed равен нулю, берем текущее время в наносекундах как затравку.
     * Любое ненулевое значение дает воспроизводимую последовательность.
     * Затравка умножается на большое простое число и складывается с
     * константой золотого сечения, чтобы развести стартовое состояние.
     *
     * @param seed начальное значение, 0 - взять время автоматически
     */
    explicit ChronoSplit(uint64_t seed = 0) {
        if (seed == 0) {
            auto now = std::chrono::high_resolution_clock::now()
                           .time_since_epoch();
            seed = static_cast<uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(now)
                    .count());
        }
        state = seed * 1099511628211ULL + 0x9E3779B97F4A7C15ULL;
    }
 
    /**
     * @brief Следующее псевдослучайное число
     *
     * Сдвигает состояние на шаг золотого сечения, затем дважды
     * перемешивает биты умножением на простые константы и сдвигами
     * вправо (лавинный эффект: один измененный бит состояния меняет
     * примерно половину бит результата).
     *
     * @return 64-битное псевдослучайное число
     */
    uint64_t next_u64() {
        state += 0x9E3779B97F4A7C15ULL;
        uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
 
    /**
     * @brief Число в заданном диапазоне
     *
     * @param range размер диапазона
     * @return число в полуинтервале [0, range)
     */
    uint32_t bounded(uint32_t range) {
        return static_cast<uint32_t>(next_u64() % range);
    }
};

int main() {
    std::cout << "=== LR3: PRNG ===\n";
    std::cout << "RANGE   = " << RANGE   << "\n";
    std::cout << "N       = " << N       << "\n";
    std::cout << "SAMPLES = " << SAMPLES << "\n";
    std::cout << "BINS    = " << BINS    << "\n";
    
    // демонстрация метода A
    ChronoSplit a(123456789ULL); // фиксированный seed -> воспроизводимо
    std::cout << "Method A (ChronoSplit), first 10 in [0, " << RANGE << "):\n  ";
    for (int i = 0; i < 10; ++i) std::cout << a.bounded(RANGE) << " ";
    std::cout << "\n";

    return 0;
}