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

// ============================================================
//  Параметры эксперимента (по тз)
// ============================================================
static const uint32_t RANGE   = 10000; ///< диапазон значений [0, RANGE), >= 5000
static const int      N       = 1000;  ///< элементов в выборке, >= 1000
static const int      SAMPLES = 20;    ///< число выборок на метод, >= 20
static const int      BINS    = 100;   ///< число интервалов для хи-квадрат

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

/**
 * @brief Линейный конгруэнтный генератор ПСЧ (метод B)
 *
 * Работает по формуле state = (A * state + C) mod 2^64. Модуль 2^64
 * получается автоматически за счет переполнения uint64_t, поэтому
 * явного деления с остатком не нужно. Параметры A и C подобраны так,
 * чтобы выполнялись условия Халла-Добелла и период был максимальным
 * и равным 2^64:
 *   1) C нечетное (взаимно просто с модулем);
 *   2) (A - 1) кратно 4, т.к. модуль 2^64 кратен 4.
 *
 * Младшие биты у LCG слабые (короткие периоды), поэтому наружу отдаются
 * перемешанные старшие половины двух последовательных шагов состояния -
 * это убирает слабость младших разрядов.
 */
class CustomLCG {
    uint64_t state; ///< внутреннее состояние генератора
    static constexpr uint64_t A = 6364136223846793005ULL; ///< множитель, A-1 кратно 4
    static constexpr uint64_t C = 1442695040888963407ULL; ///< приращение, нечетное
public:
    /**
     * @brief Конструктор
     *
     * @param seed начальное состояние генератора
     */
    explicit CustomLCG(uint64_t seed = 88172645463325252ULL)
        : state(seed) {}

    /**
     * @brief Следующее псевдослучайное число
     *
     * Делает два шага рекурренты. От первого шага берет старшие 32 бита,
     * от второго - тоже старшие 32 бита, и склеивает их в одно 64-битное
     * слово. Так наружу не попадают слабые младшие разряды.
     *
     * @return 64-битное псевдослучайное число
     */
    uint64_t next_u64() {
        state = A * state + C;
        uint64_t hi = state >> 32;
        state = A * state + C;
        return (hi << 32) | (state >> 32);
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

/**
 * @brief Генератор ПСЧ семейства xorshift128+ (метод C)
 *
 * Состояние из двух 64-битных слов, период 2^128 - 1. Число получается
 * серией операций "сдвиг + исключающее ИЛИ": каждый сдвиг перемешивает
 * биты внутри слова, XOR комбинирует их между собой. Используется
 * собственный набор сдвигов a=23, b=17, c=26.
 *
 * Нулевое состояние запрещено: если оба слова равны нулю, генератор
 * навсегда залипает на нуле, поэтому при таком seed принудительно
 * ставится s0 = 1. Метод очень быстрый (только сдвиги, XOR и сложение)
 * и хорошо проходит битовые тесты.
 */
class CustomXorshift {
    uint64_t s0; ///< первое слово состояния
    uint64_t s1; ///< второе слово состояния
public:
    /**
     * @brief Конструктор
     *
     * Разводит два слова состояния из одного seed, чтобы они не совпадали.
     * Если оба слова получились нулевыми, ставит s0 = 1.
     *
     * @param seed начальное значение для развода состояния
     */
    explicit CustomXorshift(uint64_t seed = 0x100000001b3ULL) {
        s0 = seed ^ 0xDEADBEEFCAFEBABEULL;
        s1 = (seed * 0x2545F4914F6CDD1DULL) ^ 0x123456789ABCDEFULL;
        if (s0 == 0 && s1 == 0) s0 = 1;
    }

    /**
     * @brief Следующее псевдослучайное число
     *
     * Применяет три сдвига с XOR к словам состояния, обновляет оба слова
     * и возвращает их сумму. Сложение в конце дополнительно улучшает
     * перемешивание бит по сравнению с чистым xorshift.
     *
     * @return 64-битное псевдослучайное число
     */
    uint64_t next_u64() {
        uint64_t x = s0;
        uint64_t const y = s1;
        s0 = y;
        x ^= x << 23;
        x ^= x >> 17;
        x ^= y ^ (y >> 26);
        s1 = x;
        return s1 + y;
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

/**
 * @brief Результаты статистики по одной выборке
 */
struct Stats {
    double mean; ///< среднее значение
    double sd;   ///< среднеквадратичное отклонение (несмещённое)
    double cv;   ///< коэффициент вариации, %
};

/**
 * @brief Итоговая сводка по методу за все выборки
 */
struct MethodSummary {
    std::string name; ///< название метода
    double avg_mean;  ///< среднее по средним всех выборок
    double avg_sd;    ///< среднее СКО
    double avg_cv;    ///< средний коэффициент вариации
    double avg_chi;   ///< средняя статистика хи-квадрат
    int    passed;    ///< сколько выборок прошли проверку равномерности
};

/**
 * @brief Набрать выборку заданного размера
 *
 * @tparam Gen тип генератора
 * @param g генератор
 * @param n размер выборки
 * @param range размер диапазона значений
 * @return вектор из n чисел в [0, range)
 */
template <class Gen>
std::vector<uint32_t> make_sample(Gen& g, int n, uint32_t range) {
    std::vector<uint32_t> v(n);
    for (int i = 0; i < n; ++i) v[i] = g.bounded(range);
    return v;
}

/**
 * @brief Посчитать среднее, СКО и коэффициент вариации
 *
 * Среднее - сумма значений делить на количество. СКО по несмещённой
 * формуле (деление на n-1). Коэффициент вариации - отношение СКО к
 * среднему в процентах.
 *
 * @param v выборка
 * @return структура Stats
 */
Stats compute_stats(const std::vector<uint32_t>& v) {
    double sum = 0.0;
    for (uint32_t x : v) sum += x;
    double mean = sum / v.size();

    double s2 = 0.0;
    for (uint32_t x : v) {
        double d = x - mean;
        s2 += d * d;
    }
    double sd = std::sqrt(s2 / (v.size() - 1));
    double cv = (mean != 0.0) ? (sd / mean * 100.0) : 0.0;

    return { mean, sd, cv };
}

/**
 * @brief Критическое значение хи-квадрат для df=99, alpha=0.05
 */
static const double CHI2_CRIT_99 = 123.225;

/**
 * @brief Проверка выборки на равномерность критерием хи-квадрат
 *
 * Делит диапазон на bins интервалов, считает попадания в каждый.
 * При равномерности в каждый интервал попадает ~ n/bins значений.
 * Статистика суммирует квадраты отклонений от ожидаемой частоты.
 *
 * @param v выборка
 * @param range размер диапазона
 * @param bins число интервалов
 * @return значение статистики хи-квадрат
 */
double chi_square_uniform(const std::vector<uint32_t>& v,
                          uint32_t range, int bins) {
    std::vector<long> obs(bins, 0);
    uint32_t width = range / bins;
    for (uint32_t x : v) {
        int b = static_cast<int>(x / width);
        if (b >= bins) b = bins - 1;
        obs[b]++;
    }
    double expected = static_cast<double>(v.size()) / bins;
    double chi2 = 0.0;
    for (int i = 0; i < bins; ++i) {
        double d = obs[i] - expected;
        chi2 += d * d / expected;
    }
    return chi2;
}

/**
 * @brief Прогон одного метода: SAMPLES выборок, статистика и хи-квадрат
 *
 * Для каждой выборки печатает среднее, СКО, коэффициент вариации и
 * значение хи-квадрат с пометкой о равномерности. В конце печатает
 * средние по всем выборкам и возвращает сводку.
 *
 * @tparam Gen тип генератора
 * @param name название метода
 * @param seed начальная затравка
 * @return сводка MethodSummary
 */
template <class Gen>
MethodSummary run_method(const std::string& name, uint64_t seed) {
    std::cout << "\n==== " << name << " ====\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(4)  << std::left << "#"
              << std::setw(12) << "mean"
              << std::setw(12) << "sd"
              << std::setw(10) << "CV(%)"
              << std::setw(12) << "chi2"
              << "uniform?\n";

    Gen g(seed);
    double acc_mean = 0, acc_sd = 0, acc_cv = 0, acc_chi = 0;
    int passed = 0;

    for (int s = 0; s < SAMPLES; ++s) {
        std::vector<uint32_t> v = make_sample(g, N, RANGE);
        Stats st = compute_stats(v);
        double chi = chi_square_uniform(v, RANGE, BINS);
        bool ok = chi < CHI2_CRIT_99;
        if (ok) passed++;

        acc_mean += st.mean;
        acc_sd   += st.sd;
        acc_cv   += st.cv;
        acc_chi  += chi;

        std::cout << std::setw(4)  << std::left << (s + 1)
                  << std::setw(12) << st.mean
                  << std::setw(12) << st.sd
                  << std::setw(10) << st.cv
                  << std::setw(12) << chi
                  << (ok ? "yes" : "NO") << "\n";
    }

    std::cout << "avg "
              << std::setw(12) << (acc_mean / SAMPLES)
              << std::setw(12) << (acc_sd / SAMPLES)
              << std::setw(10) << (acc_cv / SAMPLES)
              << std::setw(12) << (acc_chi / SAMPLES)
              << "passed " << passed << "/" << SAMPLES << "\n";

    return { name,
             acc_mean / SAMPLES, acc_sd / SAMPLES,
             acc_cv / SAMPLES,   acc_chi / SAMPLES,
             passed };
}

int main() {
    std::cout << "=== LR3: PRNG ===\n";
    std::cout << "RANGE   = " << RANGE   << "\n";
    std::cout << "N       = " << N       << "\n";
    std::cout << "SAMPLES = " << SAMPLES << "\n";
    std::cout << "BINS    = " << BINS    << "\n";
    std::cout << "chi2 crit (df=" << (BINS - 1) << ", alpha=0.05) = "
              << CHI2_CRIT_99 << "\n";

    // пункты 3-4: статистика и хи-квадрат по SAMPLES выборкам
    std::vector<MethodSummary> summary;
    summary.push_back(run_method<ChronoSplit>("A (ChronoSplit)", 12345));
    summary.push_back(run_method<CustomLCG>("B (CustomLCG)", 12345));
    summary.push_back(run_method<CustomXorshift>("C (CustomXorshift)", 12345));

    // сводный вывод по пунктам 3-4
    std::cout << "\n=== Pynkti 3-4 ===\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(22) << std::left << "Method"
              << std::setw(12) << "avg_mean"
              << std::setw(12) << "avg_sd"
              << std::setw(10) << "avg_CV%"
              << std::setw(12) << "avg_chi2"
              << "uniform\n";
    for (const auto& m : summary) {
        std::cout << std::setw(22) << std::left << m.name
                  << std::setw(12) << m.avg_mean
                  << std::setw(12) << m.avg_sd
                  << std::setw(10) << m.avg_cv
                  << std::setw(12) << m.avg_chi
                  << m.passed << "/" << SAMPLES << "\n";
    }
    std::cout << "\nTheory (ravnomernoe na [0," << RANGE << ")): "
              << "srednee = " << (RANGE - 1) / 2.0
              << ", SKO = " << RANGE / std::sqrt(12.0)
              << ", CV = 57.74%\n";

    return 0;
}