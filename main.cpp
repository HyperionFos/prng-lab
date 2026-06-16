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
#include <fstream>
#include <filesystem>

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

// ============================================================
//  Инфраструктура для тестов NIST (пункт 5)
// ============================================================

/**
 * @brief Превратить выборку 64-битных слов в поток бит
 *
 * Берёт n_words чисел от генератора и раскладывает каждое на 64 бита
 * (старший бит первым). Тесты NIST работают именно с битовым потоком.
 *
 * @tparam Gen тип генератора
 * @param g генератор
 * @param n_words сколько 64-битных слов сгенерировать
 * @return вектор бит (0 и 1) длиной n_words * 64
 */
template <class Gen>
std::vector<int> make_bits(Gen& g, int n_words) {
    std::vector<int> bits;
    bits.reserve(static_cast<size_t>(n_words) * 64);
    for (int i = 0; i < n_words; ++i) {
        uint64_t w = g.next_u64();
        for (int b = 63; b >= 0; --b)
            bits.push_back(static_cast<int>((w >> b) & 1ULL));
    }
    return bits;
}

double igamc(double a, double x); // предварительное объявление

/**
 * @brief Неполная гамма-функция P(a, x) (нижняя), разложение в ряд
 *
 * Вспомогательная функция для вычисления p-value. Считается рядом
 * при x < a+1, иначе через igamc.
 *
 * @param a параметр
 * @param x аргумент
 * @return значение P(a, x)
 */
double igam(double a, double x) {
    if (x <= 0.0 || a <= 0.0) return 0.0;
    if (x > a + 1.0) return 1.0 - igamc(a, x);
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 200; ++n) {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if (std::fabs(del) < std::fabs(sum) * 1e-15) break;
    }
    return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

/**
 * @brief Неполная гамма-функция Q(a, x) (верхняя, дополнительная)
 *
 * Переводит статистику хи-квадрат в p-value. Считается непрерывной
 * дробью при x >= a+1, иначе через igam. Стандартная функция из
 * библиотеки тестов NIST.
 *
 * @param a параметр (половина числа степеней свободы)
 * @param x аргумент (половина статистики)
 * @return p-value в диапазоне [0, 1]
 */
double igamc(double a, double x) {
    if (x <= 0.0 || a <= 0.0) return 1.0;
    if (x < a + 1.0) return 1.0 - igam(a, x);
    double b = x + 1.0 - a;
    double c = 1e30;
    double d = 1.0 / b;
    double h = d;
    for (int i = 1; i < 200; ++i) {
        double an = -1.0 * i * (i - a);
        b += 2.0;
        d = an * d + b; if (std::fabs(d) < 1e-30) d = 1e-30;
        c = b + an / c;  if (std::fabs(c) < 1e-30) c = 1e-30;
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < 1e-15) break;
    }
    return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}

/**
 * @brief Уровень значимости для тестов NIST
 */
static const double ALPHA = 0.01;

/**
 * @brief Тест 1 - частотный (monobit)
 *
 * Проверяет, поровну ли в потоке нулей и единиц. Каждый бит переводится
 * в +1 или -1, всё суммируется; при балансе сумма близка к нулю.
 *
 * @param bits битовый поток
 * @return p-value (тест пройден, если p >= ALPHA)
 */
double test_frequency(const std::vector<int>& bits) {
    long S = 0;
    for (int b : bits) S += (b ? 1 : -1);
    double n = static_cast<double>(bits.size());
    double s_obs = std::fabs(static_cast<double>(S)) / std::sqrt(n);
    return std::erfc(s_obs / std::sqrt(2.0));
}

/**
 * @brief Тест 2 - тест серий (runs)
 *
 * Считает число серий (переключений между 0 и 1) и сравнивает с
 * ожидаемым для случайной последовательности. Предусловие: доля единиц
 * близка к 0.5, иначе тест неприменим (p = 0).
 *
 * @param bits битовый поток
 * @return p-value (тест пройден, если p >= ALPHA)
 */
double test_runs(const std::vector<int>& bits) {
    double n = static_cast<double>(bits.size());
    long ones = 0;
    for (int b : bits) ones += b;
    double pi = ones / n;
    if (std::fabs(pi - 0.5) >= 2.0 / std::sqrt(n)) return 0.0;
    long V = 1;
    for (size_t i = 1; i < bits.size(); ++i)
        if (bits[i] != bits[i - 1]) V++;
    double num = std::fabs(V - 2.0 * n * pi * (1.0 - pi));
    double den = 2.0 * std::sqrt(2.0 * n) * pi * (1.0 - pi);
    return std::erfc(num / den);
}

/**
 * @brief Тест 3 - блочный частотный (block frequency)
 *
 * Делит поток на блоки длиной M бит и проверяет баланс единиц внутри
 * каждого. Ловит локальные перекосы при общем балансе 50/50.
 *
 * @param bits битовый поток
 * @param M длина блока в битах
 * @return p-value (тест пройден, если p >= ALPHA)
 */
double test_block_frequency(const std::vector<int>& bits, int M = 128) {
    int N = static_cast<int>(bits.size()) / M;
    if (N == 0) return 0.0;
    double chi2 = 0.0;
    for (int i = 0; i < N; ++i) {
        long ones = 0;
        for (int j = 0; j < M; ++j) ones += bits[i * M + j];
        double pi = static_cast<double>(ones) / M;
        double d = pi - 0.5;
        chi2 += d * d;
    }
    chi2 *= 4.0 * M;
    return igamc(N / 2.0, chi2 / 2.0);
}

/**
 * @brief Тест 4 - тест на пары бит (serial)
 *
 * Считает частоты четырёх 2-битных шаблонов (00, 01, 10, 11). При
 * случайном потоке они равновероятны; перекос означает зависимость
 * между соседними битами. Статистика хи-квадрат с 3 степенями свободы.
 *
 * @param bits битовый поток
 * @return p-value (тест пройден, если p >= ALPHA)
 */
double test_serial(const std::vector<int>& bits) {
    long cnt[4] = {0, 0, 0, 0};
    for (size_t i = 0; i + 1 < bits.size(); ++i) {
        int pat = (bits[i] << 1) | bits[i + 1];
        cnt[pat]++;
    }
    double total = static_cast<double>(bits.size() - 1);
    double expected = total / 4.0;
    double chi2 = 0.0;
    for (int k = 0; k < 4; ++k) {
        double d = cnt[k] - expected;
        chi2 += d * d / expected;
    }
    return igamc(3.0 / 2.0, chi2 / 2.0);
}

/**
 * @brief Тест 5 - накопленные суммы (cumulative sums)
 *
 * Биты переводятся в шаги +1/-1, строится нарастающая сумма (случайное
 * блуждание). Тест смотрит на максимальное удаление суммы от нуля -
 * большой размах означает неслучайность. P-value через нормальное
 * распределение.
 *
 * @param bits битовый поток
 * @return p-value (тест пройден, если p >= ALPHA)
 */
double test_cumulative_sums(const std::vector<int>& bits) {
    double n = static_cast<double>(bits.size());
    long S = 0, z = 0;
    for (int b : bits) {
        S += (b ? 1 : -1);
        if (std::labs(S) > z) z = std::labs(S);
    }
    if (z == 0) return 1.0;
    auto Phi = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };
    double sq = std::sqrt(n);
    double sum1 = 0.0, sum2 = 0.0;
    int start1 = static_cast<int>((-n / z + 1.0) / 4.0);
    int end1   = static_cast<int>(( n / z - 1.0) / 4.0);
    for (int k = start1; k <= end1; ++k) {
        sum1 += Phi(((4 * k + 1) * z) / sq);
        sum1 -= Phi(((4 * k - 1) * z) / sq);
    }
    int start2 = static_cast<int>((-n / z - 3.0) / 4.0);
    int end2   = static_cast<int>(( n / z - 1.0) / 4.0);
    for (int k = start2; k <= end2; ++k) {
        sum2 += Phi(((4 * k + 3) * z) / sq);
        sum2 -= Phi(((4 * k + 1) * z) / sq);
    }
    double p = 1.0 - sum1 + sum2;
    if (p < 0.0) p = 0.0;
    if (p > 1.0) p = 1.0;
    return p;
}

/**
 * @brief Прогон всех пяти тестов NIST для одного метода
 *
 * Генерирует битовый поток и прогоняет по нему пять тестов, печатая
 * для каждого p-value и вердикт PASS/FAIL (порог ALPHA).
 *
 * @tparam Gen тип генератора
 * @param name название метода
 * @param seed начальная затравка
 * @param n_words сколько 64-битных слов сгенерировать
 */
template <class Gen>
void run_nist_tests(const std::string& name, uint64_t seed, int n_words = 2000) {
    Gen g(seed);
    std::vector<int> bits = make_bits(g, n_words);

    double p1 = test_frequency(bits);
    double p2 = test_runs(bits);
    double p3 = test_block_frequency(bits, 128);
    double p4 = test_serial(bits);
    double p5 = test_cumulative_sums(bits);

    std::cout << "\n==== " << name << "  (" << bits.size() << " bits) ====\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  1. Frequency (monobit)  p=" << p1
              << "  " << (p1 >= ALPHA ? "PASS" : "FAIL") << "\n";
    std::cout << "  2. Runs                 p=" << p2
              << "  " << (p2 >= ALPHA ? "PASS" : "FAIL") << "\n";
    std::cout << "  3. Block frequency      p=" << p3
              << "  " << (p3 >= ALPHA ? "PASS" : "FAIL") << "\n";
    std::cout << "  4. Serial (2-bit)       p=" << p4
              << "  " << (p4 >= ALPHA ? "PASS" : "FAIL") << "\n";
    std::cout << "  5. Cumulative sums      p=" << p5
              << "  " << (p5 >= ALPHA ? "PASS" : "FAIL") << "\n";
}

/**
 * @brief Замерить среднее время генерации заданного числа элементов
 *
 * Прогоняет генерацию count чисел repeats раз и возвращает среднее
 * время одного прогона в миллисекундах. Сгенерированные числа
 * суммируются в sink, чтобы компилятор не выбросил цикл генерации
 * как не имеющий эффекта (защита от чрезмерной оптимизации).
 *
 * @tparam Gen тип генератора
 * @param seed начальная затравка
 * @param count сколько чисел генерировать за прогон
 * @param repeats сколько прогонов усреднить
 * @return среднее время прогона в миллисекундах
 */
template <class Gen>
double measure_ms(uint64_t seed, long count, int repeats) {
    volatile uint64_t sink = 0;
    double total_ms = 0.0;
    for (int r = 0; r < repeats; ++r) {
        Gen g(seed + r);
        auto t0 = std::chrono::high_resolution_clock::now();
        uint64_t acc = 0;
        for (long i = 0; i < count; ++i) acc += g.next_u64();
        auto t1 = std::chrono::high_resolution_clock::now();
        sink = sink + acc;
        total_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
    (void)sink;
    return total_ms / repeats;
}

/**
 * @brief Замер скорости для стандартного генератора mt19937_64
 *
 * Аналог measure_ms для стандартного генератора языка (Вихрь
 * Мерсенна) - эталон для сравнения.
 *
 * @param seed начальная затравка
 * @param count сколько чисел генерировать
 * @param repeats сколько прогонов усреднить
 * @return среднее время прогона в миллисекундах
 */
double measure_mt19937(uint64_t seed, long count, int repeats) {
    volatile uint64_t sink = 0;
    double total_ms = 0.0;
    for (int r = 0; r < repeats; ++r) {
        std::mt19937_64 g(seed + r);
        auto t0 = std::chrono::high_resolution_clock::now();
        uint64_t acc = 0;
        for (long i = 0; i < count; ++i) acc += g();
        auto t1 = std::chrono::high_resolution_clock::now();
        sink = sink + acc;
        total_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
    (void)sink;
    return total_ms / repeats;
}

/**
 * @brief Пункт 6: замер скорости генерации и запись в CSV
 *
 * Для объёмов 1000..1000000 элементов замеряет среднее время генерации
 * тремя собственными методами и стандартным mt19937_64. Печатает
 * таблицу в консоль и пишет timing.csv для построения графиков.
 *
 * @param repeats сколько прогонов усреднить на каждый замер
 */
void run_timing(int repeats = 100) {
    const long sizes[] = {1000, 10000, 100000, 1000000};
    const int n_sizes = 4;

    std::cout << "\n=== Пункт 6: скорость генерации (среднее по "
              << repeats << " прогонам, мс) ===\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << std::setw(12) << std::left << "N"
              << std::setw(14) << "A_ChronoSpl"
              << std::setw(14) << "B_LCG"
              << std::setw(14) << "C_Xorshift"
              << std::setw(14) << "std_mt19937" << "\n";

    std::filesystem::create_directories("data");
    std::ofstream csv("data/timing.csv");
    csv << "N,ChronoSplit,CustomLCG,CustomXorshift,mt19937\n";

    for (int i = 0; i < n_sizes; ++i) {
        long n = sizes[i];
        double tA = measure_ms<ChronoSplit>(12345, n, repeats);
        double tB = measure_ms<CustomLCG>(12345, n, repeats);
        double tC = measure_ms<CustomXorshift>(12345, n, repeats);
        double tM = measure_mt19937(12345, n, repeats);

        std::cout << std::setw(12) << std::left << n
                  << std::setw(14) << tA
                  << std::setw(14) << tB
                  << std::setw(14) << tC
                  << std::setw(14) << tM << "\n";

        csv << n << "," << tA << "," << tB << "," << tC << "," << tM << "\n";
    }
    csv.close();
    std::cout << "Saved data/timing.csv\n";
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
    system("chcp 65001 > nul");
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
    std::cout << "\n=== Сводка по пунктам 3-4 ===\n";
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
    std::cout << "\nТеория (равномерное на [0," << RANGE << ")): "
              << "среднее=" << (RANGE - 1) / 2.0
              << ", СКО=" << RANGE / std::sqrt(12.0)
              << ", CV=57.74%\n";

    // пункт 5: тесты NIST
    std::cout << "\n=== Пункт 5: тесты NIST (alpha=" << ALPHA
              << ", p>=alpha => PASS) ===\n";
    run_nist_tests<ChronoSplit>("A (ChronoSplit)", 12345);
    run_nist_tests<CustomLCG>("B (CustomLCG)", 12345);
    run_nist_tests<CustomXorshift>("C (CustomXorshift)", 12345);

    run_timing(10);

    return 0;
}