"""
Построение графиков сравнения скорости генерации ПСЧ (пункт 6).
Читает data/timing.csv (создаётся программой main) и строит два графика
в папку data/:
  - линейный масштаб (timing_linear.png)
  - логарифмический масштаб по обеим осям (timing_log.png)
Запуск из корня проекта:  python scripts/plot_timing.py
"""
import csv
import os
import matplotlib.pyplot as plt

DATA_DIR = "data"
CSV_PATH = os.path.join(DATA_DIR, "timing.csv")

# чтение CSV
sizes = []
data = {}
with open(CSV_PATH, newline="") as f:
    reader = csv.reader(f)
    header = next(reader)
    methods = header[1:]            # имена методов
    for m in methods:
        data[m] = []
    for row in reader:
        sizes.append(int(row[0]))
        for i, m in enumerate(methods):
            data[m].append(float(row[i + 1]))

markers = ["o", "s", "^", "d"]

# --- график 1: линейный масштаб ---
plt.figure(figsize=(8, 5))
for i, m in enumerate(methods):
    plt.plot(sizes, data[m], marker=markers[i % len(markers)], label=m)
plt.xlabel("Объём выборки, элементов")
plt.ylabel("Время генерации, мс")
plt.title("Скорость генерации ПСЧ (линейный масштаб)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(DATA_DIR, "timing_linear.png"), dpi=120)
print("saved", os.path.join(DATA_DIR, "timing_linear.png"))

# --- график 2: логарифмический масштаб ---
plt.figure(figsize=(8, 5))
for i, m in enumerate(methods):
    plt.plot(sizes, data[m], marker=markers[i % len(markers)], label=m)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Объём выборки, элементов (log)")
plt.ylabel("Время генерации, мс (log)")
plt.title("Скорость генерации ПСЧ (логарифмический масштаб)")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(DATA_DIR, "timing_log.png"), dpi=120)
print("saved", os.path.join(DATA_DIR, "timing_log.png"))