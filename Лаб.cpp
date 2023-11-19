#include <iostream>
#include <cmath>
#include <fstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Метод итераций 
double fixedPointIteration(double M, double e, double initial_guess, double tolerance, int max_iterations) {
    double E = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double next_E = M + e * sin(E);
        if (std::abs(next_E - E) < tolerance) {
            return next_E;
        }
        E = next_E;
    }
    return -1; // Возвращаем -1 в случае неудачи или превышения максимального числа итераций
}

// Метод половинного деления
double bisectionMethod(double M, double e, double a, double b, double tolerance, int max_iterations) {
    double E_a = a;
    double E_b = b;
    double f_a = E_a - e * sin(E_a) - M;
    for (int i = 0; i < max_iterations; ++i) {
        double E_mid = (E_a + E_b) / 2;
        double f_mid = E_mid - e * sin(E_mid) - M;
        if (std::abs(f_mid) < tolerance * (E_b - E_a) / 2 < tolerance) {
            return E_mid;
        }
        if ((f_mid > 0 && f_a > 0) || (f_mid < 0 && f_a < 0)) {
            E_a = E_mid;
            f_a = f_mid;
        }
        else {
            E_b = E_mid;
        }
    }
    return -1; 
}

// Метод золотого сечения
double goldenSectionMethod(double M, double e, double a, double b, double tolerance, int max_iterations) {
    const double golden_ratio = (1 + sqrt(5)) / 2;
    double x1 = b - (b - a) / golden_ratio;
    double x2 = a + (b - a) / golden_ratio;
    double f_x1 = x1 - e * sin(x1) - M;
    for (int i = 0; i < max_iterations; ++i) {
        double f_x2 = x2 - e * sin(x2) - M;
        if (std::abs(b - a) < tolerance) {
            return (a + b) / 2;
        }
        if (f_x1 < f_x2) {
            b = x2;
            x2 = x1;
            x1 = b - (b - a) / golden_ratio;
            f_x2 = f_x1;
            f_x1 = x1 - e * sin(x1) - M;
        }
        else {
            a = x1;
            x1 = x2;
            x2 = a + (b - a) / golden_ratio;
            f_x1 = f_x2;
            f_x2 = x2 - e * sin(x2) - M;
        }
    }
    return -1; 
}


// Метод Ньютона 
double newtonMethod(double M, double e, double initial_guess, double tolerance, int max_iterations) {
    double E = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double f = E - e * sin(E) - M;
        double f_prime = 1 - e * cos(E);
        double delta = f / f_prime;
        E -= delta;
        if (std::abs(delta) < tolerance) {
            return E;
        }
    }
    return -1; 
}

double calculateMeanAnomaly(double time)
{
     return ((2 * M_PI) / 5554.14) * time;
    
}

double calculateTrueAnomaly(double time,double e) {

    double M = ((2 * M_PI) / 5554.14) * time;
    double curr_ecs_anomaly = bisectionMethod(M, 0.0001112, 0, 2 * M_PI, 1e-8, 1000);
    return acos((cos(curr_ecs_anomaly) - e) / (1 - e * cos(curr_ecs_anomaly)));

}

double ecs_anomaly(double time) {
    double M = ((2 * M_PI) / 5554.14) * time;
    return bisectionMethod(M, 0.0001112, 0, 2 * M_PI, 1e-8, 1000);
}
// Функция для создания файла с данными для построения графика
void createDataFile(double res, double e)
{
    std::ofstream dataFile("trueAnomaly.txt");

    if (dataFile.is_open())
    {
        // Задаем шаг времени и пределы для построения графика
        double step = 10.0;
        double startTime = 0.0;
        double endTime = 5573.0;

        // Цикл для записи значений времени и средней аномалии в файл
        for (double t = startTime; t <= endTime; t += step)
        {
            double Anomaly = calculateTrueAnomaly(t, e);
            dataFile << t << " " << Anomaly << std::endl;
        }

        dataFile.close();
        std::cout << "Файл данных  2 успешно создан!" << std::endl;
    }
    else
    {
        std::cout << "Не удалось создать файл данных!" << std::endl;
    }
}

int main() {
    double M = ((2 * M_PI) / 5554.14) * 788832032; // Значение средней аномалии на 19.11.2023 17:37
    double e = 0.0001112; // Эксцентриситет орбиты

    // Пример использования методов для вычисления эксцентрической аномалии
    double initial_guess = M; // Начальное предположение для всех методов
    double tolerance = 1e-8; // Точность
    int max_iterations = 1000; // Максимальное количество итераций

    double result_fixed_point = fixedPointIteration(M, e, initial_guess, tolerance, max_iterations);
    double result_bisection = bisectionMethod(M, e, 0, 2 * M_PI, tolerance, max_iterations);
    double result_golden_section = goldenSectionMethod(M, e, 0, 2 * M_PI, tolerance, max_iterations);
    double result_newton = newtonMethod(M, e, initial_guess, tolerance, max_iterations);

    // Вывод результатов
    std::cout << "Метод итераций: " << result_fixed_point << std::endl;
    std::cout << "Метод половинного деления: " << result_bisection << std::endl;
    std::cout << "Метод золотого сечения: " << result_golden_section << std::endl;
    std::cout << "Метод Ньютона: " << result_newton << std::endl;
    // createDataFile(result_fixed_point, e);
    return 0;
}