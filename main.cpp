#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

#include <boost/multiprecision/gmp.hpp>


auto startTime = std::chrono::high_resolution_clock::now();


template<typename T>
boost::multiprecision::mpf_float factorial(const T &number) {
    boost::multiprecision::mpf_float num = 1;
    for (int i = 1; i <= number; i++) {
        num = num * i;
    }
    return num;
}

template<typename T, typename I>
boost::multiprecision::mpf_float power(const T &base, const I &expo) {
    boost::multiprecision::mpf_float num = base;
    for (int i = 1; i <= expo; i++) {
        num = base * num;
    }
    return num;
}

float time_elapsed() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime);
    return float(duration.count()) / 1000;
}

class chudnovsky_formula {
public:
    chudnovsky_formula() {
        boost::multiprecision::mpf_float::default_precision(10000000000);

        boost::multiprecision::mpf_float m;
        boost::multiprecision::mpf_float l;
        boost::multiprecision::mpf_float x;
        std::fstream outfile;
        outfile.open("hp_series.txt", std::ios::trunc | std::ios::out);

        std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
        outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);


        for (long unsigned i = 0; i <= 1000; ++i) {
//        m = (factorial(6 * i)) / (factorial(3 * i) * power(factorial(i), 3));
            m = calc_m_value(i);
            l = calc_l_value(i);
            x = calc_x_value(i);
            outfile << ((m * l) / x) << "\n";
            std::cout << i << "/" << 1000 << "     " << (i / 1000.000) * 100.00 << "%\n";
        }
        std::cout << time_elapsed() << std::endl;
    }

private:
    boost::multiprecision::mpf_float last_m = 0;
    boost::multiprecision::mpf_float last_l = 0;
    boost::multiprecision::mpf_float last_x = 0;

    boost::multiprecision::mpf_float calc_m_value(const long unsigned &i) {
        if (last_m == 0) {
            boost::multiprecision::mpf_float num = factorial(6 * i);
            boost::multiprecision::mpf_float den = factorial(3 * i);
            den = (den * power(factorial(i), 3));
            last_m = (num / den);
            return last_m;
        } else {
            last_m = last_m * (((12 * i - 2) * (12 * i - 6) * (12 * i - 10)) / power(i, 3));
            return last_m;
        }
    }

    boost::multiprecision::mpf_float calc_l_value(const long unsigned &i) {
        if (last_l == 0) {
            last_l = (545140134 * i) + 13591409;
        } else {
            last_l += 545140134;
        }
        return last_l;
    }

    boost::multiprecision::mpf_float calc_x_value(const long unsigned &i) {
        if (last_x == 0) {
            last_x = boost::multiprecision::pow(boost::multiprecision::mpf_float(-262537412640768000), i);
        } else {
            last_x *= -262537412640768000;
        }
        return last_x;
    }
};

void
rabinowitz_and_wagon(long long int start_iterations, long long int stop_iterations, const std::string &outfile_name) {
    boost::multiprecision::mpf_float num;
    boost::multiprecision::mpf_float den;
    boost::multiprecision::mpf_float::default_precision(1000000000);
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    std::fstream outfile;
    outfile.open(outfile_name, std::ios::out | std::ios::trunc);
    outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    for (long unsigned i = start_iterations; i <= stop_iterations; ++i) {
        std::cout << i << "/" << stop_iterations << '\n';
        if (i % 2 == 0) {
            num = 1;
        } else {
            num = -1;
        }
        den = ((2 * i) + 1);
        outfile << num / den << "\n";
    }
}


void spigot_algorithm(long long int start_iterations, long long int stop_iterations, const std::string &outfile_name) {
    boost::multiprecision::mpf_float num;
    boost::multiprecision::mpf_float den;
    boost::multiprecision::mpf_float::default_precision(1000000000);
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    std::fstream outfile;
    outfile.open(outfile_name, std::ios::out | std::ios::trunc);
    outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    for (long unsigned i = start_iterations; i <= stop_iterations; ++i) {
        std::cout << i << "/" << stop_iterations << '\n';
        num = ((boost::multiprecision::pow(factorial(i), 2)) *
               boost::multiprecision::pow(boost::multiprecision::mpf_float(2), (i + 1)));
        den = factorial((2 * i + 1));
        outfile << (num / den) << '\n';
    }
};


void chudnovsky_formula_no() {
    namespace mp = boost::multiprecision;
    boost::multiprecision::mpf_float num;
    boost::multiprecision::mpf_float den;
    boost::multiprecision::mpf_float::default_precision(10000000000);
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    std::fstream outfile;
    outfile.open("series.txt", std::ios::out | std::ios::trunc);
    outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    for (long unsigned i = 0; i < 10; ++i) {
        std::cout << i << "/" << 10000 << "     " << (i / 10000.00) * 100.00 << "%\n";
        num = pow(-1, i) * factorial(6 * i) * (545140134 * i + 13591409);
        den = (factorial(3 * i) * power(factorial(i), 3)) * mp::pow(mp::mpf_float(640320), mp::mpf_float(3 * i + 1.5));
        outfile << num / den << "\n";
    }
}


void newton_euler_formula() {
    boost::multiprecision::mpf_float pi = 0;
    boost::multiprecision::mpf_float num;
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);

    for (long unsigned i = 0; i <= 100; ++i) {
        boost::multiprecision::mpf_float den = 2 * i + 1;
        num = factorial(i);
        den = factorial(den);
        den = factorial(den);
        pi = pi + ((num / den) * 2);
        std::cout << pi << std::endl;
    }
}


void leibniz_series(long long int start_iterations, long long int stop_iterations) { // NOLINT
    boost::multiprecision::mpf_float::default_precision(1000);
    boost::multiprecision::mpf_float num = 4;
    boost::multiprecision::mpf_float pi = 0;
    boost::multiprecision::mpf_float dem = 1;
    std::fstream outfile;
    outfile.open("series.txt", std::ios::out | std::ios::trunc);
    outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    int counter = 0;
    for (long unsigned i = start_iterations; i <= stop_iterations; ++i) {
        counter++;
        boost::multiprecision::mpf_float x = num / dem;
        num *= -1;
        dem += 2;
        pi += x;
        if ((counter % 500000) == 0) {
            std::cout << pi << "\n";
        }
    }
    outfile.close();
}


boost::multiprecision::mpf_float sum_from_file(const std::string &file_name) {
    std::cout << "calculating_sum\n";
    boost::multiprecision::mpf_float::default_precision(1000);
    boost::multiprecision::mpf_float running_sum;
    std::fstream file_to_sum;
    file_to_sum.open(file_name, std::ios::in);
    std::string current_number_s;
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    int counter = 0;
    while (getline(file_to_sum, current_number_s)) {
        std::cout << counter << std::endl;
        boost::multiprecision::mpf_float x(current_number_s);
        running_sum += x;
        counter++;
    }
    return running_sum;
}


int main() { // NOLINT



//    rabinowitz_and_wagon(0, 100, "ss.txt");
//    std::cout << (4*sum_from_file("ss.txt"));

//    spigot_algorithm(0, 10, "spigot_algo.txt");
//    std::cout << sum_from_file("spigot_algo.txt");
//    newton_euler_formula();
//    chudnovsky_formula();
//    leibniz_series(1, 10000000);
//    boost::multiprecision::mpf_float::default_precision(1000000);
//    chudnovsky_formula n;
    auto pi = sum_from_file("hp_series.txt");

//    std::cout << pi << std::endl;
//    pi = boost::multiprecision::pow(pi, -1);
//    boost::multiprecision::mpf_float const_d = (426880 * boost::multiprecision::pow(boost::multiprecision::mpf_float(10005), 0.5));
//    pi = pi*const_d;
//    pi = ((426880 * pow(boost::multiprecision::mpf_float(10005), 0.5))/pi);
//    std::cout << pi << std::endl;

//    std::cout << time_elapsed() << std::endl;
//    return 0;
}
