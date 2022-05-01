#include <iostream>
#include <fstream>
#include <cmath>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <boost/math/special_functions/factorials.hpp>




template <typename T>
boost::multiprecision::mpf_float factorial(const T &number )
{
    boost::multiprecision::mpf_float num = 1;
    for (int i = 1; i <= number; i++){
        num = num * i;
    }
    return num;
}

template <typename T, typename I>
boost::multiprecision::mpf_float power(const T &base, const I &expo)
{
    boost::multiprecision::mpf_float num = base;
    for (int i = 1; i <= expo; i++){
        num = base * num;
    }
    return num;
}



void chudnovsky_formula(){
    boost::multiprecision::mpf_float  pi = 0;
    boost::multiprecision::mpf_float  m;
    boost::multiprecision::mpf_float  l;
    boost::multiprecision::mpf_float  x;
    boost::multiprecision::mpf_float  c = 426880*sqrt(100.024996875781);
    boost::multiprecision::mpf_float  num;
    boost::multiprecision::mpf_float  den;

    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);

    for(long unsigned i = 0; i <= 1000;  ++i){
        m = (factorial(6*i))/(factorial(3*i)*power(factorial(i), 3));
        l = (545140134*i + 13591409);
        x = (power(-262537412640768000, i));
        pi = pi + ((m*l)/x);
        std::cout << pi << std::endl;
    }
    pi = pow(640320, 1.5)/(pi*12);
    std::cout << pi << "\n";
}

void chudnovsky_formula_no(){
    namespace mp = boost::multiprecision;
    boost::multiprecision::mpf_float  pi = 0;
    boost::multiprecision::mpf_float  num;
    boost::multiprecision::mpf_float  den;
    boost::multiprecision::mpf_float::default_precision(50000);
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    std::fstream outfile;
    outfile.open("series.txt", std::ios::out | std::ios::trunc);
    outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    for(long unsigned i = 0; i < 10000;  ++i) {
        std::cout << i << "/" << 10000 << "     " << (i/10000.0)*100.00 << "\n";
        num = pow(-1, i) * factorial(6 * i) * (545140134 * i + 13591409);
        den = (factorial(3 * i) * power(factorial(i), 3) * mp::pow(mp::mpf_float(640320), (3 * i + 1.5)));
        outfile << num/den << '\n';
    }
//    pi *= 12;
//    pi = 1/pi;
//    std::cout << pi << "\n";
}


void newton_euler_formula(){
    boost::multiprecision::mpf_float  pi = 0;
    boost::multiprecision::mpf_float  num;
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);

    for(long unsigned i = 0; i <= 100; ++i){
        boost::multiprecision::mpf_float den = 2*i+1;
        num = factorial(i);
        den = factorial(den);
        den = factorial(den);
        pi  = pi + ((num/den)*2);
        std::cout << pi << std::endl;
    }
}


void leibniz_series(long long int start_iterations, long long int stop_iterations){ // NOLINT
    boost::multiprecision::mpf_float::default_precision(1000);
    boost::multiprecision::mpf_float  num = 4;
    boost::multiprecision::mpf_float  pi = 0;
    boost::multiprecision::mpf_float  dem = 1;
    std::fstream outfile;
    outfile.open("series.txt", std::ios::out | std::ios::trunc);
    outfile.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);
    int counter = 0;
    for(long unsigned i = start_iterations; i <= stop_iterations; ++i){
        counter ++;
        boost::multiprecision::mpf_float  x = num/dem;
        num *= -1;
        dem += 2;
        pi += x;
//        outfile << x << '\n';
        if(counter % 5000000000){
            std::cout << pi << "\n";
//            cout << counter << "/" << stop_iterations << "\n";
//            cout << float(float(counter)/float(stop_iterations)) << "\n";
        }
    }
    outfile.close();
}

boost::multiprecision::mpf_float sum_from_file(const std::string& file_name){
    std::cout << "calculating_sum\n";
    boost::multiprecision::mpf_float::default_precision(1000);
    boost::multiprecision::mpf_float running_sum;
    std::fstream file_to_sum;
    file_to_sum.open(file_name, std::ios::in);
    std::string current_number_s;
    std::cout.precision(std::numeric_limits<boost::multiprecision::mpf_float>::digits10);

    while(getline(file_to_sum, current_number_s)){
        boost::multiprecision::mpf_float x(current_number_s);
        running_sum += x;
    }
    std::cout << running_sum << "\n";
    return running_sum;
}

void factorial_test(){
    std::cout << factorial(230) << std::endl;
}

int main()  { // NOLINT
//    factorial_test();
//    newton_euler_formula();
//    chudnovsky_formula();
//    leibniz_series(1, 10000000);
    chudnovsky_formula_no();
    std::cout << 1/(12*sum_from_file("series.txt")) << std::endl;
    return 0;
}
