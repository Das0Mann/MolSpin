#include <armadillo>
#include <iostream>

class Test {
private:
    arma::vec tdFreqs;

public:
    Test() : tdFreqs(3, arma::fill::zeros) {
        std::cout << "Constructor: tdFreqs initialized" << std::endl;
    }

    Test(const Test& other) : tdFreqs(other.tdFreqs) {
        std::cout << "Copy constructor: tdFreqs copied" << std::endl;
    }

    Test& operator=(const Test& other) {
        if (this != &other) {
            tdFreqs = other.tdFreqs;
            std::cout << "Assignment operator: tdFreqs assigned" << std::endl;
        }
        return *this;
    }

    ~Test() {
        std::cout << "Destructor: tdFreqs cleaned up" << std::endl;
    }
};

int main() {
    Test t1;
    Test t2 = t1;
    t2 = t1;
    return 0;
}
