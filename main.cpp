#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <cmath>
#include <iterator>
#include <memory>
#include <algorithm>
#include <string>
#include <array>


typedef std::tuple<double, double, double> data_tuple;
typedef std::vector<data_tuple> frame;
typedef std::vector<frame> beads;


std::vector<data_tuple> coordinates_read (const std::string & name);

std::vector<data_tuple> centroid (beads & configs); // returns std::vector of centroid coordinates




int main() {

}


std::vector<data_tuple> centroid (beads & configs) {
    for (int i = 0; i < configs.size(); ++i) {
        for (int j = 0; j < configs[i].size(); ++i) {

        }
    }
}


namespace std {
    istream& operator >> (istream& in, data_tuple & coordinates) {
        double first, second, third;
        in >> first >> second >> third;
        coordinates = {first, second, third};
        return in;
    }

    ostream& operator << (ostream& out, const data_tuple & coordinates) {
        auto [first, second, third] = coordinates;
        out << first << ' ' << second << ' ' << third << ' ';
        return out;
    }
}


// Read data from columns in text-file.
std::vector<data_tuple> coordinates_read (const std::string & name) {
    std::ifstream fin(name);
    if (!fin.is_open()) throw std::runtime_error("Error opening file.");
    std::vector<data_tuple> tuples_vector;
    copy(std::istream_iterator<data_tuple> {fin},
         std::istream_iterator<data_tuple> {},
         back_inserter(tuples_vector));
    //copy(tuples_vector.begin(), tuples_vector.end(), std::ostream_iterator<data>(std::cout, "\n"));
    return tuples_vector;
}