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


const int frame_size = 698; // atoms
const int oxygen_search_list_length = 648;
const int T_O = 3; // oxygen appearance period in search list

const double group_step = 1.1; // angstroms hist step.


typedef std::tuple<double, double, double> data_tuple;
typedef std::vector<data_tuple> frame;
typedef std::vector<data_tuple> configuration;
typedef std::vector<frame> beads;
typedef std::vector<frame> frames;


std::vector<data_tuple> getting_centroid_coordinates (const std::string & coordinate_files_name, int & data_files_count);

std::vector<data_tuple> coordinates_read (const std::string & name);

std::vector<data_tuple> centroid (beads & configs); // returns std::vector of centroid coordinates

template<typename T>
T fromString (const std::string& s);

std::string exec (const std::string& str);


frames selection_oxygens (std::vector<data_tuple> & coordinates, const int & atoms_per_frame_count,
                          const int & search_area, const int & appearance_period);

std::vector<double> distances_in_frame (frame & atoms);



int main() {
    int data_files_count = fromString<int>(exec("find ./base_p-1.pos -type f | wc -l"));
    std::vector<data_tuple> centroid_coordinates = std::move(getting_centroid_coordinates("base_p-1.pos", data_files_count));
    frames oxygen_frames = std::move(selection_oxygens(centroid_coordinates, frame_size, oxygen_search_list_length, T_O));
    std::vector<std::vector<double>> distances_in_frames;

    for (auto & oxygen_frame : oxygen_frames) {
        distances_in_frames.emplace_back(distances_in_frame(oxygen_frame));
    }



}


// From all frames

/*
При написании этой функции файл содержал исключительно координаты 13М раствора солянки.
 Координаты атомов хлора были исключены из рассмотрения.
 698 - количество атомов в таком кадре
 646 - координаты атомов воды (дальше шли протоны)
 if ((j+(i)/698) % 3 == 0)
 Поскольку период не кратен трём, накапливается смещение при смещении по кадрам.
 В соответствии со стандартом, при делении целых чисел, дробная часть отбрасывается.
 Таким образом результат выражения в скобках есть не что иное, как целая часть при делении.
 Если полученное число кратно трём, то строка с указанным номером, содержащая координаты, отвечает атому кислорода.
 В общем-то, коль скоро именно кислороды нас и интересуют в этой ситуации, для 13М-раствора эту функцию можно оставить
 без изменений. На выходе можно получить вектор фреймов, или же одномерный массив координат с известном периодом.
*/
std::vector<data_tuple> find_oxygens (std::vector<data_tuple> & coordinates) {
    std::vector<data_tuple> oxygens;
    for (long i = 0; i < coordinates.size()-700; i += 698)
        for (long j = i; j < i + 646; ++j) //
            if ((j+(i)/698) % 3 == 0) // ! if ((j+(i)/698) % 3 == 0)
                oxygens.emplace_back(coordinates[j]);

    return oxygens;
}


frames selection_oxygens (std::vector<data_tuple> & coordinates, const int & atoms_per_frame_count,
                          const int & search_area, const int & appearance_period) {
    frame oxygens;
    frames oxygen_frames;
    for (long i = 0; i < coordinates.size() - atoms_per_frame_count; i += atoms_per_frame_count) {
        for (long j = i; j < i + search_area; ++j)
            if (((j + i) / atoms_per_frame_count) % appearance_period == 0) // ! unsafe conversion !
                oxygens.emplace_back(coordinates[j]);
        oxygens.clear();
        oxygen_frames.emplace_back(oxygens);
    }
    return oxygen_frames;
}


template<typename T, size_t... Is>
double distance_impl (T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return (std::sqrt((std::pow(std::get<Is>(t) - std::get<Is>(t1), 2) + ...)));
}

template <class Tuple>
double distance (const Tuple & t, const Tuple & t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return distance_impl(t, t1, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


std::vector<double> distances_in_frame (frame & atoms) {
    std::vector<double> results;
    for (int i = 0; i < atoms.size(); ++i)
        for (int j = 0; j < i; ++j)
           results.emplace_back(distance(atoms[i], atoms[j]));
    return results;
}


// Before start dist analyze, find maximum for scaling, otherwise use small scale. Of course, we need just two distances
// intervals, but we can do it more useful.
// Anyway it doesn't matter in this function.
std::vector<std::pair<int, int>> histogram (std::vector<double> & groups_borders, std::vector<double> & distances) {
    std::vector<std::pair<int, int>> result;
    for (int i = 0; i < distances.size(); ++i)
        for (int j = 0; j < groups_borders.size()-1; ++j)
            if (distances[i] >= groups_borders[j] && distances[i] <= groups_borders[j+1])
                result.emplace_back(j, i); // rewrite this shit
    return result;
}


std::vector <double> mesh (double left_border, const double & right_border, const double & step) {
    std::vector <double> xx ((right_border-left_border) / step);
    std::generate(xx.begin(), xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


frames atoms_per_frame (std::vector <data_tuple> & atoms, const int & T) { // T -- period for toms of one type in frame
    frames data;
    std::vector<data_tuple> frame;
    for (long j = 0; j < atoms.size(); j += T) {
        for (long i = j; i < j + T; ++i)
            frame.emplace_back(atoms[i]);
        data.emplace_back(frame);
        frame.clear();
    }
    return data;
}





//The function returns the terminal ans. Input - string for term.
std::string exec (const std::string& str) {
    const char* cmd = str.c_str();
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        result += buffer.data();
    result = result.substr(0, result.length()-1);
    return result;
}


// std::to_string not safe enough. It will be used everywhere instead of std::to_string.
template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


// Collects data in tuples for averaging. rename to centroid.
template<size_t Is = 0, typename... Tp>
void coordinates_average (std::tuple<Tp...>& coordinate, std::tuple<Tp...>& new_data, const double & data_count) {
    std::get<Is>(coordinate) += (std::get<Is>(new_data) / data_count);
    if constexpr(Is + 1 != sizeof...(Tp))
        coordinates_average<Is + 1>(coordinate, new_data, data_count);
}


std::vector<data_tuple> getting_centroid_coordinates (const std::string & coordinate_files_name, int & data_files_count) {
    std::vector<data_tuple> data = std::move(coordinates_read(coordinate_files_name + toString(0)));
    if (data_files_count == 1)
        return data;
    data_tuple zero_tuple = std::make_tuple(0.0, 0.0, 0.0);
    for (auto & i : data) i = zero_tuple;
    for (int i = 0; i < data_files_count; ++i)
        for (long j = 0; j < data.size(); ++j) {
            std::vector<data_tuple> new_data = std::move(coordinates_read(coordinate_files_name + toString(i)));
            coordinates_average(data[j], new_data[j], data_files_count);
        }
    return data;
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