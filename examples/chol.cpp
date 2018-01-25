// Chol
#include <iostream>
#include <cmath>
#include <random>

#include "hmat_cpp_interface.hpp"
#include "default_engine.hpp"

using namespace hmat;

template<typename T>
class TestAssemblyFunction : public SimpleAssemblyFunction<T> {
public:
  /// Point coordinates
  const DofCoordinates& points;

public:
  /** Constructor.

      \param _points Point cloud
      \param _k Wavenumber
   */
  TestAssemblyFunction(const DofCoordinates& _points)
    : SimpleAssemblyFunction<T>(), points(_points) {}
  typename Types<T>::dp interaction(int i, int j) const;
  double distanceTo(int i, int j) const;
  double diffTo(int i, int j) const;
  ScalarArray<T>* createRhs(double l) const;
};

template<typename T>
double TestAssemblyFunction<T>::diffTo(int i, int j) const
{
  double r = points.get(0, i) - points.get(0, j);
  return r;
}

template<typename T>
double TestAssemblyFunction<T>::distanceTo(int i, int j) const
{
  double r = sqrt((points.get(0, i) - points.get(0, j))*(points.get(0, i) - points.get(0, j)) );
  return r;
}

template<>
Types<S_t>::dp TestAssemblyFunction<S_t>::interaction(int i, int j) const {
  double distance = this->distanceTo(i, j) + 1.0;
  return  1.0 / distance;
}
template<>
Types<D_t>::dp TestAssemblyFunction<D_t>::interaction(int i, int j) const {
  double distance = this->distanceTo(i, j) + 1.0;
  return  1.0 / distance;
}
template<>
Types<C_t>::dp TestAssemblyFunction<C_t>::interaction(int i, int j) const {
  double distance = this->distanceTo(i, j) + 1.0;
  double diff = this->diffTo(i, j); // to have hermitian symmetry
  //double diff = std::abs(this->diffTo(i, j)); // to have symmetry but not hermitian
  Z_t result(cos(diff)/distance, sin(diff)/distance);
  return result;
}
template<>
Types<Z_t>::dp TestAssemblyFunction<Z_t>::interaction(int i, int j) const {
  double distance = this->distanceTo(i, j) + 1.0;
  double diff = this->diffTo(i, j); // to have hermitian symmetry
  //double diff = std::abs(this->diffTo(i, j)); // to have symmetry but not hermitian
  Z_t result(cos(diff)/distance, sin(diff)/distance);
  return result;
}

template<>
ScalarArray<S_t>* TestAssemblyFunction<S_t>::createRhs(double l) const {
  const int seed = 31;
  std::mt19937 generator(seed);
  std::normal_distribution<double> normal(0.0, l);

  const int n = (int) points.size();
  ScalarArray<S_t>* rhs = new ScalarArray<S_t>(n, 1);
  for (int i = 0; i < n; i++) {
    rhs->get(i, 0) = normal(generator);
  }
  return rhs;
}
template<>
ScalarArray<D_t>* TestAssemblyFunction<D_t>::createRhs(double l) const {
  const int seed = 31;
  std::mt19937 generator(seed);
  std::normal_distribution<double> normal(0.0, l);

  const int n = (int) points.size();
  ScalarArray<D_t>* rhs = new ScalarArray<D_t>(n, 1);
  for (int i = 0; i < n; i++) {
    rhs->get(i, 0) = normal(generator);
  }
  return rhs;
}
template<>
ScalarArray<C_t>* TestAssemblyFunction<C_t>::createRhs(double l) const {
  const int seed = 31;
  std::mt19937 generator(seed);
  std::normal_distribution<double> normal(0.0, l);

  const int n = (int) points.size();
  ScalarArray<C_t>* rhs = new ScalarArray<C_t>(n, 1);
  for (int i = 0; i < n; i++) {
    rhs->get(i, 0) = std::complex<float>(normal(generator),normal(generator));
  }
  return rhs;
}
template<>
ScalarArray<Z_t>* TestAssemblyFunction<Z_t>::createRhs(double l) const {
  const int seed = 31;
  std::mt19937 generator(seed);
  std::normal_distribution<double> normal(0.0, l);

  const int n = (int) points.size();
  ScalarArray<Z_t>* rhs = new ScalarArray<Z_t>(n, 1);
  for (int i = 0; i < n; i++) {
    rhs->get(i, 0) = std::complex<double>(normal(generator),normal(generator));
  }
  return rhs;
}

template<typename T, template<typename> class E> struct Configuration
{
    void configure(HMatInterface<T,E> &){}
};

hmat::StandardAdmissibilityCondition admissibilityCondition(3.);

template<typename T, template<typename> class E>
void go(const DofCoordinates& coord, double epsilon) {
  if (0 != HMatInterface<T, E>::init())
    return;
  {
    //bool use_symmetry = false;
    bool use_symmetry = true;
    ClusterTree* ct = createClusterTree(coord);
    std::cout << "ClusterTree node count = " << ct->nodesCount() << std::endl;
    TestAssemblyFunction<T> f(coord);
    HMatInterface<T, E>* hmat = nullptr;
    if (use_symmetry) {
      hmat = new HMatInterface<T, E>(ct, ct, kLowerSymmetric, &admissibilityCondition);
    } else {
      hmat = new HMatInterface<T, E>(ct, ct, kNotSymmetric, &admissibilityCondition);
    }
    std::cout << "HMatrix node count = " << hmat->nodesCount() << std::endl;
    Configuration<T, E>().configure(*hmat);
    if (use_symmetry) {
      hmat->assemble(f, kLowerSymmetric);
      hmat->factorize(hmat_factorization_chol);
    } else {
      hmat->assemble(f, kNotSymmetric);
      hmat->factorize(hmat_factorization_lu);
    }
    hmat_info_t info;
    hmat->info(info);
    std::cout << "Compression Ratio = "
              << 100 * ((double) info.compressed_size) / info.uncompressed_size
              << "%" << std::endl;
    hmat->createPostcriptFile("h_matrix.ps");

    std::cout << "Resolution...";
    double l = 0.01;
    ScalarArray<T>* rhs = f.createRhs(l);
    ScalarArray<T> rhsCopy(rhs->rows, 1);
    rhsCopy.copyMatrixAtOffset(rhs, 0, 0);
    hmat->solve(*rhs);
    std::cout << "done." << std::endl;

    std::cout << "Accuracy...";
    double rhsCopyNorm = rhsCopy.norm();
    int n = coord.size();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	rhsCopy.get(i, 0) -= static_cast<T>(f.interaction(i, j)) * rhs->get(j, 0);
      }
    }
    double diffNorm = rhsCopy.norm();
    std::cout << "Done" << std::endl;
    std::cout << "||Ax - b|| / ||b|| = " << diffNorm / rhsCopyNorm << std::endl;

    delete rhs;
  }
  HMatInterface<T, E>::finalize();
}

template<template<typename> class E>
void goA(char arithmetic, const DofCoordinates& coord, double epsilon) {
    switch (arithmetic) {
    case 'S':
        go<S_t, E>(coord, epsilon);
        break;
    case 'D':
        go<D_t, E>(coord, epsilon);
        break;
    case 'C':
        go<C_t, E>(coord, epsilon);
        break;
    case 'Z':
        go<Z_t, E>(coord, epsilon);
        break;
    default:
      std::cerr << "Unknown arithmetic code " << arithmetic << std::endl;
    }
}

int main(int argc, char **argv) {
  HMatSettings& settings     = HMatSettings::getInstance();
  settings.maxParallelLeaves = 10000;

  if (argc != 3 && argc != 4) {
    std::cout << "Usage: " << argv[0] << " n_points (S|D|C|Z) [epsilon]"
              << std::endl;
    return 0;
  }
  int n           = atoi(argv[1]);
  char arithmetic = argv[2][0];
  double epsilon  =  0.0;
  if(argc == 4 ) {
    epsilon = atof(argv[3]);
  }

  settings.compressionMethod = AcaPlus;
  if(epsilon != 0.0) {
    settings.assemblyEpsilon = epsilon;
    settings.recompressionEpsilon = epsilon;
  }
  settings.setParameters();
  settings.printSettings();
  double * x = new double[n];
  for(size_t i = 0; i < n; ++i)
  {
    x[i] = double(i);
  }
  DofCoordinates coord(x, 1, n, true);
  std::cout << "done.\n";

  goA<DefaultEngine>(arithmetic, coord, epsilon);
  return 0;
}
