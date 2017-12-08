/*


 */

#ifndef __SINGLEDATAITERATOR_H__
#define __SINGLEDATAITERATOR_H__

#include "unused.hxx"
#include <iostream>
#include <iterator>
#include <output.hxx>
#include <vector>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
int SDI_spread_work(int num_work, int thread, int max_thread);
#endif
  
#ifndef MAXREGIONBLOCKSIZE
#define MAXREGIONBLOCKSIZE 64
#endif

/// Region for the SingleDataIterator to iterate over
typedef std::vector<int> RegionIndices;

typedef std::pair<int, int> contiguousBlock;
typedef std::vector<contiguousBlock> contiguousBlocks;


#define BLOCK_REGION_LOOP_SERIAL(regionStr, index, ...)			\
  {const auto blocks = mesh->getRegion(regionStr).blocks;		\
  for(auto block = blocks.begin(); block < blocks.end() ; ++block){ \
  for(int index = (*block).first ; index <= (*block).second; ++index){ \
  __VA_ARGS__ \
    }}}

#define BLOCK_REGION_LOOP(regionStr, index, ...)			\
  {const auto blocks = mesh->getRegion(regionStr).blocks;		\
  BOUT_OMP(parallel for) \
  for(auto block = blocks.begin(); block < blocks.end() ; ++block){ \
  for(int index = (*block).first ; index <= (*block).second; ++index){ \
  __VA_ARGS__ \
    }}}

/// Helper function to create a RegionIndices, given the start and end
/// points in x, y, z, and the total y, z lengths
inline RegionIndices createRegionIndices(int xstart, int xend, int ystart, int yend,
                                         int zstart, int zend, int ny, int nz) {

  int len = (xend - xstart + 1) * (yend - ystart + 1) * (zend - zstart + 1);
  RegionIndices region(len);
  int j = 0;
  int x = xstart;
  int y = ystart;
  int z = zstart;

  bool done = false;
  j = -1;
  while (!done) {
    j++;
    region[j] = (x * ny + y) * nz + z;
    if (x == xend && y == yend && z == zend) {
      done = true;
    }
    ++z;
    if (z > zend) {
      z = zstart;
      ++y;
      if (y > yend) {
        y = ystart;
        ++x;
      }
    }
  }
  return region;
}

/*!
 * Set of indices - SingleDataIterator is dereferenced into these
 */
struct SIndices {
  SIndices(int i, int nx, int ny, int nz) : i(i), nx(nx), ny(ny), nz(nz) {}
  int i; // index of array
  int nx;
  int ny;
  int nz;

  // This is a little gross, but required so that dereferencing
  // SingleDataIterators works as expected
  SIndices *operator->() { return this; }
  const SIndices *operator->() const { return this; }

  /// Convert to a 2D index
  int i2d() const { return i / nz; }

  /// Convert to x, y, z indices
  int x() const { return (i / nz) / ny; }
  int y() const { return (i / nz) % ny; }
  int z() const { return (i % nz); }

  /*
   * Shortcuts for common offsets, one cell
   * in each direction.
   */

  /// The index one point +1 in x
  SIndices xp() const { return {i + ny * nz, nx, ny, nz}; }
  SIndices xpzp() const {
    return {(i + 1) % nz == 0 ? i + ny * nz - nz + 1 : i + ny * nz + 1, nx, ny, nz};
  }
  SIndices xpzm() const {
    return {i % nz == 0 ? i + ny * nz + nz - 1 : i + ny * nz - 1, nx, ny, nz};
  }
  /// The index one point -1 in x
  SIndices xm() const { return {i - ny * nz, nx, ny, nz}; }
  SIndices xmzp() const {
    return {(i + 1) % nz == 0 ? i - ny * nz - nz + 1 : i - ny * nz + 1, nx, ny, nz};
  }
  SIndices xmzm() const {
    return {i % nz == 0 ? i - ny * nz + nz - 1 : i - ny * nz - 1, nx, ny, nz};
  }
  /// The index one point +1 in y
  SIndices yp() const { return {i + nz, nx, ny, nz}; }
  SIndices ypp() const { return {i + 2 * nz, nx, ny, nz}; }
  /// The index one point -1 in y
  SIndices ym() const { return {i - nz, nx, ny, nz}; }
  SIndices ymm() const { return {i - 2 * nz, nx, ny, nz}; }
  /// The index one point +1 in z. Wraps around zend to zstart
  SIndices zp() const {
    return {(i + 1) % nz == 0 ? i - nz + 1 : i + 1, nx, ny, nz};
  }
  /// The index one point -1 in z. Wraps around zstart to zend
  SIndices zm() const { return {i % nz == 0 ? i + nz - 1 : i - 1, nx, ny, nz}; }

  /*!
   * Add an offset to the index for general stencils
   */
  SIndices offset(int dx, int dy, int dz) const {
    int z0 = i % nz;
    if (dz > 0) {
      int zp = z0;
      for (int j = 0; j < dz; ++j) {
        zp = (zp == nz - 1 ? 0 : zp + 1);
      }
      return {i + ny * nz * dx + nz * dy + zp - z0, nx, ny, nz};
    } else {
      int zm = z0;
      for (; dz != 0; ++dz) {
        zm = (zm == 0 ? nz - 1 : zm - 1);
      }
      return {i + ny * nz * dx + nz * dy + zm - z0, nx, ny, nz};
    }
  }
};


class Region {
public:
  RegionIndices indices;
  contiguousBlocks blocks;

  Region(int xstart, int xend, int ystart, int yend,
	 int zstart, int zend, int ny, int nz){
    indices = createRegionIndices(xstart, xend, ystart, yend, zstart, zend, ny, nz);
    blocks  = getContiguousBlocks();
  };
  Region(RegionIndices& indices) : indices(indices){
    blocks = getContiguousBlocks();
  }
  //Want to make this private to disable but think it may be needed as we put Regions into
  //maps which seems to need to be able to make "empty" objects.
  Region(){};

  /* 
   * Returns a vector of all contiguous blocks contained in the passed region.
   * Limits the maximum size of any contiguous block to maxBlockSize.
   * A contiguous block is described by the inclusive start and the exclusive end
   * of the contiguous block.
   */
  contiguousBlocks getContiguousBlocks() const{
    const int lastPoint = indices.size();
    contiguousBlocks result;
    int index = 0;
    
    while (index < lastPoint){
      const int startIndex = indices[index];
      int count = 1; //We will always have at least startPair in the block so count starts at 1
      
      //Consider if the next point should be added to this block
      for(index++; count<MAXREGIONBLOCKSIZE; index++){
	if((indices[index]-indices[index-1])==1){
	  count++;
	}else{//Reached the end of this block so break
	  break;
	}
      }
      
      //Add pair to output, denotes inclusive start and exclusive end
      result.push_back({startIndex, indices[index-1]});
    }
    
    return result;
  }
private:

};

/// Provides range-based iteration over indices.
/// If OpenMP is enabled, then this divides work between threads.
class SingleDataIterator {
private:
#ifdef _OPENMP
  /// This initialises OpenMP threads if enabled, and
  /// divides iteration index ranges between threads
  void omp_init();
#endif

  SingleDataIterator(); // Disable null constructor

  /// The internal start, end indices of region_iter, needed for OpenMP
  int icountstart, icountend;
  /// Vector of indices this SDI should iterate over
  std::vector<int>::const_iterator region_iter;

public:
  // iterator traits
  using difference_type = int;
  using value_type = SIndices;
  using pointer = const SIndices *;
  using reference = const SIndices &;
  using iterator_category = std::random_access_iterator_tag;

  // SIndexRange needs to be a friend so that it can modify the
  // private members to create the begin and end iterators
  friend class SIndexRange;
  /*!
   * Constructor. This sets index ranges.
   * If OpenMP is enabled, the index range is divided
   * between threads using the omp_init method.
   */
  SingleDataIterator(int nx, int ny, int nz, RegionIndices &region)
      : icountstart(0), icountend(region.size()), region_iter(region.cend()),
        nx(nx), ny(ny), nz(nz), region(region) {
#ifdef _OPENMP
    omp_init();
#endif
  }

  /// Size of the x, y, z dimensions
  const int nx, ny, nz;
  /// Set of indices to iterate over
  const RegionIndices &region;

  /// Pre-increment operator. Use this rather than post-increment when possible
  SingleDataIterator &operator++() {
    ++region_iter;
    return *this;
  }

  /// Post-increment operator
  SingleDataIterator operator++(int) {
    SingleDataIterator tmp(*this);
    ++(*this);
    return tmp;
  }

  /// Pre-decrement operator
  SingleDataIterator &operator--() {
    --region_iter;
    return *this;
  }

  /// Post-decrement operator
  SingleDataIterator operator--(int) {
    SingleDataIterator tmp(*this);
    --(*this);
    return tmp;
  }

  /*!
   * Dereference operators
   */
  SIndices operator*() const { return {*region_iter, nx, ny, nz}; }
  // This is a little gross, but required so that we don't need to do:
  //     (*iter).i
  // in order to access the elements of a SIndices, where iter is a
  // SingleDataIterator
  SIndices operator->() const { return {*region_iter, nx, ny, nz}; }

  /// Arithmetic operators
  SingleDataIterator &operator+=(int n) {
    this->region_iter += n;
    return *this;
  }

  SingleDataIterator &operator-=(int n) { return *this += -n; }

  /// Indexing operator
  SIndices operator[](int n) const { return {*(region_iter + n), nx, ny, nz}; }

  // The following three operators are friend functions so that they
  // have access to the private region_iter

  /// Equivalence operator
  friend bool operator==(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
    return lhs.region_iter == rhs.region_iter;
  }

  /// Less-than operator
  friend bool operator<(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
    return lhs.region_iter < rhs.region_iter;
  }

  /// Iterator difference operator
  friend difference_type operator-(const SingleDataIterator &lhs,
                                   const SingleDataIterator &rhs) {
    return lhs.region_iter - rhs.region_iter;
  }
};

/// Relational operators
inline bool operator!=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator==(lhs, rhs);
}

inline bool operator>(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return operator<(rhs, lhs);
}

inline bool operator<=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator>(lhs, rhs);
}

inline bool operator>=(const SingleDataIterator &lhs, const SingleDataIterator &rhs) {
  return !operator<(lhs, rhs);
}

/// Arithmetic operators with integers, needed for constant-time random access
inline SingleDataIterator operator+(SingleDataIterator lhs, int n) { return lhs += n; }
inline SingleDataIterator operator+(int n, SingleDataIterator rhs) { return rhs += n; }
inline SingleDataIterator operator-(SingleDataIterator lhs, int n) { return lhs -= n; }

#ifdef _OPENMP
inline int SDI_spread_work(int work, int cp, int np) {
  // Spread work between threads. If number of points do not
  // spread evenly between threads, put the remaining "rest"
  // points on the threads 0 ... rest-1.
  int pp = work / np;
  int rest = work % np;
  int result = pp * cp;
  if (rest > cp) {
    result += cp;
  } else {
    result += rest;
  }
  return result;
};

inline void SingleDataIterator::omp_init() {
  // In the case of OPENMP we need to calculate the range
  int threads = omp_get_num_threads();
  if (threads > 1) {
    int work = region.size();
    int current_thread = omp_get_thread_num();
    icountstart = SDI_spread_work(work, current_thread, threads);
    icountend = SDI_spread_work(work, current_thread + 1, threads);
  }
};
#endif

/*!
 * Specifies a vector of indices which can be iterated over
 * and begin() and end() methods for range-based for loops
 *
 * Example
 * -------
 *
 * Index ranges can be defined manually:
 *
 *     RegionIndices region {0, 2, 4, 8};
 *     SIndexRange r(2, 2, 2, region);
 *
 * then iterated over using begin() and end()
 *
 *     for (auto i = r.begin(); i < r.end(); i++ ) {
 *       output.write("%d,%d,%d\n", i.x, i.y, i.z);
 *     }
 *
 * or the more convenient range for loop:
 *
 *     for (auto i : r) {
 *       output.write("%d,%d,%d\n", i->x(), i->y(), i->z());
 *     }
 *
 * A common use for this class is to loop over
 * regions of a field:
 *
 *     Field3D f(0.0);
 *     for (auto i : f.region(REGION_NOBNDRY)) {
 *       f[i] = 1.0;
 *     }
 *
 * where REGION_NOBNDRY specifies a region not including
 * boundary/guard cells.
 */
class SIndexRange {
public:
  SIndexRange(int nx, int ny, int nz, Region &region)
      : nx(nx), ny(ny), nz(nz), region(region) {}

  int nx, ny, nz;
  Region &region;

  /*!
   * Resets DataIterator to the start of the range
   */
  SingleDataIterator begin() const {
    SingleDataIterator iter{nx, ny, nz, region.indices};
    iter.region_iter = region.indices.cbegin();
    std::advance(iter.region_iter, iter.icountstart);
    return iter;
  }

  /*!
   * Sets DataIterator to one index past the end of the range
   */
  SingleDataIterator end() const {
    SingleDataIterator iter{nx, ny, nz, region.indices};
    iter.region_iter = region.indices.cbegin();
    std::advance(iter.region_iter, iter.icountend);
    return iter;
  }

  /* 
   * Returns a vector of all contiguous blocks contained in the passed region.
   * Limits the maximum size of any contiguous block to maxBlockSize.
   * A contiguous block is described by the inclusive start and the exclusive end
   * of the contiguous block.
   */
  contiguousBlocks getContiguousBlocks() const{
    return region.blocks;
  }

};

#endif // __SINGLEDATAITERATOR_H__
