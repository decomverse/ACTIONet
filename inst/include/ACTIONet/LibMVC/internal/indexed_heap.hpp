#ifndef INDEXED_HEAP_INCLUDED
#define INDEXED_HEAP_INCLUDED

#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

namespace libmvc {
namespace internal {

template <typename T, class Compare = std::less<T>>
class Indexed_Heap {
  static_assert(std::is_integral<T>::value, "Type has to be integral type");

 private:
  using Container = std::vector<T>;
  using Index = std::vector<T>;

 public:
  using container_type = Container;
  using value_compare = Compare;
  using value_type = typename Container::value_type;
  using size_type = typename Container::size_type;
  using reference = typename Container::reference;
  using const_reference = typename Container::const_reference;

 private:
  Container data;
  Index index;
  Compare comp;

 private:
  void element_swap(const size_type& a, const size_type& b) {
    std::swap(data[a], data[b]);
    index[data[a]] = a;
    index[data[b]] = b;
  }

  bool is_leaf(const size_type& pos) const {
    if ((pos >= data.size() / 2) && (pos < data.size())) return true;
    return false;
  }

  size_type left_child(const size_type& pos) const { return 2 * pos + 1; }

  size_type right_child(const size_type& pos) const { return 2 * pos + 2; }

  size_type parent(const size_type& pos) const { return (pos - 1) / 2; }

  void shift_down(size_type pos) {
    while (!is_leaf(pos)) {
      // std::cout << "shift _ down pos :" << pos << std::endl;
      auto min = pos;
      auto lc = left_child(pos);
      auto rc = right_child(pos);
      if ((lc < data.size()) && (!comp(data[lc], data[min]))) min = lc;
      if ((rc < data.size()) && (!comp(data[rc], data[min]))) min = rc;
      if (min == pos) return;
      element_swap(pos, min);
      pos = min;
    }
  }

 public:
  explicit Indexed_Heap(const Compare& compare = Compare()) : comp(compare) {}

  explicit Indexed_Heap(const size_type& size,
                        const Compare& compare = Compare())
      : index(size, std::numeric_limits<T>::max()), comp(compare) {
    data.reserve(size);
  }

 public:
  void resize(const size_type& size) {
    data.reserve(size);
    index.resize(size);
  }

  void clear() {
    data.clear();
    std::fill(index.begin(), index.end(), std::numeric_limits<T>::max());
  }

  const_reference top() const { return data[0]; }

  bool empty() const { return data.size() == 0; }

  size_type size() const { return data.size(); }

  void push(const value_type& value) {
    int curr = data.size();
    index[value] = curr;
    data.push_back(value);

    while (curr != 0 && !comp(data[curr], data[parent(curr)])) {
      element_swap(curr, parent(curr));
      curr = parent(curr);
    }
  }

  void pop() {
    index[0] = std::numeric_limits<T>::max();
    std::swap(data[0], data.back());
    data.pop_back();
    index[data[0]] = 0;
    if (data.size() != 0) shift_down(0);
  }

  /**
   * returns 1 if element is in heap, 0 otherwise
   */
  size_type count(const value_type& value) const {
    return index[value] != std::numeric_limits<T>::max();
  }

  /**
   * remove element at given position and restore heap
   */
  void erase(const size_type& pos) {
    if (pos == (data.size() - 1)) {
      index[data.back()] = std::numeric_limits<T>::max();
      data.pop_back();
    } else {
      element_swap(pos, data.size() - 1);
      index[data.back()] = std::numeric_limits<T>::max();
      data.pop_back();
      auto idx = pos;
      while ((idx != 0) && (!comp(data[idx], data[parent(idx)]))) {
        element_swap(idx, parent(idx));
        idx = parent(idx);
      }
      if (data.size() != 0) shift_down(idx);
    }
  }

  /**
   * return position in heap of given value
   */
  size_type operator[](const value_type& value) const { return index[value]; }
};

} // namespace internal
} // namespace libmvc

#endif

