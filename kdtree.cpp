/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <cmath>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
    if (first[curDim] < second[curDim]) return true;
    if (first[curDim] == second[curDim]) return first < second;
    return false;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
     double distCurBest = 0;
     double distPot = 0;

     for (int i = 0; i < Dim; i++) {
       distCurBest += ((target[i] - currentBest[i]) * (target[i] - currentBest[i]));
       distPot += ((target[i] - potential[i]) * (target[i] - potential[i]));
     }

     if (distPot < distCurBest) return true;
     if (distPot == distCurBest) return potential < currentBest;
     return false;
}

// helper function for quickselect (inspired by https://en.wikipedia.org/wiki/Quickselect)
// sorts list into two parts, those less than median value, and those more than median value
// places these two parts to the left and right of median value respectively
template <int Dim>
int KDTree<Dim>::Partition(vector<Point<Dim>>& list, int left, int right, int medianIndex, int dim) {

  double medianValue = list.at(medianIndex)[dim];

  // swap list[medianIndex] and list[end]
  // move median to end of list (i'm copying point values, not actually moving the points)
  Point<Dim> temp = list.at(medianIndex);
  list.at(medianIndex) = list.at(right);
  list.at(right) = temp;

  int storeIndex = left; // storeIndex holds what will become the final median index
  // pushes all points less than the median value to the left side of the list
  for (int i = left; i < right; i++) {
    if (list.at(i)[dim] < medianValue || (list.at(i)[dim] == medianValue && list[i] < temp)) {
      // swap list[storeIndex] and list[i]
      Point<Dim> temp2 = list.at(storeIndex);
      list.at(storeIndex) = list.at(i);
      list.at(i) = temp2;
      storeIndex++;
    }
  }

  // swap list[end] and list[storeIndex] so that median point is now in the median index
  Point<Dim> temp3 = list.at(right);
  list.at(right) = list.at(storeIndex);
  list.at(storeIndex) = temp3;

  return storeIndex;
}

// helper function for ctor
// finds the k-th smallest point (median) in resepct to dimension dim in the range [left, right]
// also partitions list
template <int Dim>
Point<Dim> KDTree<Dim>::quickSelect(vector<Point<Dim>>& list, int left, int right, int k, int dim) { // should i return Point&?
  if (left == right) return list.at(left); // if list only contains one element, return that element
  int pivotIndex = ((left + right) / 2); // selects a random pivotIndex between start and end
  pivotIndex = Partition(list, left, right, pivotIndex, dim);

  if (k == pivotIndex) return list.at(k);
  else if (k < pivotIndex) return quickSelect(list, left, pivotIndex - 1, k, dim);
  else return quickSelect(list, pivotIndex + 1, right, k, dim);
}

// recursive helper for constructor
// builds and organizes tree
template <int Dim>
void KDTree<Dim>::constructorHelper(vector<Point<Dim>>& list, KDTreeNode*& subroot, int left, int right, int curDim) {
    //base case
    if (left > right) return;

    int medianIndex = (right + left) / 2; // floor?

    // put median of current dimension in "middle" of sublist
    subroot = new KDTreeNode(quickSelect(list, left, right, medianIndex, curDim));

    // build left subtree
    constructorHelper(list, subroot->left, left, medianIndex - 1, (curDim + 1) % Dim);

    // build right subtree
    constructorHelper(list, subroot->right, medianIndex + 1, right, (curDim + 1) % Dim);
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */

     // sortedPoints will be semi-sorted after quickSelect call
     root = NULL;
     vector<Point<Dim>> sortedPoints; // doesn't need to be private member since the data is stored in nodes

     // copy newPoints into a vector
     for (unsigned long i = 0; i < newPoints.size(); i++) {
       sortedPoints.push_back(newPoints.at(i));
     }
     constructorHelper(sortedPoints, root, 0, sortedPoints.size() - 1, 0);
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree& other) {
  /**
   * @todo Implement this function!
   */
   root = NULL;
   copyConstructorHelper(root, other.root);
   size = other.size;
}


template <int Dim>
void KDTree<Dim>::copyConstructorHelper(KDTreeNode*& copy, KDTreeNode* original) {
    if (original == NULL) return;

    copy = new KDTreeNode(original->point);
    copyConstructorHelper(copy->left, original->left);
    copyConstructorHelper(copy->right, original->right);
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree& rhs) {
  /**
   * @todo Implement this function!
   */
  destructorHelper(this->root);
  copyConstructorHelper(root, rhs.root);
  return *this;
}


template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
   destructorHelper(root);
}

template <int Dim>
void KDTree<Dim>::destructorHelper(KDTreeNode*& subroot) {
  if (subroot == NULL) return;
  destructorHelper(subroot->left);
  destructorHelper(subroot->right);
  delete subroot;
  subroot = NULL;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    Point<Dim> ret = root->point;
    double rad = -1;
    //nearestNeighborHelper(query, root, ret, 0, rad);
    nearestNeighborHelper(query, root, ret, 0, rad);
    return ret;
}

template <int Dim>
void KDTree<Dim>::nearestNeighborHelper(const Point<Dim>& target, KDTreeNode* subroot, Point<Dim>& closest, int curDim, double& radius) const
{
  // base case if we reached a leaf
  if (subroot == NULL) return;

  // recursive call to find nearest neighbor
  if (target[curDim] < subroot->point[curDim] || (target[curDim] == subroot->point[curDim] && target < subroot->point))
    // if target belongs in left subtree based on the current dimension
    nearestNeighborHelper(target, subroot->left, closest, (curDim + 1) % Dim, radius);
  else
    // if target belongs in right subtree based on the current dimension (target[dim] > )
    nearestNeighborHelper(target, subroot->right, closest, (curDim + 1) % Dim, radius);

    // these lines should be able to do everything, but they fail three tests. Probably need to check other side of tree, but below implementation is easier
//  if (shouldReplace(target, closest, subroot->point))
//    closest = subroot->point;

  // find radius between target and subroot
  double dist = 0;
  for (int i = 0; i < Dim; i++) {
    dist += (target[i] - subroot->point[i]) * (target[i] - subroot->point[i]);
  }
  dist = sqrt(dist);

  // update radius if subroot is closer than previous closer
  if (radius == -1 || (dist == radius && subroot->point < closest) || dist < radius) {
    radius = dist;
    closest = subroot->point;
  }

  // check to see if the current splitting plane's distance from search node is within current radius. If so, must search opposite subtree
  if (radius >= abs(subroot->point[curDim] - target[curDim])) {
    if (target[curDim] < subroot->point[curDim] || (target[curDim] == subroot->point[curDim] && target < subroot->point))
      nearestNeighborHelper(target, subroot->right, closest, (curDim + 1) % Dim, radius);
    else
      nearestNeighborHelper(target, subroot->left, closest, (curDim + 1) % Dim, radius);
  }
  
}
